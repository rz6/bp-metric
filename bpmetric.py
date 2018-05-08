import intervaltree as itr
import pandas as pd
import numpy as np
import os
import warnings
from math import log
from argparse import ArgumentParser


###################################

def checkColumns(df):
    assert {'start','end'} & set(df.columns.tolist()) == {'start','end'} , "TADs file must have start, end columns!"
    assert all(df['start'] >= 0) , "TAD start must be >= 0!"
    assert all(df['end'] > 0) , "TAD end must be > 0!"

def parseTads(path,sep='\t'):
    df = pd.read_csv(path,sep=sep)
    checkColumns(df)
    if 'chromosome' not in df.columns:
        df['chromosome'] = pd.Series([1] * len(df))
    df.sort_values(['chromosome','start','end'],inplace=True)
    if 'id' not in df.columns:
        df['id'] = df.index + 1
    tads = {}
    for c in set(df.chromosome):
        itree = itr.IntervalTree()
        cdf = df.loc[df.chromosome == c,]
        for idx, row in cdf.iterrows():
            itree.addi(row['start'],row['end'],row['id'])
        # check if file do not contain overlaps
        with_overlaps = itree.copy()
        with_overlaps.split_overlaps()
        if len(with_overlaps) != len(itree):
            print "Found overlapping tads on chromosome {0} in {1} file:".format(c,path)
            for t in itree - with_overlaps:
                print t
            print "Discarding overlaps ..."
            with_overlaps.merge_equals()
            itree = with_overlaps
        tads[c] = itree
    # if there were any overlaps new id need to be assigned
    counter = max(map(lambda x: max(x,key=lambda y: y.data).data,tads.values())) + 1
    renumerated = {}
    for c,v in tads.iteritems():
        rt = itr.IntervalTree()
        for t in v:
            if t.data:
                rt.add(t)
            else:
                rt.addi(t.begin,t.end,counter)
                counter += 1
        renumerated[c] = rt
    return renumerated

def fillGaps(tads,N,counter):
    # THIS FUNCTION WORKS FOR SINGLE CHROMOSOME
    # chromosome is [0,N] interval
    gaps = itr.IntervalTree([itr.Interval(0,N,1)])
    for tad in tads:
        gaps.chop(tad.begin,tad.end)
    # create new tads object with filled gaps
    filled = tads.copy()
    for gap in gaps:
        filled.addi(gap.begin,gap.end,counter)
        counter += 1
    return counter, filled

def align(tads1,tads2):
    atads1 = {}
    atads2 = {}
    c1 = max(map(lambda x: max(x,key=lambda y: y.data),tads1.values())).data + 2
    c2 = max(map(lambda x: max(x,key=lambda y: y.data),tads2.values())).data + 2
    for c in set(tads1) & set(tads2):
        t1 = tads1[c]
        t2 = tads2[c]
        N = max(max(t,key=lambda x: x.end).end for t in [t1,t2])
        c1, ft1 = fillGaps(t1,N,c1)
        c2, ft2 = fillGaps(t2,N,c2)
        atads1[c] = ft1
        atads2[c] = ft2
    return atads1, atads2

def removeUnmappable(tads1,tads2):
    # unmappable regions are marked with negative ids
    u1 = itr.IntervalTree([t for t in tads1 if t.data < 0])
    u2 = itr.IntervalTree([t for t in tads2 if t.data < 0])
    u = u1 | u2
    u.split_overlaps()
    u.merge_equals()
    for gap in u:
        tads1.chop(gap.begin,gap.end)
        tads2.chop(gap.begin,gap.end)
    return tads1, tads2

def removeUnmappableAll(tads1,tads2):
    atads1 = {}
    atads2 = {}
    for c in set(tads1) & set(tads2):
        t1, t2 = removeUnmappable(tads1[c],tads2[c])
        atads1[c] = t1
        atads2[c] = t2
    return atads1, atads2

def it2csv(tads):
    rows = []
    for c, it in tads.iteritems():
        for s,e,i in it:
            rows.append([i,c,s,e])
    return pd.DataFrame(rows,columns=['id','chromosome','start','end']).sort_values(['chromosome','start'])

def overlaps(tads1,tads2):
    o = []
    for c in set(tads1) & set(tads2):
        for s1,e1,id1 in tads1[c]:
            for s2,e2,id2 in tads2[c][s1:e1]:
                o.append([c,s1,e1,id1,s2,e2,id2])
    return pd.DataFrame(o,columns=['chromosome','start1','end1','id1','start2','end2','id2'])

def JI(tads1,tads2):
    jiscore = {}
    for c in set(tads1) & set(tads2):
        boundaries1 = set([t.begin for t in tads1[c]]) | set([t.end for t in tads1[c]])
        boundaries2 = set([t.begin for t in tads2[c]]) | set([t.end for t in tads2[c]])
        jiscore[c] = 1 - 1.0 * len(boundaries1 & boundaries2) / len(boundaries1 | boundaries2)
    return jiscore

def overlap(s1,s2,e1,e2):
    return 1.0 * min(e1,e2) - max(s1,s2)

def BP(tads1,tads2,saveLocalBPscore=''):
    bpscore = {}
    o = overlaps(tads1,tads2)
    o['localScore'] = o.apply(lambda x: overlap(x.start1,x.start2,x.end1,x.end2)/max(x.end1-x.start1,x.end2-x.start2),axis=1)
    if saveLocalBPscore:
        lbp = o.copy()
        lbp.localScore = 1.0 - lbp.localScore
        p,ext = os.path.splitext(saveLocalBPscore)
        lbp.sort_values(['chromosome','start1','start2']).to_csv('{0}-localBP{1}'.format(p,ext),sep='\t',index=False)
    for c in set(o.chromosome):
        oc = o.loc[o.chromosome == c,]
        ostarts = oc[['start1','start2']].max(axis=1)
        oends = oc[['end1','end2']].min(axis=1)
        L = (oends - ostarts).sum()
        bpscore[c] = 1.0 - (oc.localScore * (oc.apply(lambda x: 1.0*min(x.end1,x.end2)-max(x.start1,x.start2),axis=1) / L)).sum()
    return bpscore

def localMI(tads1,tads2,L):
    def lMI(y,l):
        ovl = overlap(y.start1,y.start2,y.end1,y.end2)
        return -1.*(ovl / l) * log(ovl * l / ((y.end1-y.start1) * (y.end2-y.start2)),2)
    return overlaps(tads1,tads2) \
            .groupby('chromosome') \
            .apply(lambda x: pd.concat(
                [x,pd.Series(x.apply(lambda y: lMI(y,L[x.name]),axis=1), name='localScore')],
                axis=1))

def VI(tads1,tads2,saveLocalMIscore=''):
    L = {}
    H = {}
    for c in set(tads1) & set(tads2):
        l1 = np.array([t.end-t.begin for t in tads1[c]])
        l2 = np.array([t.end-t.begin for t in tads2[c]])
        assert l1.sum() == l2.sum()
        L[c] = 1.0 * l1.sum()
        p1 = l1 / L[c]
        p2 = l2 / L[c]
        H[c] = [(-1.0 * p1 * np.log2(p1)).sum(), (-1.0 * p2 * np.log2(p2)).sum()]
    lmi = localMI(tads1,tads2,L)
    if saveLocalMIscore:
        p,ext = os.path.splitext(saveLocalMIscore)
        lmi.sort_values(['chromosome','start1','start2']).to_csv('{0}-localMI{1}'.format(p,ext),sep='\t',index=False)
    return {k : sum(H[k]) + 2 * v.localScore.sum() for k,v in lmi.groupby('chromosome')}

###################################

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument("--tads1",dest="t1path",required=True,
            help="Path to file containing first set of TADs")
    parser.add_argument("--tads2",dest="t2path",required=True,
            help="Path to file containing second set of TADs")
    parser.add_argument("-m","--metric",dest="metric",choices=[1,2,3],required=True,type=int,
            help="Which distance to use: 1 - JI, 2 - BP, 3 - VI")
    parser.add_argument("-u","--unmappable",dest="removeu",type=int,choices=[0,1],default=0,
            help="Whether to remove unmappable regions.")
    parser.add_argument("-p","--preprocessed",dest="savep",type=int,choices=[0,1],default=0,
            help="Whether to save preprocessed file, i.e. with filled gaps and removed unmappable regions.")
    parser.add_argument("-l","--local_score",dest="savel",type=int,choices=[0,1],default=0,
            help="Whether to save file with local score.")
    parser.add_argument("-s","--separator",dest="sep",default='\t',
            help="Column separator for TAD files.")
    args = parser.parse_args()
    # args:
    # tads1 in tab separated file (required)
    # tads2 in tab separated file (required)
    # remove unmappable (optional) --> unmappable tads should be marked by user with negative id
    # return parsed (optional) --> save preprocessed domains to file
    # return matching (optional) --> return file with contribution score

    # tads is csv file with following columns:
    # [id] [chromosome] start end
    # if there are 2 consecutive TADs, say with ids 5 and 6, which share boundary
    # then tad5.end = tad6.start

    # - domains which belong to the same chromosome in given file can't overlap
    # - domains don't have to be sorted --> they will be sorted during preprocessing
    # - gaps between any pair of neighboring domains will be filled with artificial domain
    #   this artificial tads will start their id from last tad id +1, so in order to
    #   find artificial tads all you need to do is sort tad ids ascending and find a tad
    #   s.t. its id is equal to previous tad id + 1 --> this will be first artificial tad

    # test
    t1 = parseTads(args.t1path,sep=args.sep)
    t2 = parseTads(args.t2path,sep=args.sep)

    name1 = os.path.splitext(os.path.split(args.t1path)[1])[0]
    name2 = os.path.splitext(os.path.split(args.t2path)[1])[0]

    tads1,tads2 = align(t1,t2)
    if args.removeu:
        tads1,tads2 = removeUnmappableAll(tads1,tads2)
    if args.savep:
        it2csv(tads1).to_csv('{0}-preprocessed.csv'.format(name1),sep='\t',index=False)
        it2csv(tads2).to_csv('{0}-preprocessed.csv'.format(name2),sep='\t',index=False)

    if args.metric == 1:
        if args.savel:
            warnings.warn("Specified flag savel, but for selected metric (JI) there is no local score! Omitting ...")
        result = JI(tads1,tads2)
    elif args.metric == 2:
        if args.savel:
            result = BP(tads1,tads2,saveLocalBPscore='{0}-{1}-localBP.csv'.format(name1,name2))
        else:
            result = BP(tads1,tads2)
    elif args.metric == 3:
        if args.savel:
            result = VI(tads1,tads2,saveLocalMIscore='{0}-{1}-localMI.csv'.format(name1,name2))
        else:
            result = VI(tads1,tads2)

    # print result to std out
    print "chromosome\tdistance"
    for c, val in result.iteritems():
        print "{0}\t{1}".format(c,val)
