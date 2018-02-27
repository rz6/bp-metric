import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import get_cmap
from math import ceil
from argparse import ArgumentParser
from mpl_toolkits.axes_grid1 import make_axes_locatable
import config as cnf


def segmentation(path,sep='\t',binSize=1):
    df = pd.read_csv(path,sep=sep)
    for column in ['start1','end1','start2','end2']:
        df[column] /= binSize
    return df.sort_values(['chromosome','start1','start2'])

def segments(sgn):
    ss = sgn.apply(lambda x: pd.Series({'start' : max(x.start1,x.start2), 'end' : min(x.end1,x.end2)}),axis=1)
    ssWithInt = pd.concat([sgn.chromosome, ss[['start','end']], sgn.localScore],axis=1)
    # fill gaps with nan
    ssContinous = []
    for chromosome, df in ssWithInt.groupby('chromosome'):
        if df.iloc[0].start > 0:
            ssContinous.append([chromosome,0,df.iloc[0].start,np.nan])
        ssContinous.append(df.iloc[0].tolist())
        for i in range(1,len(df.index)):
            if df.iloc[i].start - df.iloc[i-1].end > 0:
                ssContinous.append([chromosome, df.iloc[i-1].end, df.iloc[i].start, np.nan])
            ssContinous.append(df.iloc[i].tolist())
    return pd.DataFrame(ssContinous, columns=['chromosome','start','end','localScore'])

def segments2mtx(seg):
    return {c : sum([[row.localScore] * int((row.end - row.start)) for idx, row in lbp.iterrows()],[]) for c,lbp in seg.groupby('chromosome')}

def cut(tends,bs,be):
    return (tends[(bs <= tends) & (tends <= be)] - bs).tolist()

def plot_segmentation(ax,data,**kwargs):
    cmap = get_cmap(kwargs['heatmap_color'])
    cmap.set_bad(color=kwargs['um_color'])
    if kwargs['data_range'] is None:
        im = ax.imshow(data,
                cmap=cmap,
                interpolation='none',
                extent=[0,data.shape[1],0,10])
    else:
        im = ax.imshow(data,
                vmin=kwargs['data_range'][0],
                vmax=kwargs['data_range'][1],
                cmap=cmap,
                interpolation='none',
                extent=[0,data.shape[1],0,10])
    ax.autoscale(False)

    tads_ends1 = kwargs['tads_ends1']
    tads_ends2 = kwargs['tads_ends2']

    all_tad_ends1 = set(tads_ends1) & set(tads_ends2)
    all_tad_ends2 = set(tads_ends1) ^ set(tads_ends2)
    for x in all_tad_ends1:
        ax.axvline(x=x,ymin=0,ymax=10,color=kwargs['ovb_color'],linestyle='--')
    for x in all_tad_ends2:
        ax.axvline(x=x,ymin=0,ymax=10,color=kwargs['novb_color'],linestyle='--')
    if tads_ends1:
        ax.plot(tads_ends1, [0] * len(tads_ends1), color=kwargs['ovb_color'], marker='^', ms=10)
    if tads_ends2:
        ax.plot(tads_ends2, [10] * len(tads_ends2), color=kwargs['ovb_color'], marker='v', ms=10)

    ax.yaxis.set_visible(False)
    if kwargs['labels']:
        ll,lc = map(np.array,zip(*kwargs['labels']))
        ax.set_xticks(lc,minor=False)
        ax.set_xticklabels(ll,rotation = 45,fontsize=10)
        ax.tick_params(axis='x',direction='out',which='major',length=5,top='on',bottom='off',labeltop='on',labelbottom='off')
    else:
        ax.tick_params(axis='x',which='major',top='off',bottom='off',labeltop='off',labelbottom='off')

    ax.set_xticks(np.arange(data.shape[1]) + 1,minor=True)
    ax.set_xticklabels(kwargs['xcoords'], fontsize=10, fontname='monospace', weight='bold', minor=True)
    ax.tick_params(axis='x',which='minor',length=5,direction="out",top='off',bottom='on')

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if kwargs['put_label']:
        ax.set_xlabel("$i$",fontsize=20)
    ax.xaxis.set_label_position('top')
    if 'title' in kwargs:
        ax.set_title(kwargs['title'] + '\n' * 4)
    return im

##########################################################################################

def guard():
    if not cnf.HEATMAP_COLOR in plt.colormaps():
        raise Exception("Given heatmap color not found in pyplot cmaps! Available cmaps are: " , plt.colormaps())
    if not all(isinstance(c,basestring) for c in [cnf.UNMAPPABLE_COLOR,cnf.OVERLAPPING_BOUNDARIES_COLOR,cnf.NONOVERLAPPING_BOUNDARIES_COLOR]):
        raise Exception("Colors (unmappable and both boundaries) must be of type string!")
    if not all(isinstance(i,int) for i in [cnf.BIN_LABEL_FREQ,cnf.ROW_SIZE]):
        raise Exception("bin lable frequency and row size must be integer!")
    if not all(isinstance(i,float) or isinstance(i,int) for i in [cnf.FIG_HEIGHT,cnf.FIG_WIDTH]):
        raise Exception("figure width and height must be float or int!")
    if not all(i > 0 for i in [cnf.BIN_LABEL_FREQ,cnf.FIG_HEIGHT,cnf.FIG_WIDTH,cnf.ROW_SIZE]):
        raise Exception("bin lable frequency, width, height and row size must larger than 0!")

if __name__ == "__main__":
    # parse parameters
    parser = ArgumentParser()
    parser.add_argument("-s","--segmentation",dest="infile",required=True,
            help="File containing segmentation to be plotted.")
    parser.add_argument("-r","--resolution",dest="resolution",type=int,default=1,
            help="Resolution of Hi-C map, i.e. number of bp in 1 bin.")
    parser.add_argument("-o","--out",dest="outfile",default='segmentation.pdf',
            help="File to save resulting segmentation plot.")
    parser.add_argument("-c","--chromosomes",nargs='+',type=int,dest="chromosomes",default=[],
            help="Chromosomes to plot, if not specified all chromosomes will be plotted.")
    args = parser.parse_args()
    guard()

    segs = segmentation(args.infile,binSize=args.resolution)
    properties = {}

    with PdfPages(args.outfile) as pdf:
        for chromosome, cs in segs.groupby('chromosome'):
            if args.chromosomes and not (chromosome in args.chromosomes):
                continue
            csegments = segments(cs)
            labels = []
            counter = 1
            for _, row in csegments.reset_index().iterrows():
                if np.isnan(row.localScore):
                    continue
                labels.append((counter,row.start + (row.end-row.start)/2.))
                counter += 1
            nbins = sum(csegments.end - csegments.start)
            nlines = int(ceil(nbins / cnf.ROW_SIZE))
            mtx = segments2mtx(csegments)[chromosome]
            if cnf.BIN_LABEL_FREQ is not None:
                xcoords = [str(i) if i % cnf.BIN_LABEL_FREQ == 0 else '' for i in range(1,len(mtx))]
            else:
                xcoords = []
            fig = plt.figure()
            for i in range(nlines):
                submtx = mtx[i * cnf.ROW_SIZE : (i+1) * cnf.ROW_SIZE]
                ax = plt.subplot2grid((nlines+1, cnf.ROW_SIZE), (i, 0), colspan=len(submtx))
                # subset must be sorted
                tends1 = cut(cs.end1.drop_duplicates().values, i*cnf.ROW_SIZE, (i+1)*cnf.ROW_SIZE)
                tends2 = cut(cs.end2.drop_duplicates().values, i*cnf.ROW_SIZE, (i+1)*cnf.ROW_SIZE)
                if i == 0:
                    if 0 not in tends1:
                        tends1 = [0] + tends1
                    if 0 not in tends2:
                        tends2 = [0] + tends2
                sublabels = [(l,c-(i*cnf.ROW_SIZE)) for l,c in labels if i*cnf.ROW_SIZE<= c <= (i+1)*cnf.ROW_SIZE]
                properties = {
                        'tads_ends1' : tends1,
                        'tads_ends2' : tends2,
                        'labels' : sublabels,
                        'put_label' : False if i else True,
                        'xcoords' : xcoords[i*cnf.ROW_SIZE:(i+1)*cnf.ROW_SIZE],
                        'data_range' : cnf.DATA_RANGE,
                        'heatmap_color' : cnf.HEATMAP_COLOR,
                        'ovb_color' : cnf.OVERLAPPING_BOUNDARIES_COLOR,
                        'novb_color' : cnf.NONOVERLAPPING_BOUNDARIES_COLOR,
                        'um_color' : cnf.UNMAPPABLE_COLOR
                    }
                if i == 0:
                    properties['title'] = 'segmentation of chromosome {0}'.format(chromosome)
                im = plot_segmentation(ax, np.array([submtx]), **properties)

            ax = plt.subplot2grid((nlines+1, cnf.ROW_SIZE), (nlines, 0), colspan=cnf.ROW_SIZE)
            cb = fig.colorbar(im, cax=ax, orientation='horizontal')
            #cb.set_label(r"$d_{A,B}^{\mathbf{BP}}(i) / 0.5$",size=15)
            cb.set_label(r"$d_{A,B}(i)$",size=15)
            ax.set_aspect(0.05)

            fig.set_size_inches(cnf.FIG_WIDTH, (nlines+1) * cnf.FIG_HEIGHT)
            pdf.savefig(fig,bbox_inches='tight')
