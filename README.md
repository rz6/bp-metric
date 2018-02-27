# BP metric #
-------------

## Quick summary ##
-------------------

This is a python based package for comparing chromosome segmentations, in particular those produced during calling of Topologically Associating Domains (TADs) from Hi-C data. Given 2 partitions of chromosome(s) (i.e. 2 sets of TADs) this package enables calculation of 3 different distances:

* Jaccard Index (JI) distance,
* BP distance,
* Variation of Information (VI) distance.

In case of BP and VI distances it is also possible to keep local score and illustrate it for deeper insight into similarity of subchromosomal regions. The details of calculation and meaning of all mentioned distances is precisely disscussed in a paper. Here only information about scripts usage can be found.

## Prerequisites ##
-------------------

To run this package following dependencies are required (versions for which tests were performed are given in parenthesis):

* python (2.7.12)
* numpy (1.11.0)
* pandas (0.17.1)
* matplotlib (1.5.1)
* intervaltree (2.0)

if you are not sure whether your system satisfies all dependencies run checkDep.py script:

```sh
python checkDep.py
```

## Configuration ##
-------------------

Just clone this repo.

## Usage ##
-----------

There are 2 scripts to use:

* bpmetric.py - calculation of distance between 2 chromosome partitions.
* plotSegmentation.py - illustrate partitions (dis)similarity.

### Distance calculation ###
----------------------------

Mandatory flags are:

* --tads1, --tads2 --> path to tab separated files with a list of TADs (partitions) from 2 experiments
* -m --> metric to calculate: 1 - JI, 2 - BP, 3 - VI

Both TADs files must contain start and end columns. There must be no (self)overlaps in both TAD files. Additional columns may also be specified:

* chromosome
* id

If user want to discard some regions (for example unmappable regions) from comparison then this regions must have negative id. In such case script should be executed with -u flag turned on (i.e. equal to 1). Both TAD sets will be parsed:

* negative id regions will be cut out from TAD set,
* both sets will be aligned (i.e. some regions in one set will be trimmed) to produce matching intervals.

Additionally user may specify -p flag to save those parsed sets to file for later analysis. Finally user can keep local scores by specifying -l flag (this is only available for BP and VI distances). An example of input files containing TADs can be found in sample-data directory:

* sample-data/complete-tads --> with complete set of TADs produced using Armatus software (https://github.com/kingsfordgroup/armatus)
* sample-data/trimmed-tads --> reduced set of TADs

An example command to produce distance between 2 sample TAD sets is shown below:

```sh
python bpmetric.py --tads1 sample-data/complete-tads/IPF-hESC-all-HindIII-40k-TADs-Armatus --tads2 sample-data/complete-tads/IPF-hMSC-all-HindIII-40k-TADs-Armatus -m 2 -u 1 -l 1 > distance-ESC-vs-MSC
```

The above command will keep local BP scores, which can be later used to illustrate local segmentation. The global score is printed to std out in 2 column format: chromosome	distance, so remember to redirect output to file.

### (Dis)similarity illustration ###
------------------------------------

To plot local scores you need to specify the input file and typically you will also want to specify resolution (as long as your TADs didn't have their coordinates in bins):

```sh
python plotSegmentation.py -s sample-data/results/IPF-hESC-all-HindIII-40k-TADs-Armatus-IPF-hMSC-all-HindIII-40k-TADs-Armatus-localBP-localBP.csv -r 40000 -c 1 2
```

-c flag is optional and indicates, which chromosomes should be plotted. If not specified all of them will be illustrated.

There are additional options in config.py file to configure plotting area:

* ROW_SIZE --> a maximum number of bins, which can fit in single row,
* FIG_WIDTH --> page width in inches,
* FIG_HEIGHT --> height of single row in inches, so page height is FIG_HEIGHT times number of rows,
* BIN_LABEL_FREQ --> how often (every BIN_LABEL_FREQ) to put chromosome coordinate (bin number) in bottom x axis,
* DATA_RANGE --> range for colorscale, note that for BP score it is better to keep this as [0,1], however for VI score it should be None in general,
* HEATMAP_COLOR --> a pallette for segments coloring according to their local score. For BP local score it is better to set this as symmetric pallette whereas for VI local score a pallette with larger number of sequential color should be picked,
* UNMAPPABLE_COLOR --> color for unmappable segments,
* OVERLAPPING_BOUNDARIES_COLOR --> color to mark overlapping boundaries,
* NONOVERLAPPING_BOUNDARIES_COLOR --> color to mark nonoverlapping boundaries.

An example of (dis)similarity illustration can be found in: sample-data/results/segmentation.pdf, which is a illustration of IPF-hESC-all-HindIII-40k-TADs-Armatus-IPF-hMSC-all-HindIII-40k-TADs-Armatus-localBP-localBP.csv of chromosomes 1 and 2.

Briefly speaking this plot illustrates (dis)similarity between each pair of overlapping domains. The power of this effect is indicated with colorscale. Bottom x-axis with triangles represent first partition (with triangles being TAD boundaries), top x-axis is second partition. Vertical lines help locate segments boundaries. Additionally vertical lines are coloured depending on whether they indicate overlap or not. Bottom x-axis contain ticks which indicate chromosome bining (with some Hi-C resolution) and top x-axis ticks represent segment index. Grey colored segments are discarded regions and therefore they are not indexed. Each page contain separate chromosome and at the bottom of each page a colorscale may be found.
