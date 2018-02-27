# max row size in bins
ROW_SIZE = 80

# page (row) width in inches
FIG_WIDTH = 10

# row height in inches (so page height is FIG_HEIGHT * number of rows)
FIG_HEIGHT = 2.5

# how often to put bin coordinate label, every BIN_LABEL_FREQ bins
BIN_LABEL_FREQ = 10

# max and min segmentation, if None this will be infered from data
# for BP score it is recommended to fix this as [0,1]
DATA_RANGE = [0,1] #None

# heatmap pallete, for set of available palletes see:
# https://matplotlib.org/examples/color/colormaps_reference.html
# if you want to reverse color append '_r' at the end
HEATMAP_COLOR = 'PiYG_r' #'jet'

# what color should unmappable regions be colored with
UNMAPPABLE_COLOR = 'gray'

# color for overlapping boundaries
OVERLAPPING_BOUNDARIES_COLOR = 'black'

# color for non-overlapping boundaries
NONOVERLAPPING_BOUNDARIES_COLOR = 'blue'
