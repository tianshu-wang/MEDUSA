import numpy as np
import pylab
from cycler import cycler

FIG_WIDTH_PT   = 720.0 #360.0
FONT_SIZE      = 14

LINEWIDTH      = cycler('lw',range(2,4))
LINESTYLES     = cycler('ls',['-','--','-.',':'])
LINECOLORS     = cycler('color', ['Blue','Red','Green','Magenta','Goldenrod','DeepSkyBlue','DarkViolet',
                  'OrangeRed','BlueViolet','SaddleBrown'])
STYLES         = LINECOLORS*LINEWIDTH

def set_style():
    inches_per_pt = 1.0/72.27                   # Convert pt to inches
    golden_mean   = (np.sqrt(5)-1.0)/2.0        # Aesthetic ratio
    fig_width     = FIG_WIDTH_PT*inches_per_pt  # width in inches
    fig_height    = fig_width*golden_mean       # height in inches
    fig_size      = [fig_width,fig_height]

    params = {'backend'              : 'ps',
              'font.family'          : 'serif',
              'axes.prop_cycle'      : STYLES,
              'axes.labelsize'       : FONT_SIZE,
              'font.size'            : FONT_SIZE,
              'legend.fontsize'      : FONT_SIZE - 2,
              'legend.borderpad'     : 0.2,
              'legend.handletextpad' : 0.5,
              'legend.labelspacing'  : 0.1,
              'lines.linewidth'      : 2,
              'xtick.labelsize'      : FONT_SIZE - 2,
              'ytick.labelsize'      : FONT_SIZE - 2,
              'savefig.bbox'         : 'tight',
              'figure.figsize'       : fig_size}

    pylab.rcParams.update(params)
