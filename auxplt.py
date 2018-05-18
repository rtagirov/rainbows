import matplotlib.pyplot as plt

from pylab                           import rcParams
from matplotlib                      import rc
from matplotlib.backends.backend_pdf import PdfPages

#import seaborn as sns
import numpy as np

def figpar(xtick_maj_pad = 15, ytick_maj_pad = 10, fontsize = 10):

#    sns.set(rc={"figure.figsize": figsize})
#    np.random.seed(sum(map(ord, "palettes")))

#    seaborn.palplot(seaborn.color_palette('cubehelix', 11))

#    seaborn.palplot(seaborn.diverging_palette(255, 133, l=60, n=11, center='dark'))

#    sns.palplot(sns.color_palette('dark'))

#    sns.set(font_scale = 1.5)

#    sns.set_style(style = 'white')

    rcParams['xtick.major.pad'] = xtick_maj_pad
    rcParams['ytick.major.pad'] = ytick_maj_pad

    rc('font', size = fontsize)
#    rc('font', family = 'serif')
#    rc('font', serif = 'Times')
#    rc('text', usetex = True)

# initiate figure
def initfig(nr, nc, size):

#   erase previosuly invoked figures from the memory
    plt.close('all')

#   invoke a figure with nrows x ncols subplots grid and given overall figure size
#   ax (for axes) is an object to which various methods can be applied (e.g. set_xlim, set_ylim, etc.)
    fig, ax = plt.subplots(nrows = nr, ncols = nc, figsize = size)

#   minimize unused areas in the figure
    fig.tight_layout()

    return fig, ax

# save figure as pdf in 'figdir' directory
# with 'filename' name
def savepdf(filename, figdir = './fig/'):

    pp = PdfPages(figdir + filename + '.pdf')

    pp.savefig(bbox_inches = 'tight')
    pp.close()
