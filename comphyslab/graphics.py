import os, sys
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# update fonts
FONTSIZE = 12
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : FONTSIZE}
mp.rc('font', **font)

# use latex if available on system, otherwise set usetex=False
mp.rc('text', usetex=True)

# use JavaScript for rendering animations
mp.rc('animation', html='jshtml')
# ---------------------------------------------------------------------
def plot_central_axes(ax, 
                      xmin, xmax, nxticks, xlabel,
                      ymin, ymax, nyticks, ylabel, 
                      ftsize=16):
    
    # location of xtick marks 
    xticks = np.linspace(xmin, xmax, nxticks)
    
    # get rid of the central tick mark count is odd
    if nxticks % 2 == 1:
        xticks = np.delete(xticks, nxticks//2)
        
    # location of ytick marks
    yticks = np.linspace(ymin, ymax, nyticks)
    if nyticks % 2 == 1:
        yticks = np.delete(yticks, nyticks//2)
        
    # define graph domain, tick marks, and labels
    ax.set_xlim(xmin, xmax)
    ax.set_xticks(ticks=xticks)
    ax.set_xlabel(xlabel, loc='right', fontsize=ftsize)

    # define graph range, tick marks, and labels  
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(ticks=yticks)
    ax.set_ylabel(ylabel, loc='top', fontsize=ftsize)

    # move left y-axis and bottom x-axis to center of plot
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    # eliminate upper and right axes
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
