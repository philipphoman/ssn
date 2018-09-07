#!/usr/bin/env python
import os
import sys
import mayavi
import numpy as np
import nibabel as nib
from surfer import Brain
#from matplotlib import rcParams
#import matplotlib.font_manager

# parse atlas input argument
import argparse
parser = argparse.ArgumentParser()
#parser.add_argument("a")
parser.add_argument("-a", nargs='?', default="dk")
args = parser.parse_args()
atlas = args.a

#from sys import exit
#exit(0)

  
import matplotlib.pylab as plt
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = ['Arial']
plt.rc('font', family='Arial')
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.cm as cm


# Create file stamp
from subprocess import call
call(["touch", "../output/python/fe_brains.pyout"])

def make_pysurf_fig(hemi, data, subid='fsaverage', surf='pial',
                    bg='white', cmap="Blues", cbarview=None,
                    lower=0, upper=1, savefile=True, th=0,
                    emptyval=-99, atlas='dk'):

    if atlas == 'dk':
        annot = '.aparc.annot'
    else:
        annot = '.aparc.a2009s.annot'
     

    aparc_file = os.path.join(os.environ["SUBJECTS_DIR"],
                              subid, "label",
                              hemi + annot)
    labels, ctab, names = nib.freesurfer.read_annot(aparc_file)


    vtx_data = data[labels]
    vtx_data[labels==-1] = emptyval 
    brain = Brain(subid, hemi, surf, background=bg)
    #brain = Brain(subid, 'split', surf, views=['lat', 'med'])
    brain.add_data(vtx_data, lower, upper,
                   thresh=th, colormap=cmap, alpha=.8)
    #brain.save_image("%s_lat.png" % subid)
    if (savefile==True):
        prefix = subid + '_' + hemi
        print prefix
        brain.save_imageset(prefix=prefix,
                        colorbar=cbarview,
                        views=['med', 'lat'], filetype='png')
    return;

def add_colorbar(grid, fig, cmap_name,
                 y_min=0, y_max=1, cbar_min=0,
                 cbarlabelsize=18,
                 cbarminlabel=0, cbarmaxlabel=1,
                 cbar_max=1, vert=False, label='Percentile',
                 show_ticks=False, pad=-18.5):
    '''
    Add a colorbar to the big_fig in the location defined by grid 
    
    grid       :  grid spec location to add colormap
    fig        :  figure to which colorbar will be added
    cmap_name  :  name of the colormap
    x_min      :  the minimum value to plot this colorbar between
    x_max      :  the maximum value to plot this colorbar between
    cbar_min   :  minimum value for the colormap (default 0)
    cbar_max   :  maximum value for the colormap (default 1)
    vert       :  whether the colorbar should be vertical 
    label      :  the label for the colorbar (default: None)
    ticks      :  whether to put the tick values on the colorbar 
    pad        :  how much to shift the colorbar label by (default: 0)
    '''
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap

    # Add an axis to the big_fig
    ax_cbar = plt.Subplot(fig, grid)
    fig.add_subplot(ax_cbar)
    
    # Normalise the colorbar so you have the correct upper and
    # lower limits and define the three ticks you want to show
    norm = mpl.colors.Normalize(vmin=cbar_min, vmax=cbar_max)

    if show_ticks:
        #ticks = [y_min, np.average([y_min, y_max]), y_max]
        ticks = [y_min, y_max]
    else:
        ticks=[]
            
    # Figure out the orientation
    if vert:
        orientation='vertical'
        rotation=270
    else:
        orientation='horizontal'
        rotation=0
        
    # Add in your colorbar:
    cb = mpl.colorbar.ColorbarBase(ax_cbar, cmap=cmap_name, norm=norm,
                                   orientation=orientation,
                                   ticks=ticks,
                                   boundaries=np.linspace(y_min,
                                                          y_max, 100))
         
    cb.outline.set_visible(False)                              
    cb.ax.set_xticklabels([cbarminlabel, cbarmaxlabel])

    cb.ax.tick_params(labelsize=cbarlabelsize, length=0) 
    #cb.set_family("Arial")
    #ax = cb.ax
    #text = ax.xaxis.label
    #font = plt.font_manager.FontProperties(family='Arial',
    #text.set_font_properties(font)
    #show()

    if label:
        cb.set_label(label, rotation=rotation, labelpad=pad,
                     size=cbarlabelsize, weight='bold', family='Arial')
    return fig

#-----------------------------------------------------------------------
# From Kristie Withaker
#-----------------------------------------------------------------------
def combine_pngs(fname, subid='fsaverage', cmap="Blues",
                 cbarlabel="PLS score percentile", y_min=0, y_max=1,
                 cbarminlabel=0, cbarmaxlabel=1,
                 cbarlabelsize=18):
    '''
    Find four images and combine them into one nice picture
    '''
    figsize = (5, 5)
    fig = plt.figure(figsize=figsize, facecolor='white')

    grid = gridspec.GridSpec(2, 2)
    grid.update(left=0.01, right=0.99, top=0.99, bottom=0.2, wspace=0,
                hspace=0)

    f_list = [ '_'.join([subid, 'lh', 'lat.png']),
               '_'.join([subid, 'rh', 'lat.png']),
               '_'.join([subid, 'lh', 'med.png']),
               '_'.join([subid, 'rh', 'med.png']) ]

    # Plot each figure in turn
    for g_loc, f in zip(grid, f_list):
        ax = plt.Subplot(fig, g_loc)
        fig.add_subplot(ax)
        img = mpimg.imread(f)

        # Crop the figures appropriately NOTE: this can change depending
        # on which system you've made the images on originally - it's a
        # bug that needs to be sorted out!
        if 'lat' in f:
            img_cropped = img[130:(-130), 5:(-40), :]
        else:
            img_cropped = img[120:(-100), :, :]
        ax.imshow(img_cropped, interpolation='none')
        #ax.imshow(img, interpolation='none')
        ax.set_axis_off()

    # Add the bottom of one of the images as the color bar
    # at the bottom of the combo figure
    grid_cbar = gridspec.GridSpec(1, 1)
    grid_cbar.update(left=0.3, right=0.7, top=0.18, bottom=0.15,
                     wspace=0, hspace=0)

    # Save the figure
    #filename = subid + '_' + 'combined.pdf'
        
    fig = add_colorbar(grid_cbar[0], fig, cmap, show_ticks=True,
                       cbarminlabel=cbarminlabel,
                       cbarmaxlabel=cbarmaxlabel,
                       cbarlabelsize=cbarlabelsize,
                       pad=3, label=cbarlabel, y_max=y_max, y_min=y_min)
    print fig
    fig.savefig(fname, bbox_inches=0, dpi=300)


def plot_dk_parc(subid='fsaverage', atlas='dk'):
    import os
    from os.path import join as pjoin
    from surfer import Brain

    subid = 'fsaverage'
    hemi = 'lh'
    surf = 'pial'
    view = 'lateral'

    if atlas == 'dk':
        annot = 'aparc'
    else:
        annot = 'aparc.a2009s'
        

    """
    Bring up the visualization
    """
    brain = Brain(subid, hemi, surf, views=view,
                  cortex="bone", background="white")
    brain.add_annotation(annot, borders=False)
    prefix = 'fe_freesurfer_dkparc_lh'
    brain.save_imageset(prefix=prefix,
                        views=['med', 'lat'], filetype='png')
    return;


# Load percentile scores for each region
plsfiles = ['../output/tables/fe_freesurfer_pyinput_plsweights1_'
            + atlas + '.csv',
            '../output/tables/fe_freesurfer_pyinput_plsweights2_'
            + atlas + '.csv']
plspdfs = ['../output/figures/fe_freesurfer_pls1brain_' 
           + atlas + '.png',
           '../output/figures/fe_freesurfer_pls2brain_'
           + atlas + '.png']
cmap = "Blues"
for f, p  in zip(plsfiles, plspdfs):
    from numpy import genfromtxt
    mydata = genfromtxt(f, delimiter=',')
    x = np.array([0])
    ml = len(mydata)
    rhdata = np.insert(mydata[ml/2:ml], 3, x)
    rhdata = np.insert(rhdata, 0, x)
    lhdata = np.insert(mydata[0:ml/2], 3, x)
    lhdata = np.insert(lhdata, 0, x)
    make_pysurf_fig(hemi='lh', surf='pial', data=lhdata, cmap=cmap,
                    atlas=atlas)
    make_pysurf_fig(hemi='rh', surf='pial', data=rhdata, cmap=cmap,
                    atlas=atlas)
    combine_pngs(fname=p, cmap=cmap, cbarlabel="PLS score percentile")
    plot_dk_parc(atlas='destrieux')

    
# Create t map off group differences in cortical thickness
from numpy import genfromtxt
mydata = genfromtxt("../output/tables/fe_freesurfer_pyinput_roitvals_"
                    + atlas + ".csv", delimiter=',')
x = np.array([0])
ml = len(mydata)
rhdata = np.insert(mydata[ml/2:ml], 3, x)
rhdata = np.insert(rhdata, 0, x)
lhdata = np.insert(mydata[0:ml/2], 3, x)
lhdata = np.insert(lhdata, 0, x)

# This will always use DK atlas for now
make_pysurf_fig(hemi='lh', surf='pial', data=lhdata, cmap="RdBu_r",
                th=min(mydata), atlas='dk',
                lower=min(mydata), upper=abs(min(mydata)), savefile=True)

make_pysurf_fig(hemi='rh', surf='pial', data=rhdata, cmap="RdBu_r",
                th=min(mydata), atlas='dk',
                lower=min(mydata), upper=abs(min(mydata)), savefile=True)


combine_pngs(fname="../output/figures/fe_freesurfer_tvals_brain_" +
             atlas + ".pdf",
             cmap="RdBu_r",
             y_min=0, y_max=1, cbarminlabel=round(min(mydata), 1),
             cbarmaxlabel=round(abs(min(mydata)), 1), cbarlabel="t",
             cbarlabelsize=12)


# Create  map of group differences in nodal degree
from numpy import genfromtxt
mydata = genfromtxt("../output/tables/fe_freesurfer_pyinput_deltadegrees_sig_" + atlas + ".csv",
                    delimiter=',')
x = np.array([0])
ml = len(mydata) 
rhdata = np.insert(mydata[ml/2:ml], 3, x)
rhdata = np.insert(rhdata, 0, x)
lhdata = np.insert(mydata[0:ml/2], 3, x)
lhdata = np.insert(lhdata, 0, x)

make_pysurf_fig(hemi='lh', surf='pial', data=lhdata, cmap="RdBu_r",
                th=-max(lhdata), atlas=atlas,
                lower=-max(mydata), upper=max(mydata), savefile=True)

make_pysurf_fig(hemi='rh', surf='pial', data=rhdata, cmap="RdBu_r",
                th=-max(rhdata), atlas=atlas,
                lower=-max(mydata), upper=max(mydata), savefile=True)

combine_pngs(fname="../output/figures/fe_freesurfer_deltadegrees_brain_" +
             atlas + ".pdf",
             cmap="RdBu_r",
             y_min=0, y_max=1, cbarminlabel=round(-max(mydata), 1),
             cbarmaxlabel=round(max(mydata), 1), cbarlabel="Degree",
             cbarlabelsize=12)
