import os
import numpy as np
import matplotlib.pyplot as plt
from figure_funcs import ff
from matplotlib import patches
import matplotlib
matplotlib.use('Agg')
import nibabel as nib


def fig_intro_anatomy():
    # Import and calculate mean images
    T1 = ff.load_volume(volume='T1', subjs=[1])[0]
    #T1 = np.mean(np.array(T1), axis=0)

    # Create axes and subplots
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(7, 2))

    # Define which slices we're going to plot in each view (tuple of slices, corresponding to tuple of views)
    slices = [102, 95, 70]
    views = ['coronal', 'coronal', 'sagittal']

    # To compute xlims and ylims, we define a central and viewbox around that centroid
    centroid = [92, 110, 80]
    viewbox = 80  # voxels
    viewbox_narrow = 35  # voxels

    def get_xlim_ylim(centroid, viewbox):
        xlims = [(centroid[0]-viewbox, centroid[0]+viewbox),  # x-axis
                 (centroid[0]-viewbox, centroid[0]+viewbox),  # x-axis
                 (centroid[0]+viewbox, centroid[0]-viewbox)]  # y-axis, reversed
        ylims = [(centroid[2]-viewbox, centroid[2]+viewbox),  # z-axis
                 (centroid[2]-viewbox, centroid[2]+viewbox),  # z-axis
                 (centroid[2]-viewbox, centroid[2]+viewbox)]  # z-axis
        return xlims, ylims

    xlims, ylims = get_xlim_ylim(centroid, viewbox)
    xlims_narrow, ylims_narrow = get_xlim_ylim(centroid, viewbox_narrow)

    # Loop through and plot big viewbox
    for slice, view, ax, xlim, ylim, xlim_narrow, ylim_narrow, viewbox_color \
        in zip(slices, views, axs, xlims, ylims, xlims_narrow, ylims_narrow, [np.array([199, 0, 255])/255,
                                                                              np.array([255, 185, 0])/255,
                                                                              np.array([32, 183, 0])/255]):
        # Plot T1 on bottom row
        ff.plot_slice_oo(ax, T1, slice, view, lims_x=xlim, lims_y=ylim)
        print(view)
        # Plot red dashed line to indicate other slice positions
        if view == 'coronal':
            ax.plot([slices[2], slices[2]], [0, 182], linestyle='dashed', linewidth=1.0, color=np.array([32, 183, 0])/255)
        elif view == 'sagittal':
            ax.plot([slices[0], slices[0]], [0, 218], linestyle='dashed', linewidth=1.0, color=np.array([199, 0, 255])/255)
            ax.plot([slices[1], slices[1]], [0, 218], linestyle='dashed', linewidth=1.0, color=np.array([255, 185, 0])/255)
        # Plot red box around narrow viewbox
        rect = patches.Rectangle((xlim_narrow[0], ylim_narrow[0]), xlim_narrow[1]-xlim_narrow[0],
                                 ylim_narrow[1]-ylim_narrow[0], linewidth=1, edgecolor=viewbox_color, facecolor='none')
        ax.add_patch(rect)

    # Save this figure out to disk
    plt.savefig(os.path.join(ff.dir_plots, 'fig_intro_anatomy_1.png'), dpi=300)

    # Loop through and plot small viewbox
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(7, 2))
    for slice, view, ax, xlim_narrow, ylim_narrow in zip(slices, views, axs, xlims_narrow, ylims_narrow):
        # Plot T1 on bottom row
        ff.plot_slice_oo(ax, T1, slice, view, lims_x=xlim_narrow, lims_y=ylim_narrow)

    # Save this figure out to disk
    plt.savefig(os.path.join(ff.dir_plots, 'fig_intro_anatomy_2.png'), dpi=300)

if __name__ == '__main__':
    fig_intro_anatomy()


