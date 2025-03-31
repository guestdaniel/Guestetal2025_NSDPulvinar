import os
import numpy as np
import nibabel as nib
import matplotlib as mpl
import matplotlib.pyplot as plt
from figure_funcs import ff
from matplotlib import patches
import matplotlib
matplotlib.use('Agg')
plt.ion()
from figs_gen.fig_corr_cor_to_sub.fig_corr_cor_to_sub import clean_data


def fig_corr_methods_make_fake_data(n_samp=10, N=30):
    colors = ['red', 'green', 'blue', 'gray']
    for color in colors:
        for sample in range(n_samp):
            plt.figure(figsize=(2, 0.5))
            x = np.random.normal(0, 1, N)
            plt.plot(x, color=color)
            plt.axis('off')
            plt.tight_layout(pad=0.0)
            plt.savefig(os.path.join(ff.dir_plots, 'fig_corr_methods_fake_' + color + '_' + str(sample) + '.png'), dpi=300)
            plt.close()


def fig_corr_methods_plot_anatomy():
    # Set parameters
    shave_x = 40
    shave_y = 40
    # Plot example coronal T1 slice
    T1 = ff.load_volume(subjs=[1], space='mni', volume='T1')
    plt.figure(figsize=(4, 4))
    ff.plot_slice(T1[0], 95, 'coronal', lims_x=(0+shave_x, 182-shave_x), lims_y=(0+shave_y-10, 182-shave_y-10))
    plt.gca().set_aspect('equal')
    plt.tight_layout(0.0)
    plt.savefig(os.path.join(ff.dir_plots, 'fig_corr_methods_anatomy_subcortical.png'), dpi=300, inetrpolation=None)


def fig_corr_methods_plot_example_outputs():
    # Set parameters
    shave_x = 40
    shave_y = 40
    # Plot example coronal T1 slice
    T1 = ff.load_volume(volume='T1', subjs=[1])
    THA = ff.load_volume(volume='postthalamus')
    T1 = T1[0]

    # Indicate which maps you want to load
    maps_to_load = [0, 5]

    # Load in and store average data
    maps = []
    for label in maps_to_load:
        temp = []
        for subj in range(8):
            # Load data
            temp.append(clean_data(nib.load(os.path.join(ff.dir_data, 'subj0' + str(subj + 1), 'mni',
                                                         'corr_cor_to_sub' '_hemi_' + str(0 + 1) + '_label_' + str(
                                                             label + 1) + '_method_2.nii.gz')).get_fdata(),
                                   THA[subj]))
        maps.append(temp)

    # Set up colors
    label_colors = list(matplotlib.cm.plasma(np.linspace(0, 1, 14)))
    cortical_roi_labels = ['V1', 'V2', 'V3', 'hV4', 'OFA', 'FFA', 'aTL-f', 'EBA', 'FBA', 'PPA', 'VWFA', 'MT', 'IPS', 'FEF']

    # Create plot
    for idx_map, map in enumerate(maps):
        # Create figure
        plt.figure(figsize=(4, 4))
        # Average across subjects
        prob_map = np.sum(np.array(map) > 0.01, axis=0)
        # Plot slice
        ff.plot_slice_with_overlay(T1, prob_map, 95, 'coronal', mask=(prob_map <= 0), lims_x=(0+shave_x, 182-shave_x), lims_y=(0+shave_y-10, 182-shave_y-10), cmap2='plasma', clim2=(0, 8))
        # Handle aspect and layout
        plt.gca().set_aspect('equal')
        plt.tight_layout(0.0)
        # Save figure to disk
        plt.savefig(os.path.join('../figures', 'fig_corr_methods_example_outputs_' + str(idx_map) + '.png'), dpi=300,
                    interpolation='none')
        # Close figure
        plt.close()