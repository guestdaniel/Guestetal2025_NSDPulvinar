import os
import numpy as np
import nibabel as nib
import matplotlib as mpl
import matplotlib.pyplot as plt
from figure_funcs import ff
from skimage import measure

# Set non-interactive backend
mpl.use('Agg')

# Import data volumes
anatomical_data = dict()
for data_type in ['T1', 'T2', 'EPI', 'SWI', 'TOF']:
    anatomical_data[data_type] = ff.load_volume(volume=data_type)
functional_data = dict()
for data_type in ['contrastNEW', 'faceauto', 'bodyauto']:
    functional_data[data_type] = ff.load_volume(volume=data_type + '_' + 'R2')

# Set parameters
slices = list(range(88, 108))      # which coronal slices to view
shave_x = 53                       # how many voxels to shave in x-dir from both sides
shave_y = 48                       # how many voxels to shave in y-dir from both sides
shave_y_top = 40                   # how much additional to shave off top of y
clims = dict()                     # color limits of each anatomical data type
clims['T1'] = (0, 1700)
clims['T2'] = (0, 900)
clims['EPI'] = (0, 1700)
clims['SWI'] = (0, 900)
clims['TOF'] = (0, 1500)


def plot_slice_with_contours(ax, anatomy, functional, shave_x, shave_y, shave_y_top, clim, cutoffs, outline_color=None):
    # Plot functional data
    ff.plot_slice_oo(ax, anatomy, slice, 'coronal', lims_x=(0 + shave_x, 182 - shave_x),
                  lims_y=(0 + shave_y, 182 - shave_y - shave_y_top), clim1=clim)
    plt.axis('off')
    # Plot R2 contours
    for idx_cutoff, R2_cutoff in enumerate(cutoffs):
        R2_slice = functional > R2_cutoff
        R2_slice = np.flip(np.squeeze(R2_slice[:, slice, :]), axis=0)
        contour = measure.find_contours(R2_slice, 0.5)
        if len(contour) > 0:
            for cont in contour:
                if outline_color is not None:
                    ax.plot(cont[:, 0], cont[:, 1], linewidth=1, color=outline_color)  # fix colormap
                else:
                    ax.plot(cont[:, 0], cont[:, 1], linewidth=1,
                            color=mpl.cm.hot(np.linspace(0, 1, len(cutoffs)))[idx_cutoff])  # fix colormap


# ===== Version 1 =====
# In this version, we plot four R2 contours and every subject in a single plot
R2_cutoffs = [0.5, 1.0, 1.5, 2.0]  # which R2 cutoffs to draw contours at
# Loop through functional data types
for idx_func, dt_func in enumerate(['contrastNEW', 'faceauto', 'bodyauto']):
    # Loop through slices
    for slice in slices:
        f, axs = plt.subplots(5, 8, figsize=(12 * 1.4, 6.5 * 1.1))
        # Loop through subjects
        for subj in range(8):
            # Loop through anatomical data types
            for idx_anat, dt_anat in enumerate(['T1', 'T2', 'EPI', 'SWI', 'TOF']):
                plot_slice_with_contours(axs[idx_anat][subj], anatomical_data[dt_anat][subj],
                                         functional_data[dt_func][subj], shave_x, shave_y, shave_y_top, clims[dt_anat],
                                         R2_cutoffs)
        plt.tight_layout()
        plt.title('Slice ' + str(slice))
        plt.savefig(os.path.join('../figures_bulk', 'activation_vs_anatomy_' + 'slice_' + str(slice) + '_' + dt_func + '.png'))
        plt.close()

# ===== Version 2 =====
# In this version, we plot a single R2 contours at 1% and every subject in a single plot
R2_cutoffs = [1]  # which R2 cutoffs to draw contours at
# Loop through functional data types
for idx_func, dt_func in enumerate(['contrastNEW', 'faceauto', 'bodyauto']):
    # Loop through slices
    for slice in slices:
        f, axs = plt.subplots(5, 8, figsize=(12 * 1.4, 6.5 * 1.1))
        # Loop through subjects
        for subj in range(8):
            # Loop through anatomical data types
            for idx_anat, dt_anat in enumerate(['T1', 'T2', 'EPI', 'SWI', 'TOF']):
                plot_slice_with_contours(axs[idx_anat][subj], anatomical_data[dt_anat][subj], functional_data[dt_func][subj], shave_x,
                                         shave_y, shave_y_top, clims[dt_anat], R2_cutoffs, outline_color='red')
        plt.tight_layout()
        plt.title('Slice ' + str(slice))
        plt.savefig(os.path.join('../figures_bulk', 'activation_vs_anatomy_' + 'slice_' + str(slice) + '_' + dt_func + 'simple.png'))
        plt.close()

# ===== Version 3 =====
# In this version, we plot a single prob map contours and average across subjects
prob_map_cutoffs = [4, 5, 6, 7, 8]  # which probmap cutoffs to draw contours at
prob_map_criterion = 0.1            # R2 criterion to use for constructing prob maps
# Loop through functional data types
for idx_func, dt_func in enumerate(['contrastNEW', 'faceauto', 'bodyauto']):
    # Construct probmap
    prob_map = np.sum(np.array(functional_data[dt_func]) > prob_map_criterion, axis=0)
    # Loop through slices
    for slice in slices:
        f, axs = plt.subplots(1, 5, figsize=(17, 6))
        for idx_anat, dt_anat in enumerate(['T1', 'T2', 'EPI', 'SWI', 'TOF']):
            # Construct avg anat
            avg_anat = np.mean(np.array(anatomical_data[dt_anat]), axis=0)
            plot_slice_with_contours(axs[idx_anat], avg_anat, prob_map, shave_x,
                                     shave_y, shave_y_top, clims[dt_anat], prob_map_cutoffs)
        plt.tight_layout()
        plt.title('Slice ' + str(slice))
        plt.savefig(os.path.join('../figures_bulk', 'activation_vs_anatomy_' + 'slice_' + str(slice) + '_' + dt_func + 'avg.png'))
        plt.close()
