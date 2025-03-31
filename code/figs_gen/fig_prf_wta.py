import os
import numpy as np
import nibabel as nib
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from figure_funcs import ff


def fig_prf_wta_mip_maps():
    # Load in all data and store it in a list in order
    data = list()
    maps = ['contrastNEW', 'salience', 'foregroundauto', 'faceauto', 'bodyauto']
    for map in maps:
        data.append(np.mean(np.array(ff.load_volume(volume=map + '_R2')), axis=0))
    T1 = np.mean(np.array(ff.load_volume(volume='T1')), axis=0)
    roi = ff.load_volume(volume='thalamus')
    postthal = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus' + '.nii.gz')).get_fdata()

    plt.figure(figsize=(8, 2))
    for idx_map, map in enumerate(data):
        plt.subplot(1, len(maps), idx_map+1)
        ff.plot_mip_roi_overlay(plt.gcf().gca(), roi, 'coronal', linewidth=0.7, linestyle=(0, (3, 1)))
        map[postthal != 1] = 0
        ff.plot_mip(map, 'coronal', (0 + 50, 182 - 50), (0 + 55, 182 - 80), cmap1='hot', clim1=(0, 1))
        #plt.tight_layout()
    plt.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.1)
    plt.savefig(os.path.join('../figures', 'fig_prf_wta_mip.png'), dpi=300)


def fig_prf_wta_contextual_anatomy():
    # Load anatomical data
    T1 = np.mean(np.array(ff.load_volume(volume='T1')), axis=0)
    roi = ff.load_volume(volume='thalamus')

    # Set slices
    slices = [100, 97, 95, 93]

    # Plot sagittal T1 with markers at slice positions
    plt.figure(figsize=(4, 4))
    ff.plot_slice(T1, 105, 'sagittal', lims_x=(145, 65), lims_y=(35, 120))
    for slice in slices:
        plt.plot([slice, slice], [0, 182], color='red')
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join('../figures', 'fig_prf_wta_anatomy.png'), dpi=300)

    # Plot coronal T1 with ROI labels
    # Construct figure 1
    slices = [100, 97, 95]
    plt.figure()
    for idx_slice, slice in enumerate(slices):
        plt.subplot(len(slices), 1, idx_slice+1)
        ff.plot_slice(T1, slice, 'coronal', (0 + 50, 182 - 50), (0 + 55, 182 - 80))
        ff.plot_roi_overlay(roi, 'coronal', slice, outline_width=2.0)
    plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.25)
    plt.savefig(os.path.join('../figures', 'fig_prf_wta_anatomy_coronal.png'), dpi=300)


def fig_prf_wta_maps():
    # Load in all data and store it in a list in order
    data = list()
    maps = ['contrastNEW', 'salience', 'foregroundauto', 'bodyauto', 'faceauto']
    for map in maps:
        data.append(np.mean(np.array(ff.load_volume(volume=map + '_R2')), axis=0))
    T1 = np.mean(np.array(ff.load_volume(volume='T1')), axis=0)
    THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()

    # Compute WTA map
    WTA = np.argmax(np.array(data), axis=0)
    R2_max = np.max(np.array(data), axis=0)

    # Construct figure 1
    slices = [100, 97, 95, 93]
    n_slice = len(slices)
    cutoff = 0.2
    plt.figure(figsize=(12, 2.5))
    for idx_slice, slice in enumerate(slices):
        # R2
        plt.subplot(1, n_slice, idx_slice+1)
        ff.plot_slice_with_overlay(T1, WTA, slice, 'coronal', np.logical_or(R2_max < cutoff, THA != 1),
                                   (0 + 50, 182 - 50), (0 + 55, 182 - 80), cmap2=ff.cmap_feat2, clim2=(0, 4))
    plt.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.1)
    plt.savefig(os.path.join('../figures', 'fig_prf_wta_maps.png'), dpi=300)