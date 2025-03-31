import os
import numpy as np
import nibabel as nib
import matplotlib as mpl
import matplotlib.pyplot as plt
from figure_funcs import ff
from matplotlib import patches
import matplotlib
#matplotlib.use('Agg')
matplotlib.use("Qt5Cairo")


def fig_prf_contrast_anatomy():
    # Import and calculate mean images
    T1 = ff.load_volume(volume='T1')
    T1 = np.mean(np.array(T1), axis=0)
    ROI = ff.load_volume(volume='thalamus')
    ROI = np.array(ROI)

    # Construct figure anatomy
    slices = [100, 97, 95]
    n_slice = len(slices)
    cutoff = 0.2
    plt.figure(figsize=(4, 10.5))
    for idx_slice, slice in enumerate(slices):
        plt.subplot(n_slice, 1, idx_slice+1)
        ff.plot_slice(T1, slice, 'coronal', (0 + 30, 182 - 30), (0 + 30, 182 - 40))
        ff.plot_roi_overlay(ROI, 'coronal', slice)
        rect = patches.Rectangle(xy=(50, 55), width=((182-50)-50), height=((182-80)-55), edgecolor='black', facecolor='none')
        plt.gca().add_patch(rect)
        plt.gca().set_aspect('equal')
    # Save to disk
    plt.tight_layout()
    plt.savefig(os.path.join('../figures', 'fig_contrast_anatomy_coronal.png'))

    # Construct figure anatomy
    slices = [100, 97, 95]
    n_slice = len(slices)
    cutoff = 0.2
    plt.figure(figsize=(4, 4))
    # R2
    ff.plot_slice(T1, 85, 'sagittal', (218 - 70, 0 + 70), (0 + 35, 182 - 60))
    for slice in slices:
        plt.plot([slice, slice], [55, 55+47], color='black')
    # Save to disk
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join('../figures', 'fig_contrast_anatomy_sagittal.png'))
    plt.close('all')
    

def fig_prf_contrast_maps():
    # Import data volumes
    T1 = ff.load_volume(volume='T1')
    R2 = ff.load_volume(volume='contrastNEW_R2')
    AN = ff.load_volume(volume='contrastNEW_angle')
    EC = ff.load_volume(volume='contrastNEW_eccentricity')
    SZ = ff.load_volume(volume='contrastNEW_size')
    ROI = ff.load_volume(volume='thalamus')

    # Transform data into np.arrays and transform as necessary
    ROI = np.array(ROI)
    T1 = np.mean(np.array(T1), axis=0)
    R2 = np.median(np.array(R2), axis=0)
    AN_orig = np.real(np.angle(np.nanmedian(np.exp(1j * np.copy(AN)*np.pi/180), axis=0))*180/np.pi)
    AN = np.abs(np.angle(np.nanmedian(np.exp(1j * (np.array(AN)*np.pi/180 + np.pi/2)), axis=0)))*180/np.pi
    EC = np.median(np.array(EC), axis=0)
    SZ = np.median(np.array(SZ), axis=0)
    THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()

    # Construct figure
    ff.plot_prf_parameter_sequence(T1, R2, AN, AN_orig, EC, SZ, THA, ROI)
    plt.savefig(os.path.join('../figures', 'fig_contrast_mainmaps_updated.png'))


def fig_prf_contrast_rf_coverage():
    # Import and calculate mean images
    T1 = ff.load_volume(volume='T1')
    R2 = ff.load_volume(volume='contrastNEW_R2')
    AN = ff.load_volume(volume='contrastNEW_angle')
    EC = ff.load_volume(volume='contrastNEW_eccentricity')
    SZ = ff.load_volume(volume='contrastNEW_size')
    NSD = {'T1': T1, 'R2': R2, 'AN': AN, 'EC': EC, 'SZ': SZ}
    THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus' + '.nii.gz')).get_fdata()

    # Construct means
    NSD_mean = dict()
    NSD_mean['R2'] = np.median(np.array(NSD['R2']), axis=0)
    NSD_mean['AN'] = np.angle(np.nanmean(np.exp(1j * (np.array(NSD['AN']) * np.pi / 180)), axis=0)) * 180 / np.pi
    NSD_mean['SZ'] = np.median(np.array(NSD['SZ']), axis=0)
    NSD_mean['EC'] = np.median(np.array(NSD['EC']), axis=0)

    # Append NSD_mean to beginning of NSD
    for data_type in ['R2', 'AN', 'SZ', 'EC']:
        NSD[data_type] = [NSD_mean[data_type]] + NSD[data_type]

    # Plot + save to disk
    ff.plot_visual_field_coverage_horizontal(NSD['R2'], NSD['AN'], NSD['EC'], NSD['SZ'], include=THA, max_vox=100)
    plt.gcf().set_size_inches((12 * 1.1, 3.5 * 1.1))
    plt.savefig(os.path.join('../figures', 'fig_contrast_rf_coverage.png'))
    plt.close('all')

    # Add some extra info print out
    ff.summarize_centrality_rf_coverage(NSD['R2'], NSD['EC'], include=THA, cutoff=0.1)


if __name__ == '__main__':
    # test1.py executed as script
    # do something
    fig_prf_contrast_anatomy()
    fig_prf_contrast_maps()
    fig_prf_contrast_rf_coverage()
