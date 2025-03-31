import os
from re import L
import numpy as np
import nibabel as nib
import matplotlib as mpl
import matplotlib.pyplot as plt
from figure_funcs import ff
from matplotlib import patches
import matplotlib
matplotlib.use('Agg')


def fig_prf_anatomy():
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
    plt.savefig(os.path.join('../figures', 'fig_prf_anatomy_coronal.png'))

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
    plt.savefig(os.path.join('../figures', 'fig_prf_anatomy_sagittal.png'))
    plt.close('all')

    

def fig_prf_maps(map_name='contrastNEW'):
    # Import data volumes
    T1 = ff.load_volume(volume='T1')
    R2 = ff.load_volume(volume=map_name + '_R2')
    AN = ff.load_volume(volume=map_name + '_angle')
    EC = ff.load_volume(volume=map_name + '_eccentricity')
    SZ = ff.load_volume(volume=map_name + '_size')
    ROI = ff.load_volume(volume='thalamus')

    # Transform data into np.arrays and transform as necessary
    ROI = np.array(ROI)
    T1 = np.mean(np.array(T1), axis=0)
    R2 = np.median(np.array(R2), axis=0)
    if map_name == 'salience':
        AN_orig = np.real(np.angle(np.nanmedian(np.exp(1j * np.mod(270-np.copy(AN), 360)*np.pi/180), axis=0))*180/np.pi)
        AN = np.abs(np.angle(np.nanmedian(np.exp(1j * (np.mod(270-np.array(AN), 360)*np.pi/180 + np.pi/2)), axis=0)))*180/np.pi
    else:
        AN_orig = np.real(np.angle(np.nanmedian(np.exp(1j * np.copy(AN)*np.pi/180), axis=0))*180/np.pi)
        AN = np.abs(np.angle(np.nanmedian(np.exp(1j * (np.array(AN)*np.pi/180 + np.pi/2)), axis=0)))*180/np.pi
    EC = np.median(np.array(EC), axis=0)
    SZ = np.median(np.array(SZ), axis=0)
    THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()

    # Construct figure
    ff.plotsave_prf_parameter_sequence_reduced(T1, R2, AN, AN_orig, EC, SZ, THA, ROI, os.path.join('../figures', 'fig_prfall_' + map_name + '_'))


if __name__ == '__main__':
    # test1.py executed as script
    # do something
    fig_prf_anatomy()
    map_names = ['contrastNEW', 'bodyauto', 'faceauto', 'backgroundauto', 'foregroundauto', 'salience', 'wordauto']
    for map_name in map_names:
        fig_prf_maps(map_name)
