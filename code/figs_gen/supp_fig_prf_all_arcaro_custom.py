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


def fig_prf_maps(map_name='contrastNEW', cutoff=0.1, flag=''):
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
    ff.plotsave_prf_parameter_sequence_reduced_arcaro_custom(T1, R2, AN, AN_orig, EC, SZ, THA, ROI, os.path.join('../figures', 'fig_prfall_arcaro_custom_' + flag + '_' + map_name + '_'), outline_width=5.0, cutoff=cutoff)


if __name__ == '__main__':
    # test1.py executed as script
    # do something
#    fig_prf_anatomy()
    map_names = ['contrastNEW', 'bodyauto', 'faceauto', 'backgroundauto', 'foregroundauto', 'salience', 'wordauto']
    for map_name in map_names:
        fig_prf_maps(map_name)
    fig_prf_maps(cutoff=5.0, flag='nodata')
