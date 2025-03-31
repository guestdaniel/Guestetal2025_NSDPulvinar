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


def fig_prf_maps(map_name='contrastNEW', cutoff=0.1):
    # Import data volumes
    T1 = ff.load_volume(volume='T1')
    R2 = ff.load_volume(volume=map_name + '_R2')
    ROI = ff.load_volume(volume='thalamus')

    # Transform data into np.arrays and transform as necessary
    ROI = np.array(ROI)
#    T1 = np.mean(np.array(T1), axis=0)
#    R2 = np.median(np.array(R2), axis=0)
    THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()
    cutoff = 0.1

    slices = [100, 97, 95]
    n_slice = len(slices)
    plt.figure(figsize=(12, 4))
    for idx_slice, slice in enumerate(slices):
        for subj in [1, 2, 3, 4, 5, 6, 7, 8]:
            # R2
            plt.subplot(n_slice, 8, (idx_slice * 8) + subj)
            ff.plot_slice_with_overlay(T1[subj-1], R2[subj-1], slice, 'coronal', np.logical_or(R2[subj-1] < cutoff, THA != 1), (0 + 50, 182 - 50),
                                    (0 + 55, 182 - 80), cmap2=matplotlib.cm.get_cmap('hot'), clim2=(0, 1))
            ff.plot_roi_overlay(ROI, 'coronal', slice, 0.5)
    plt.tight_layout(h_pad=0.1, w_pad=0.1)
    plt.savefig('/home/daniel/Guestetal2021/../figures/subject_consistency_arcaro_' + map_name + '.png', dpi=600)

if __name__ == '__main__':
    # test1.py executed as script
    # do something
#    fig_prf_anatomy()
    map_names = ['contrastNEW', 'bodyauto', 'faceauto', 'backgroundauto', 'foregroundauto', 'salience', 'wordauto']
    for map_name in map_names:
        fig_prf_maps(map_name)
