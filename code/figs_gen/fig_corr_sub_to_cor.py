import numpy as np
import os
import matplotlib.pyplot as plt
from figure_funcs import ff
import matplotlib
matplotlib.use('Agg')


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def plot_view_sequence(data, point, coords_start, clim=None, cmap='gray', roi=None, marker_color='black'):
    # Set up figure
    plt.figure(figsize=(8, 3))

    # Define everything we're going to loop over --- which slices, marker angles and positions, limits, views
    slices = [int(coords_start[1] + point[1] - 2),
              int(coords_start[0] + point[0] - 2),
              int(coords_start[2] + point[2] - 2)]
    angles = [5*np.pi/4, 1*np.pi/4, np.pi/4]
    markers = [[(coords_start[0] + point[0] - 2), (coords_start[2] + point[2] - 2)],
               [(coords_start[1] + point[1] - 2), (coords_start[2] + point[2] - 2)],
               [(coords_start[0] + point[0] - 2), (coords_start[1] + point[1] - 2)]]
    xlims = [(30, 80), (110, 60), (30, 80)]
    ylims = [(50, 100), (50, 100), (50, 100)]
    views = ['coronal', 'sagittal', 'axial']

    # Loop through everything and visualize
    for ax, slice, angle, marker, xlim, ylim, view in zip(range(3), slices, angles, markers, xlims, ylims, views):
        plt.subplot(1, 3, ax+1)
        ff.plot_slice(data, slice, view, lims_x=xlim, lims_y=ylim, cmap1=cmap, clim1=clim, marker=marker,
                      marker_angle=angle, marker_color=marker_color, marker_length=10)
        if roi is not None:
            ff.plot_roi_overlay(roi, slice, view)

    # Tighten up the layout
    plt.tight_layout()


def plot_fig_corr_sub_to_cor_r2_and_anatomy_maps():
    # Load T1, PRF R2, and ROI labels
    T1 = ff.load_volume(subjs=[1], space='func1mm', volume='T1')[0]
    contrast_temp = ff.load_volume(subjs=[1], space='func1mm', volume='contrastNEW_R2')[0]
    body_temp = ff.load_volume(subjs=[1], space='func1mm', volume='bodyauto_R2')[0]
    roi = ff.load_volume(subjs=[1], space='func1mm', volume='thalamus')[0]

    # Manually define a few things
    coords_start = [45, 66, 58]  # manually define the LPI corner of the extent of the PRF data in the func1mm volume
    coords_end = [102, 96, 89]  # manually define the RAS corner of the extent of the PRF data in the func1mm volume
    contrast_peak = [11, 15, 12]  # manually define location of contrast peak in LH Subj 1 relative to LPI corner
    body_peak = [18, 14, 15]  # manually define location of body peak in LH Subj 1 relative to LPI corner

    # Embed contrast/body PRF data in func1mm-sized volumes
    contrast = np.zeros_like(T1)
    body = np.zeros_like(T1)
    contrast[(coords_start[0] - 1):coords_end[0], (coords_start[1] - 1):coords_end[1],
    (coords_start[2] - 1):coords_end[2]] = contrast_temp
    body[(coords_start[0] - 1):coords_end[0], (coords_start[1] - 1):coords_end[1],
    (coords_start[2] - 1):coords_end[2]] = body_temp

    # Plot view sequences
    plot_view_sequence(T1, contrast_peak, coords_start)
    plt.savefig(os.path.join('../figures', 'fig_corr_sub_to_cor_anat_contrast.png'), dpi=300)
    plot_view_sequence(contrast, contrast_peak, coords_start, cmap='hot', clim=(0, 2.5), marker_color='cyan')
    plt.savefig(os.path.join('../figures', 'fig_corr_sub_to_cor_r2_contrast.png'), dpi=300)

    plot_view_sequence(T1, body_peak, coords_start)
    plt.savefig(os.path.join('../figures', 'fig_corr_sub_to_cor_anat_body.png'), dpi=300)
    plot_view_sequence(body, body_peak, coords_start, cmap='hot', clim=(0, 2.5), marker_color='cyan')
    plt.savefig(os.path.join('../figures', 'fig_corr_sub_to_cor_r2_body.png'), dpi=300)

    # Save dummy colorbars
    a = np.array([[0, 2.5]])
    plt.figure(figsize=(1*0.75, 2.25*0.75))
    img = plt.imshow(a, cmap='hot')
    plt.gca().set_visible(False)
    cax = plt.axes([0.2, 0.1, 0.2, 0.8])
    plt.colorbar(cax=cax)
    plt.savefig(os.path.join('../figures', 'fig_corr_sub_to_cor_r2_colormap.png'), dpi=300)
