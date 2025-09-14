import numpy as np
import os
import matplotlib.pyplot as plt
from figure_funcs import ff
import matplotlib
# matplotlib.use('Agg')


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


def plot_fig_corr_sub_to_cor_pRF_data():
    """ Plot pRF data for seed voxels in Figs 4/5

    Takes the seed voxel locations (hardcoded) used in the correlation figures 
    """
    # First, hard code thje locations of the seed voxels in the contrast and body pRF data.
    # One location for each subject. These come from `f_identify_candidate_voxels_subject_native.m`.
    # These are relative to LPI corner of the subset of the func1mm volume used in the subcortical
    # pRF fitting.
    contrast_peak_lh = [
        [11, 15, 12],
        [10, 19, 9],
        [10, 17, 11],
        [8, 21, 10],
        [7, 19, 10],
        [9, 21, 10],
        [9, 17, 13],
        [6, 18, 10]
    ]

    contrast_peak_rh = [
        [48, 18, 16],
        [49, 22, 9],
        [51, 21, 16],
        [53, 20, 12],
        [47, 20, 9],
        [52, 18, 12],
        [51, 18, 11],
        [51, 17, 11]
    ]

    body_peak_lh = [
        [18, 14, 15],
        [19, 15, 14],
        [9, 18, 10],
        [16, 16, 15],
        [15, 15, 15],
        [16, 16, 13],
        [17, 13, 17],
        [15, 14, 11],
    ]

    body_peak_rh = [
        [39, 16, 15],
        [42, 17, 14],
        [44, 17, 20],
        [45, 16, 17],
        [40, 16, 15],
        [44, 14, 15],
        [46, 14, 16],
        [44, 14, 14],
    ]

    # Create figure
    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_axes(111)

    # Extract pRF data for the contrast peak and add to plot
    AN = ff.load_volume(volume='contrastNEW_angle', space='func1mm')
    EC = ff.load_volume(volume='contrastNEW_eccentricity', space='func1mm')
    SZ = ff.load_volume(volume='contrastNEW_size', space='func1mm')
    for subj in range(8):
        contrast_peak = contrast_peak_lh[subj]
        contrast_peak = [coord - 1 for coord in contrast_peak]  # Convert to zero-based indexing
        ang = AN[subj][contrast_peak[0], contrast_peak[1], contrast_peak[2]]  
        ecc = EC[subj][contrast_peak[0], contrast_peak[1], contrast_peak[2]]
        siz = SZ[subj][contrast_peak[0], contrast_peak[1], contrast_peak[2]]
        add_pRF_outline(ax, ang, ecc, siz, color=[0.8, 0.2, 0.6, 1.0], include_decorations=(subj == 0))

    # Extract pRF data for the body peak and add to plot
    AN = ff.load_volume(volume='bodyauto_angle', space='func1mm')
    EC = ff.load_volume(volume='bodyauto_eccentricity', space='func1mm')
    SZ = ff.load_volume(volume='bodyauto_size', space='func1mm')
    for subj in range(8):
        body_peak = body_peak_lh[subj]
        body_peak = [coord - 1 for coord in body_peak]  # Convert to zero-based indexing
        ang = AN[subj][body_peak[0], body_peak[1], body_peak[2]]
        ecc = EC[subj][body_peak[0], body_peak[1], body_peak[2]]
        siz = SZ[subj][body_peak[0], body_peak[1], body_peak[2]]
        add_pRF_outline(ax, ang, ecc, siz, color=[0.9, 0.8, 0.0, 1.0], include_decorations=(subj == 0))

    fig.show()
    plt.savefig(os.path.join('../figures', 'fig_corr_sub_to_cor_seed_pRF_properties.png'), dpi=300)


def add_pRF_outline(ax, ang, ecc, siz, color=[1, 0, 0], include_decorations=True):
    """ Add pRF outline to the given axis """
    # Lay out basic markings on plot
    if include_decorations:
        ax.plot([-10, 10], [0, 0], linestyle='dashed', color='gray')
        ax.plot([0, 0], [-10, 10], linestyle='dashed', color='gray')
        rect = plt.Rectangle(xy=(-4.2, -4.2), width=8.4, height=8.4, color='gray', fill=False)
        ax.add_artist(rect)

        # Set limits
        ax.set_xlim((-8.4, 8.4))
        ax.set_ylim((-8.4, 8.4))
        ax.set_xticks([])
        ax.set_yticks([])

        # Force square axes
        ax.set(adjustable='box', aspect='equal')

    # Add pRF outline
    ax.plot(ecc * np.cos(ang * np.pi / 180), ecc * np.sin(ang * np.pi / 180), color=color, markersize=3, marker='o')
    circ = plt.Circle((ecc * np.cos(ang * np.pi / 180), ecc * np.sin(ang * np.pi / 180)), siz,
                      color=color, fill=False, linewidth=2.5, alpha=1.0)
    ax.add_artist(circ)