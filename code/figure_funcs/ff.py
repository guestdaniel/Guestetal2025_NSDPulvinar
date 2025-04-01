import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from skimage import measure
import nibabel as nib
import os

# Handle some general sizing across figures
SMALL_SIZE = 18
MEDIUM_SIZE = 22
BIGGER_SIZE = 20
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Set up several colormaps
def create_cmap_arcaro():
    """ Creates a colormap similar to that from Arcaro et al. (2015)

    Returns:
        cmap_arcaro (ListedColormap): ListedColormap object containing the colormap
    """
    cmap_arcaro = np.array([[255, 0, 0],
                            [255, 68, 0],
                            [255, 153, 0],
                            [255, 255, 0],
                            [0, 255, 0],
                            [50, 205, 50],
                            [0, 255, 255],
                            [0, 153, 255],
                            [0, 105, 255],
                            [0, 68, 255]])
    cmap_arcaro = cmap_arcaro / 255
    cmap_arcaro = np.append(cmap_arcaro, np.ones((10, 1)), 1)
    cmap_arcaro = np.flip(cmap_arcaro, axis=0)
    cmap_arcaro = matplotlib.colors.ListedColormap(cmap_arcaro)
    return cmap_arcaro


cmap_angle = create_cmap_arcaro()                                    # colormap for all angle maps
cmap_rsqr = matplotlib.cm.get_cmap('hot')                            # colormap for variance explained / R^2 maps
cmap_size = matplotlib.cm.get_cmap('plasma_r')                       # colormap for size and eccentricity maps
cmap_lat  = matplotlib.cm.get_cmap('RdBu_r')                         # colormap for lateralization maps
cmap_cons = matplotlib.cm.get_cmap('winter')                         # colormap for consistency maps
#cmap_thal = matplotlib.colors.ListedColormap(matplotlib.cm.get_cmap('Set1')(np.linspace(0, 0.85, 5)))
                                                                     # colormap for hand-drawn ROI labels
cmap_thal = matplotlib.colors.ListedColormap([[0.0, 0.5, 0.4, 1.0],
                                             [0.3, 0.8, 1.0, 1.0],
                                             [0.2, 0.3, 0.7, 1.0],
                                             [0.4, 0.2, 0.7, 1.0],
                                             [0.7, 0.1, 0.6, 1.0]])                                                                     
#cmap_feat = matplotlib.colors.ListedColormap(matplotlib.cm.Accent(np.linspace(0, 1, 7)))
                                                                     # colormap for feature sets
# cmap_feat = matplotlib.colors.ListedColormap([[0.8, 0.2, 0.6, 1.0],
#                                               [0.5, 0.1, 0.3, 1.0],
#                                               [0.8, 0.1, 0.2, 1.0],
#                                               [0.9, 0.4, 0.0, 1.0],
#                                               [0.9, 0.8, 0.0, 1.0],
#                                               [0.6, 0.8, 0.2, 1.0],
#                                               [0.1, 0.6, 0.2, 1.0]])

cmap_feat = matplotlib.colors.ListedColormap([[0.8, 0.2, 0.6, 1.0],
                                              [0.8, 0.1, 0.2, 1.0],
                                              [0.9, 0.4, 0.0, 1.0],
                                              [0.9, 0.8, 0.0, 1.0],
                                              [0.1, 0.6, 0.2, 1.0]])


cmap_feat2 = matplotlib.colors.ListedColormap([[27/255, 158/255, 119/255, 1.0],
                                              [217/255, 95/255, 2/255, 1.0],
                                              [117/255, 112/255, 179/255, 1.0],
                                              [231/255, 41/255, 138/255, 1.0],
                                              [230/255, 171/255, 2/255, 1.0]])

cmap_corr = matplotlib.colors.ListedColormap([[0.0, 0.6, 0.5, 1.0],
                                              [0.0, 0.9, 0.5, 1.0],
                                              [0.4, 0.8, 0.3, 1.0],
                                              [0.9, 0.8, 0.0, 1.0],
                                              [0.0, 0.6, 0.8, 1.0],
                                              [0.1, 0.4, 0.8, 1.0]])

# Configure directories
dir_data = 'C:\\Users\\daniel\\Desktop\\Guestetal2025_NSDPulvinar\\data\\prepared'
dir_plots = 'C:\\Users\\daniel\\Desktop\\Guestetal2025_NSDPulvinar\\figures'

def load_volume(subjs=None, space='mni', volume='T1'):
    """ Load NIFTI data of multiple subjects from Guestetal2021_data

    Args:
        subjs (list): List of integers indicating which subjects' data to load. If None, defaults to all eight subjects.
            The contents of the list are not explicitly checked, and strange inputs will raised equally strange errors.
        space (string): a space to load from, either 'mni' or 'func1mm'.
        volume (string): name of datafile to load (e.g., 'T1')
    """
    # Set default subjs
    if subjs is None:
        subjs = [1, 2, 3, 4, 5, 6, 7, 8]
    elif not isinstance(subjs, list):
        raise TypeError('subject should be a list of integers!')
    # Handle space
    if space not in ['mni', 'func1mm']:
        raise ValueError('space must be either mni or func1mm')
    # Loop over subjs and load data
    storage = []
    for subj in subjs:
        storage.append(nib.load(os.path.join(dir_data, 'subj0' + str(subj), space, volume + '.nii.gz')).get_fdata())
    # Return
    return storage


def calc_areas(angle, size, ecc):
    # Convert angle to radians
    angle = np.pi/180 * angle
    # Calculate x, y, and radius
    x = np.cos(angle) * ecc
    y = np.sin(angle) * ecc
    radius = size*2
    # Calculate r, the distance from the center to the midline
    r = np.abs(x)
    # Calculate h, the distance from the edge of the circle in towards the midline
    h = radius - np.abs(r)
    # Calculate total area
    A = 2*np.pi*radius**2
    # If h is less than or equal to zero, that means the circle is entirely on one side... return now
    if h <= 0:
        if x >= 0:
            return 1
        else:
            return -1
    # Calculate R, the length of the line going from the midpoint to the intersection of the midline and circle
    R = h + r
    # Calculate a, the length of the segment
    a = 2 * np.sqrt(R**2 - r**2)
    # Calculate theta, the angle of the triangle formed by R and a
    theta = 2 * np.arctan2(a, 2*r)
    # Calculate the area of the sector and the area of the triangle
    A_sector = A * theta / (2 * np.pi)
    A_triangle = a * r
    # Calculate final area of the segment
    A_segment = A_sector - A_triangle
    if x >= 0:
        return ((A-A_segment)/A - 1/2)*2
    else:
        return (A_segment/A - 1/2)*2

    # plt.figure()
    # plt.gca().add_patch(plt.Circle((x, y), radius, fill=False))
    # plt.plot([0, 0], [-10, 10], color='black')
    # plt.plot([-10, 10], [0, 0], color='black')
    # plt.gca().set_aspect('equal')
    # plt.plot([x, x+r], [y, y], color='red')
    # plt.plot([x+r, x+r+h], [y, y], color='green')
    # plt.plot([x, x+r], [y, y+np.sin(theta/2)*R], color='orange')
    # plt.plot([x, x+r], [y, y-np.sin(theta/2)*R], color='orange')
    # plt.show()


def calc_lat(angle, size, ecc):
    """ Calculates a laterality metric for a set of angle, size, and eccentricity RF measurements

    Args:
        angle (ndarray): angle measurements
        size (ndarray): size measurements
        ecc (ndarray): eccentricity measurements

    Returns:
        lat (ndarray): laterality metric
    """
    lat = np.zeros_like(angle)
    for ii in range(angle.shape[0]):
        for jj in range(angle.shape[1]):
            for kk in range(angle.shape[2]):
                lat[ii, jj, kk] = calc_areas(angle[ii, jj, kk], size[ii, jj, kk], ecc[ii, jj, kk])
    return lat


def plot_prf_parameter_sequence(t1, r2, an, an_orig, ec, sz, thal, roi, slices=None):
    """ Plots variance explained, angle maps, eccentricity maps, and laterality maps on T1-weighted anatomy

    Args:
        t1 (ndarray): T1-weighted anatomy, one subject or average across subjects, full volume
        r2 (ndarray): pRF variance explained, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 180, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 360, ...
        ec (ndarray): pRF eccentricity maps, one subject or average across subjects, full volume
        sz (ndarray): pRF size maps, one subject or average across subjects, full volume
        thal (ndarray): ROI label indicating the extent of the thalamus
        roi (ndarray): ROI labels indicating different subcortical structures
        slices (list): list of slices to plot
    """
    # Construct figure
    if slices is None:
        slices = [100, 97, 95]
    n_slice = len(slices)
    cutoff = 0.1
    plt.figure(figsize=(10.5, 12))
    lat = calc_lat(an_orig, sz, ec)
    for idx_slice, slice in enumerate(slices):
        # R2
        plt.subplot(4, n_slice, idx_slice + 1)
        plot_slice_with_overlay(t1, r2, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_rsqr, clim2=(0, 1))
        plot_roi_overlay(roi, 'coronal', slice)
        # Angle
        plt.subplot(4, n_slice, idx_slice + 1 + n_slice * 1)
        plot_slice_with_overlay(t1, an, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_angle, clim2=(0, 180))
        # Eccentricity
        plt.subplot(4, n_slice, idx_slice + 1 + n_slice * 2)
        plot_slice_with_overlay(t1, ec, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_size, clim2=(0, 6))
        # Laterality
        plt.subplot(4, n_slice, idx_slice + 1 + n_slice * 3)
        plot_slice_with_overlay(t1, lat, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_lat, clim2=(-1, 1))
        # RF size
        #plt.subplot(4, n_slice, idx_slice + 1 + n_slice * 3)
        #plot_slice_with_overlay(t1, sz, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
        #                        (0 + 55, 182 - 80), cmap2=cmap_size, clim2=(0, 6))
    plt.tight_layout(pad=0.0, h_pad=0.1, w_pad=0.1)


def plotsave_prf_parameter_sequence_reduced(t1, r2, an, an_orig, ec, sz, thal, roi, savename, slices=None, outline_width=1.0):
    """ Plots variance explained, angle maps, and laterality maps on T1-weighted anatomy

    Args:
        t1 (ndarray): T1-weighted anatomy, one subject or average across subjects, full volume
        r2 (ndarray): pRF variance explained, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 180, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 360, ...
        ec (ndarray): pRF eccentricity maps, one subject or average across subjects, full volume
        sz (ndarray): pRF size maps, one subject or average across subjects, full volume
        thal (ndarray): ROI label indicating the extent of the thalamus
        roi (ndarray): ROI labels indicating different subcortical structures
        slices (list): list of slices to plot
    """
    # Construct figure
    if slices is None:
        slices = [100, 97, 95]
    n_slice = len(slices)
    cutoff = 0.1
    lat = calc_lat(an_orig, sz, ec)
    # Plot R2
    plt.figure(figsize=(7, 12))
    for idx_slice, slice in enumerate(slices):
        # R2
        plt.subplot(n_slice, 1, idx_slice + 1)
        plot_slice_with_overlay(t1, r2, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_rsqr, clim2=(0, 1))
        plot_roi_overlay(roi, 'coronal', slice, outline_width)
    plt.tight_layout()
    plt.savefig(savename + 'R2.png')
    # Plot Angle
    plt.figure(figsize=(7, 12))
    for idx_slice, slice in enumerate(slices):
        plt.subplot(n_slice, 1, idx_slice + 1)
        plot_slice_with_overlay(t1, an, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_angle, clim2=(0, 180))
    plt.tight_layout()
    plt.savefig(savename + 'angle.png')
    # Plot laterality
    plt.figure(figsize=(7, 12))
    for idx_slice, slice in enumerate(slices):
        plt.subplot(n_slice, 1, idx_slice + 1)
        plot_slice_with_overlay(t1, lat, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_lat, clim2=(-1, 1))
    plt.tight_layout()
    plt.savefig(savename + 'laterality.png')
    # Plot eccentricity
    plt.figure(figsize=(7, 12))
    for idx_slice, slice in enumerate(slices):
        plt.subplot(n_slice, 1, idx_slice + 1)
        plot_slice_with_overlay(t1, ec, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_size, clim2=(0, 6))
    plt.tight_layout()
    plt.savefig(savename + 'eccentricity.png')



def plotsave_prf_parameter_sequence_reduced_arcaro_custom(t1, r2, an, an_orig, ec, sz, thal, roi, savename, slices=None, outline_width=1.0, cutoff=0.1):
    """ Plots variance explained, angle maps, and laterality maps on T1-weighted anatomy

    Args:
        t1 (ndarray): T1-weighted anatomy, one subject or average across subjects, full volume
        r2 (ndarray): pRF variance explained, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 180, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 360, ...
        ec (ndarray): pRF eccentricity maps, one subject or average across subjects, full volume
        sz (ndarray): pRF size maps, one subject or average across subjects, full volume
        thal (ndarray): ROI label indicating the extent of the thalamus
        roi (ndarray): ROI labels indicating different subcortical structures
        slices (list): list of slices to plot
    """
    # Construct figure
    if slices is None:
        slices = [100, 97, 95]
    n_slice = len(slices)
    lat = calc_lat(an_orig, sz, ec)
    # Plot R2
    plt.figure(figsize=(7, 12))
    for idx_slice, slice in enumerate(slices):
        # R2
        plt.subplot(n_slice, 1, idx_slice + 1)
        plot_slice_with_overlay(t1, r2, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_rsqr, clim2=(0, 1))
        plot_roi_overlay(roi, 'coronal', slice, outline_width)
    plt.tight_layout()
    plt.savefig(savename + 'R2.png')


def plot_mip_roi_overlay(ax, ROI, view, color=None, linestyle='solid', linewidth=2, selected_rois=None):
    # Handle which ROIs to plot
    if selected_rois is None:
        selected_rois = [1, 2, 3, 4, 5]
    # Extract and overlay ROIs
    for roi_val in selected_rois:
        # Look only for one ROI val
        roi_present = np.array(ROI) == roi_val
        # Average over subjects
        roi_present = np.mean(roi_present, axis=0) > 2/8
        # Branch based on requested view
        if view == 'coronal':
            roi_present = np.max(roi_present, axis=1)
        elif view == 'sagittal':
            roi_present = np.max(roi_present, axis=0)
        elif view == 'axial':
            roi_present = np.max(roi_present, axis=2)
        contour = measure.find_contours(roi_present, 0.8)
        if len(contour) > 0:
            for cont in contour:
                if color is None:
                    ax.plot(cont[:, 0], cont[:, 1], linewidth=linewidth, linestyle=linestyle,
                             color=cmap_thal(np.linspace(0, 1, 5))[roi_val-1])
                elif isinstance(color, list):
                    ax.plot(cont[:, 0], cont[:, 1], linewidth=linewidth, linestyle=linestyle, color=color[roi_val-1])
                else:
                    ax.plot(cont[:, 0], cont[:, 1], linewidth=linewidth, linestyle=linestyle, color=color)


def plot_roi_overlay(ROI, view, idx, outline_width=1.0):
    # Extract and overlay ROIs
    for roi_val in [1, 2, 3, 4, 5]:
        # Look only for one ROI val
        roi_present = np.array(ROI) == roi_val
        # Average over subjects
        roi_present = np.mean(roi_present, axis=0) > 2/8
        # Branch based on requested view
        if view == 'coronal':
            roi_present = roi_present[:, idx, :]
        elif view == 'sagittal':
            roi_present = roi_present[idx, :, :]
        elif view == 'axial':
            roi_present = roi_present[:, :, idx]
        contour = measure.find_contours(roi_present, 0.8)
        if len(contour) > 0:
            for cont in contour:
                 plt.plot(cont[:, 0], cont[:, 1], linewidth=outline_width,
                          color=cmap_thal(np.linspace(0, 1, 5))[roi_val-1])


def plot_slice(img1, idx, view, lims_x=(0, 200), lims_y=(0, 200), cmap1='gray', clim1=None, marker=None,
               marker_color='black', marker_angle=0, marker_length=8):
    """ Plots slice of 3D brick of data.

    Args:
        img1 (ndarray): Volume of data, 3D
        idx (int): Index of the slice to plot
        view (str): either 'sagittal', 'coronal', or 'axial'
        lims_x (tuple): Limits in the plotted left-right axis
        lims_y (tuple): Limits in the plotted up-down axis
        cmap1 (str or np.array): Colormap to use for img1
        marker (tuple, list): a tuple or list of coordinates, an arrow is placed at these points
        marker_color (ndarray, string): some matplotlib-acceptable representation of a color to make the plotted arrow
        marker_angle (float): the angle of the marker in radians
        marker_length (float): the length of the marker from tip to tail
    """
    # Calculate marker start and stop positions
    x = marker_length * np.cos(marker_angle)
    y = marker_length * np.sin(marker_angle)
    # Plot marker
    if marker is not None:
        plt.arrow(marker[0]+x, marker[1]+y, -x, -y, color=marker_color, width=0.75, head_width=2, head_length=2,
                  length_includes_head=True)
    # Plot slice
    if view == 'sagittal':
        plt.imshow(np.squeeze(img1[idx, :, :]).T, cmap=cmap1)
    elif view == 'coronal':
        plt.imshow(np.squeeze(img1[:, idx, :]).T, cmap=cmap1)
    elif view == 'axial':
        plt.imshow(np.squeeze(img1[:, :, idx]).T, cmap=cmap1)
    if clim1 is not None:
        plt.clim(clim1)
    plt.xlim(lims_x)
    plt.ylim(lims_y)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)


def plot_slice_oo(ax, img1, idx, view, lims_x=(0, 200), lims_y=(0, 200), cmap1='gray', clim1=None, marker=None,
               marker_color='black', marker_angle=0, marker_length=8):
    """ Plots slice of 3D brick of data.

    Args:
        img1 (ndarray): Volume of data, 3D
        idx (int): Index of the slice to plot
        view (str): either 'sagittal', 'coronal', or 'axial'
        lims_x (tuple): Limits in the plotted left-right axis
        lims_y (tuple): Limits in the plotted up-down axis
        cmap1 (str or np.array): Colormap to use for img1
        marker (tuple, list): a tuple or list of coordinates, an arrow is placed at these points
        marker_color (ndarray, string): some matplotlib-acceptable representation of a color to make the plotted arrow
        marker_angle (float): the angle of the marker in radians
        marker_length (float): the length of the marker from tip to tail
    """
    # Calculate marker start and stop positions
    x = marker_length * np.cos(marker_angle)
    y = marker_length * np.sin(marker_angle)
    # Plot marker
    if marker is not None:
        ax.arrow(marker[0]+x, marker[1]+y, -x, -y, color=marker_color, width=0.75, head_width=2, head_length=2,
                  length_includes_head=True)
    # Plot slice
    if view == 'sagittal':
        ax.imshow(np.squeeze(img1[idx, :, :]).T, cmap=cmap1)
    elif view == 'coronal':
        ax.imshow(np.squeeze(img1[:, idx, :]).T, cmap=cmap1)
    elif view == 'axial':
        ax.imshow(np.squeeze(img1[:, :, idx]).T, cmap=cmap1)
    if clim1 is not None:
        ax.set_clim(clim1)
    ax.set_xlim(lims_x)
    ax.set_ylim(lims_y)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)


def plot_slice_token_markers(img1, idx, view, lims_x=(0, 200), lims_y=(0, 200), cmap1='gray', clim1=None, marker=None,
                             marker_color='black', marker_size=5):
    """ Plots slice of 3D brick of data.

    Args:
        img1 (ndarray): Volume of data, 3D
        idx (int): Index of the slice to plot
        view (str): either 'sagittal', 'coronal', or 'axial'
        lims_x (tuple): Limits in the plotted left-right axis
        lims_y (tuple): Limits in the plotted up-down axis
        cmap1 (str or np.array): Colormap to use for img1
        marker (tuple, list): a tuple or list of coordinates, an arrow is placed at these points
        marker_color (ndarray, string): some matplotlib-acceptable representation of a color to make the plotted arrow
        marker_angle (float): the angle of the marker in radians
        marker_length (float): the length of the marker from tip to tail
    """
    # Plot marker
    if marker is not None:
        if isinstance(marker, list):
            for _m, _c in zip(marker, marker_color):
                plt.plot(_m[0], _m[1], color=_c, marker='o', markersize=5)
        else:
            plt.plot(marker[0], marker[1], color=marker_color, marker='o', markersize=5)
    # Plot slice
    if view == 'sagittal':
        plt.imshow(np.squeeze(img1[idx, :, :]).T, cmap=cmap1)
    elif view == 'coronal':
        plt.imshow(np.squeeze(img1[:, idx, :]).T, cmap=cmap1)
    elif view == 'axial':
        plt.imshow(np.squeeze(img1[:, :, idx]).T, cmap=cmap1)
    if clim1 is not None:
        plt.clim(clim1)
    plt.xlim(lims_x)
    plt.ylim(lims_y)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)


def plot_slice_with_overlay(img1, img2, idx, view, mask, lims_x=(0, 200), lims_y=(0, 200), cmap1='gray', cmap2=cmap_rsqr,
                            clim1=None, clim2=None, marker=None, marker_color='black', marker_angle=0, marker_length=8):
    """ Plots slice of 3D brick of data with a corresponding slice of an overlay.

    Args:
        img1 (ndarray): Volume of data, 3D
        img2 (ndarray): Volume of data, 3D, plotted as overlay
        idx (int): Index of the slice to plot
        view (str): either 'sagittal', 'coronal', or 'axial'
        lims_x (tuple): Limits in the plotted left-right axis
        lims_y (tuple): Limits in the plotted up-down axis
        cmap1 (str or np.array): Colormap to use for img1
        cmap2 (str or np.array): Colormap to use for img2
        clim1 (tuple, list): limits of the first colormap. If None, default limits are applied.
        clim2 (tuple, list): limits of the second colormap. If None, default limits are applied.
        marker (tuple, list): a tuple or list of coordinates, an arrow is placed at these points
        marker_color (ndarray, string): some matplotlib-acceptable representation of a color to make the plotted arrow
        marker_angle (float): the angle of the marker in radians
        marker_length (float): the length of the marker from tip to tail
    """
    # Mask overlay image where mask indicates
    img2 = np.ma.masked_where(mask, img2)
    # Plot underlay
    plot_slice(img1, idx, view, lims_x, lims_y, cmap1, clim1, marker, marker_color, marker_angle, marker_length)
    # Plot overlay
    plot_slice(img2, idx, view, lims_x, lims_y, cmap2, clim2)


def plot_mip(img1, view, lims_x=(0, 200), lims_y=(0, 200), cmap1='gray', clim1=None, marker=None,
            marker_color='black', marker_angle=0, marker_length=8):
    """ Plots maximum intensity projection of 3D volume of data

    Args:
        img1 (ndarray): Volume of data, 3D
        view (str): either 'sagittal', 'coronal', or 'axial'
        lims_x (tuple): Limits in the plotted left-right axis
        lims_y (tuple): Limits in the plotted up-down axis
        cmap1 (str or np.array): Colormap to use for img1
        marker (tuple, list): a tuple or list of coordinates, an arrow is placed at these points
        marker_color (ndarray, string): some matplotlib-acceptable representation of a color to make the plotted arrow
        marker_angle (float): the angle of the marker in radians
        marker_length (float): the length of the marker from tip to tail
    """
    # Calculate marker start and stop positions
    x = marker_length * np.cos(marker_angle)
    y = marker_length * np.sin(marker_angle)
    # Plot marker
    if marker is not None:
        plt.arrow(marker[0]+x, marker[1]+y, -x, -y, color=marker_color, width=0.75, head_width=2, head_length=2,
                  length_includes_head=True)
    # Plot slice
    if view == 'sagittal':
        plt.imshow(np.squeeze(np.nanmax(img1, axis=0)).T, cmap=cmap1)
    elif view == 'coronal':
        plt.imshow(np.squeeze(np.nanmax(img1, axis=1)).T, cmap=cmap1)
    elif view == 'axial':
        plt.imshow(np.squeeze(np.nanmax(img1, axis=2)).T, cmap=cmap1)
    if clim1 is not None:
        plt.clim(clim1)
    plt.xlim(lims_x)
    plt.ylim(lims_y)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)


def plot_rainbow_slice(img1, img2, idx_slice, mask, lims_x=(50, 132), lims_y=(40, 102)):
    '''
    Plots base layer of data and overlay layer of data with 3D color
    Args:
        img1 (np.array): Volume of data, plotted as base layer
        img2 (np.array): Volume of data, plotted as overlay
        idx_slice (int): Index of the coronal slice to plot
        mask (np.array): Mask of which voxels to *exclude* in overlay
        lims_x (tuple): Limits in the left-right axis
        lims_y (tuple): Limits in the anterior-posterior axis
        cmap1 (str or np.array): Colormap to use for img1
        cmap2 (str or np.array): Colormap to use for img2
        transp_mask (np.array): Mask of voxels to visually suppress in img2 data via gray transparency overlay
    '''
    img2 = np.concatenate([img2, np.expand_dims(255*np.logical_not(mask), axis=3)], axis=3)
    #mask = np.stack([mask, mask, mask], axis=3)  # stack mask to match the fact that the data has color dimension
    #img2 = np.ma.masked_where(mask, img2)
    plt.imshow(np.squeeze(img1[:, idx_slice, :]).T, cmap='gray')
    plt.imshow(np.transpose(np.squeeze(img2[:, idx_slice, :, :]), [1, 0, 2]))
    plt.axis('off')
    plt.xlim(lims_x)
    plt.ylim(lims_y)


def plot_visual_field_coverage(R2, AN, EC, SZ, include, cutoff=0.25):
    '''
    Plots visual fields in stimulus space
    Args:
        R2 (list): list of 3D array of size (182, 218, 182) with variance explained metric
        AN (list): list of 3D array of size (182, 218, 182) with angle parameter estimates
        EC (list): list of 3D array of size (182, 218, 182) with eccentricity parameter estimates
        SZ (list): list of 3D array of size (182, 218, 182) with size parameter estimates
        include (ndarray): 3D array of size (182, 218, 182) with boolean values indicating which voxels to plot
        cutoff (float): cutoff of which voxels to include based on R2
    Returns:
        None

    '''
    fig, ax = plt.subplots(8, 2, figsize=(3.75, 16))
    for subj in range(8):
        # Prepare variables to calculate coverage maps
        maps = []
        x = np.linspace(-8, 8, 200)
        y = np.linspace(-8, 8, 200)
        # Design plots
        for row in range(2):
            # Lay out basic markings on plot
            ax[subj][row].plot([-10, 10], [0, 0], 'r--')
            ax[subj][row].plot([0, 0], [-10, 10], 'r--')
            rect = plt.Rectangle(xy=(-4.2, -4.2), width=8.4, height=8.4, color='r', fill=False)
            ax[subj][row].add_artist(rect)
            # Set limits
            ax[subj][row].set_xlim((-8.4, 8.4))
            ax[subj][row].set_ylim((-8.4, 8.4))
            # Set ticks
            ax[subj][row].set_xticks([-8.4, -4.2, 0, 4.2, 8.4])
            ax[subj][row].set_xticks([-6.3, -2.1, 2.1, 6.3], minor=True)
            ax[subj][row].set_yticks([-8.4, -4.2, 0, 4.2, 8.4])
            ax[subj][row].set_yticks([-6.3, -2.1, 2.1, 6.3], minor=True)
            # Force square axes
            ax[subj][row].set(adjustable='box', aspect='equal')
        # Remove unnecessary x/y-axes
        if subj != 7:
            ax[subj][0].get_xaxis().set_visible(False)
            ax[subj][1].get_xaxis().set_visible(False)
        ax[subj][1].get_yaxis().set_visible(False)
        for hemi in range(2):
            # Set all R2 to zero where voxels are not in include
            R2[subj][include != 1] = 0
            # Loop through hemispheres and plot RFs
            if hemi == 0:  # left hemisphere
                R2_subset = R2[subj][0:93, :, :]
                AN_subset = AN[subj][0:93, :, :]
                EC_subset = EC[subj][0:93, :, :]
                SZ_subset = SZ[subj][0:93, :, :]
            else:          # right hemisphere
                R2_subset = R2[subj][93:, :, :]
                AN_subset = AN[subj][93:, :, :]
                EC_subset = EC[subj][93:, :, :]
                SZ_subset = SZ[subj][93:, :, :]
            n_vox = len(R2_subset[R2_subset > cutoff])
            # Plot fields
            for vox in range(n_vox):
                # Extract parameters
                ang = AN_subset[R2_subset > cutoff][vox]
                ecc = EC_subset[R2_subset > cutoff][vox]
                size = SZ_subset[R2_subset > cutoff][vox]
                # Plot center
                if hemi == 0:
                    ax[subj][0].plot(ecc*np.cos(ang*np.pi/180), ecc*np.sin(ang*np.pi/180), color=[0, 1, 0], markersize=3, marker='o')
                else:
                    ax[subj][0].plot(ecc*np.cos(ang*np.pi/180), ecc*np.sin(ang*np.pi/180), color=[0, 0, 1], markersize=3, marker='o')
                # Plot surround
                circ = plt.Circle((ecc*np.cos(ang*np.pi/180), ecc*np.sin(ang*np.pi/180)), size, color='k', fill=False, linewidth=0.5, alpha=0.1)
                ax[subj][0].add_artist(circ)
            # Calculate coverage maps
            X, Y = np.meshgrid(x, y)
            z = np.zeros((200, 200))
            for vox in range(n_vox):
                # Extract parameters
                ang = AN_subset[R2_subset > cutoff][vox]
                ecc = EC_subset[R2_subset > cutoff][vox]
                size = SZ_subset[R2_subset > cutoff][vox]
                in_circle = np.sqrt(
                     (X - ecc * np.cos(ang * np.pi / 180)) ** 2 + (Y - ecc * np.sin(ang * np.pi / 180)) ** 2) < size
                z = z + in_circle
            maps.append(z)
        # Combine coverage maps
        coverage_map_oppositional = -1*maps[0] + maps[1]
        rgb_lh = np.zeros((200, 200, 3))
        rgb_lh[:, :, 1] = maps[0] / max(np.max(maps[0]), np.max(maps[1]))
        rgb_rh = np.zeros((200, 200, 3))
        rgb_rh[:, :, 2] = maps[1] / max(np.max(maps[0]), np.max(maps[1]))
        coverage_map = rgb_lh + rgb_rh
        # Plot
        #ax[subj][1].pcolormesh(x, y, coverage_map/np.max(np.abs(coverage_map)), cmap='GnBu', clim=(-1, 1))
        ax[subj][1].imshow(coverage_map, extent=(-8.4, 8.4, -8.4, 8.4), origin='lower')
        #ax[subj][1].set_ylim((-8.4, 8.4))  # todo: investigate
        # Calculate LH peak
        #matcoord_peak_LH = np.unravel_index(np.argmax(-1*coverage_map), coverage_map.shape)
        #matcoord_peak_RH = np.unravel_index(np.argmax(coverage_map), coverage_map.shape)
        #ax[subj][1].plot(X[matcoord_peak_LH], Y[matcoord_peak_LH], marker='x', color='r')
        #ax[subj][1].plot(X[matcoord_peak_RH], Y[matcoord_peak_RH], marker='x', color='r')
    plt.tight_layout(w_pad=0)


def summarize_centrality_rf_coverage(R2, EC, include, cutoff=0.1):
    '''
    Companion function to `plot_visual_field_coverage_horizontal` below. Prints to console
    the proportions of each set of eccentricity estimates that fall within inclusion 
    criteria that are within vs without 3 degrees of central visual field.

    Args:
        R2 (list): list of 3D array of size (182, 218, 182) with variance explained metric
        EC (list): list of 3D array of size (182, 218, 182) with eccentricity parameter estimates
        include (array): 3D array of size (182, 218, 182) with boolean values indicating which voxels to plot
        cutoff (float): cutoff of which voxels to include based on R2
    '''
    for subj in range(9):
        # Set all R2 to zero where voxels are not in include
        R2[subj][include != 1] = 0
        # Temporarily reassign R2[subj] and EC[subj]
        R2_temp = R2[subj]
        EC_temp = EC[subj]
        # Identify how many voxels exceed cutoff
        n_vox = len(R2_temp[R2_temp > cutoff])
        # Identify how many voxels exceed cutoff and also are within 3 degrees of center
        n_within = len(R2_temp[np.logical_and(EC_temp <= 3.0, R2_temp > cutoff)])
        n_without = len(R2_temp[np.logical_and(EC_temp > 3.0, R2_temp > cutoff)])
        assert(n_within+n_without == n_vox)
        print("Subj " + str(subj) + ": " + str(n_within/n_vox*100) + "%")


def plot_visual_field_coverage_horizontal(R2, AN, EC, SZ, include, cutoff=0.25, max_vox=100):
    '''
    Plots visual fields in stimulus space
    Args:
        AN (list): list of 3D array of size (182, 218, 182) with angle parameter estimates
        EC (list): list of 3D array of size (182, 218, 182) with eccentricity parameter estimates
        SZ (list): list of 3D array of size (182, 218, 182) with size parameter estimates
        R2 (list): list of 3D array of size (182, 218, 182) with variance explained metric
        include (array): 3D array of size (182, 218, 182) with boolean values indicating which voxels to plot
        cutoff (float): cutoff of which voxels to include based on R2
    Returns:
        None

    '''
    fig, ax = plt.subplots(2, 9, figsize=(18, 4.875))  # make a figure with two rows (different maps) and 9 columns (different subjects and mean subject)
    for subj in range(9):
        # Prepare variables to calculate coverage maps
        maps = []
        x = np.linspace(-8, 8, 200)
        y = np.linspace(-8, 8, 200)
        # Design plots
        for row in range(2):
            # Lay out basic markings on plot
            ax[row][subj].plot([-10, 10], [0, 0], linestyle='dashed', color='gray')
            ax[row][subj].plot([0, 0], [-10, 10], linestyle='dashed', color='gray')
            rect = plt.Rectangle(xy=(-4.2, -4.2), width=8.4, height=8.4, color='gray', fill=False)
            ax[row][subj].add_artist(rect)
            # Set limits
            ax[row][subj].set_xlim((-8.4, 8.4))
            ax[row][subj].set_ylim((-8.4, 8.4))
            # Set ticks
            #ax[row][subj].set_xticks([-8.4, -4.2, 0, 4.2, 8.4])
            #ax[row][subj].set_xticklabels(['', '-4.2', '0', '4.2', ''])
            #ax[row][subj].set_xticklabels(['', '', '', '', ''])
            ax[row][subj].set_xticks([])
            #ax[row][subj].set_xticks([-6.3, -2.1, 2.1, 6.3], minor=True)
            #ax[row][subj].set_yticks([-8.4, -4.2, 0, 4.2, 8.4])
            #ax[row][subj].set_yticks([-6.3, -2.1, 2.1, 6.3], minor=True)
            #ax[row][subj].set_yticklabels(['', '', '', '', ''])
            ax[row][subj].set_yticks([])
            # Force square axes
            ax[row][subj].set(adjustable='box', aspect='equal')
        # Remove unnecessary x/y-axes
        if subj > 0:
            ax[0][subj].get_yaxis().set_visible(False)
            ax[1][subj].get_yaxis().set_visible(False)
        ax[0][subj].get_xaxis().set_visible(False)
        for hemi in range(2):
            # Set all R2 to zero where voxels are not in include
            R2[subj][include != 1] = 0
            # Loop through and plot RFs
            if hemi == 0:  # left hemisphere
                R2_subset = R2[subj][0:93, :, :]
                AN_subset = AN[subj][0:93, :, :]
                EC_subset = EC[subj][0:93, :, :]
                SZ_subset = SZ[subj][0:93, ::, :]
            else:          # right hemisphere
                R2_subset = R2[subj][93:, :, :]
                AN_subset = AN[subj][93:, :, :]
                EC_subset = EC[subj][93:, :, :]
                SZ_subset = SZ[subj][93:, :, :]
            n_vox = len(R2_subset[R2_subset > cutoff])
            R2vals = R2_subset[R2_subset > cutoff]   # grab vector of R2 values that exceed cutoff
            idxs_sort = np.flip(np.argsort(R2vals))  # grab vector to sort R2 values in descending order
            # Plot fields
            for vox in range(min(n_vox, max_vox)):
                # Extract parameters (note that we index into voxels sorted in descending order of R2)
                ang = AN_subset[R2_subset > cutoff][idxs_sort][vox]
                ecc = EC_subset[R2_subset > cutoff][idxs_sort][vox]
                size = SZ_subset[R2_subset > cutoff][idxs_sort][vox]
                # Plot center & surround
                if hemi == 0:
                    ax[0][subj].plot(ecc*np.cos(ang*np.pi/180), ecc*np.sin(ang*np.pi/180), color=[1, 0, 0], markersize=3, marker='o')
                    circ = plt.Circle((ecc * np.cos(ang * np.pi / 180), ecc * np.sin(ang * np.pi / 180)), size,
                                      color=[1, 0, 0], fill=False, linewidth=0.5, alpha=0.1)
                    ax[0][subj].add_artist(circ)
                else:
                    ax[0][subj].plot(ecc*np.cos(ang*np.pi/180), ecc*np.sin(ang*np.pi/180), color=[0, 0, 1], markersize=3, marker='o')
                    circ = plt.Circle((ecc * np.cos(ang * np.pi / 180), ecc * np.sin(ang * np.pi / 180)), size,
                                      color=[0, 0, 1], fill=False, linewidth=0.5, alpha=0.1)
                    ax[0][subj].add_artist(circ)
            # Calculate coverage maps
            X, Y = np.meshgrid(x, y)
            z = np.zeros((200, 200))
            for vox in range(min(n_vox, max_vox)):
                # Extract parameters (note that we index into voxels sorted in descending order of R2)
                ang = AN_subset[R2_subset > cutoff][idxs_sort][vox]
                ecc = EC_subset[R2_subset > cutoff][idxs_sort][vox]
                size = SZ_subset[R2_subset > cutoff][idxs_sort][vox]
                in_circle = np.sqrt(
                     (X - ecc * np.cos(ang * np.pi / 180)) ** 2 + (Y - ecc * np.sin(ang * np.pi / 180)) ** 2) < size
                z = z + in_circle
            maps.append(z)
        # Combine coverage maps
        coverage_map_oppositional = -1*maps[0] + maps[1]
        rgb_lh = np.zeros((200, 200, 3))
        rgb_lh[:, :, 0] = maps[0] / max(np.max(maps[0]), np.max(maps[1]))
        rgb_rh = np.zeros((200, 200, 3))
        rgb_rh[:, :, 2] = maps[1] / max(np.max(maps[0]), np.max(maps[1]))
        coverage_map = rgb_lh + rgb_rh
        # Plot
        #ax[1][subj].pcolormesh(x, y, coverage_map/np.max(np.abs(coverage_map)), cmap='GnBu', clim=(-1, 1))
        ax[1][subj].imshow(coverage_map, extent=(-8.4, 8.4, -8.4, 8.4), origin='lower')
        #ax[subj][1].set_ylim((-8.4, 8.4))  # todo: investigate
        # Calculate LH peak
        #matcoord_peak_LH = np.unravel_index(np.argmax(-1*coverage_map), coverage_map.shape)
        #matcoord_peak_RH = np.unravel_index(np.argmax(coverage_map), coverage_map.shape)
        #ax[subj][1].plot(X[matcoord_peak_LH], Y[matcoord_peak_LH], marker='x', color='r')
        #ax[subj][1].plot(X[matcoord_peak_RH], Y[matcoord_peak_RH], marker='x', color='r')
    plt.tight_layout(w_pad=0)