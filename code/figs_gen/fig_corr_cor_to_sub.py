import os
import numpy as np
import nibabel as nib
import matplotlib
import matplotlib.pyplot as plt
from figure_funcs import ff
from matplotlib.colors import LightSource
from matplotlib import cm
from skimage import measure
matplotlib.use('Agg')
from scipy.io import loadmat
from scipy import ndimage
#from adjustText import adjust_text


def clean_data(data, tha):
    """ Cleans pulvino-cortical correlation data as it is loaded in

    Params:
        data (np.ndarray): Array of data of shape (182, 218, 182)
        tha (np.ndarray): Array of ROI labels indicating whether a given voxel is included in the posterior thalamus, of
            shape (182, 218, 182)
    Returns:
        data (np.ndarray): Array of data of shape (182, 218, 182). Data has been half-wave rectified, nans have been
            removed, and data outside the posterior thalamus has been eliminated
    """
    # Half-wave rectify
    data[data <= 0] = 0
    # Eliminate nans
    data[np.isnan(data)] = 0
    # Zero out data outside of the thalamus
    data[tha == 0] = 0
    # Return
    return data


def plot_correlation_contours_with_com(ax, data, roi, labels, colors, xlims, ylims, view='sagittal',
                                       contour=True, max=True, project=False, max_mode='max',
                                       project_side_y='left', project_side_x='left',
                                       annotate=False, annotation_labels=None, annotation_angles=None,
                                       hemi_constrain='lh'):
    # Define a function to handle checking the perimeter of different outlines
    def calc_perimeter(outline):
        perimeter = 0
        for ii in range(len(outline)):
            if ii == 0:
                continue
            else:
                perimeter += np.sqrt((outline[ii, 1]-outline[ii-1, 1])**2 + (outline[ii, 0]-outline[ii-1, 0])**2)
        return perimeter
    # Plot ROI overlay
    for selected_roi, linestyle in zip([1, 2, 3, 4], ['solid', 'dotted', 'dashed', 'dashdot']):
        ff.plot_mip_roi_overlay(ax, roi, view, linestyle=linestyle, color='gray', selected_rois=[selected_roi])

    # Loop through labels and handle contour
    for label, color in zip(labels, colors):
        # Select data we want
        selection = np.squeeze(data[label, :, :, :])
        # Constrain search to certain hemisphere
        if hemi_constrain == 'lh':
            selection[92:, :, :] = 0.0
        else:
            selection[0:92, :, :] = 0.0
        if view == 'sagittal':
            temp = np.max(selection, axis=0)
        elif view == 'coronal':
            temp = np.max(selection, axis=1)
        elif view == 'axial':
            temp = np.max(selection, axis=2)
        # If contour, plot contour
        if contour:
            contour = measure.find_contours(temp, np.max(temp) * 0.5)  # find contour at half-max
            if len(contour) > 0:
                perims = [calc_perimeter(cont) for cont in contour]
                largest_contour = contour[np.argmax(perims)]
                ax.plot(largest_contour[:, 0], largest_contour[:, 1], linewidth=2, color=color)
        # If we're looking at a coronal or axial slice, zero out RH
        if view == 'coronal' or view == 'axial':
            temp[92:, :] = 0
        # If max, plot max/com
        if max:
            if max_mode == 'max':
                idx_max = np.unravel_index(np.argmax(temp), temp.shape)
            else:
                idx_max = ndimage.center_of_mass(temp)
            idx_max = idx_max + np.random.uniform(-0.25, 0.25, (2,))
            plt.scatter(idx_max[0], idx_max[1], marker='x', color=color, zorder=2*len(labels)-label)
            # If annotate, we label each point with its cortical ROI label
            if annotate:
                angle = annotation_angles[label]
                r = 5
                x = idx_max[0] + np.cos(angle) * r
                y = idx_max[1] + np.sin(angle) * r
                dx = idx_max[0] - x
                dy = idx_max[1] - y
                plt.arrow(x, y, dx, dy, length_includes_head=True, head_length=0.05, head_width=0.05, width=0.01,
                          color=color)
                plt.text(x, y, annotation_labels[label], size=4, color=color, bbox=dict(facecolor='white', edgecolor=color))
                #texts.append(plt.text(idx_max[0], idx_max[1], annotation_labels[label]))
            # If project, plot tick indicating position of the max along axis
            if project:
                # Plot x-axis projection first
                if project_side_x == 'left':
                    plt.plot([idx_max[0], idx_max[0]], [ylims[0]-1, ylims[0]+1], color=color)
                else:
                    plt.plot([idx_max[0], idx_max[0]], [ylims[1]-1, ylims[1]+1], color=color)
                # Plot y-axis projection next
                if project_side_y == 'left':
                    plt.plot([xlims[0]-1, xlims[0]+1], [idx_max[1], idx_max[1]], color=color)
                else:
                    plt.plot([xlims[1]-1, xlims[1]+1], [idx_max[1], idx_max[1]], color=color)
    # Handle limits
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)


def plot_fig_corr_cor_to_sub_ventral_stream_avg_group_contours():
    # Load anatomy and ROI labels
    T1 = ff.load_volume(volume='T1')
    ROI = ff.load_volume(volume='thalamus')
    THA = ff.load_volume(volume='postthalamus')
    T1 = np.mean(np.array(T1), axis=0)

    # Load in and store average data
    maps_mean = []
    for label in range(14):
        temp = []
        for subj in range(8):
            # Load data
            temp.append(clean_data(nib.load(os.path.join(ff.dir_data, 'subj0' + str(subj + 1), 'mni',
                                                         'corr_cor_to_sub' '_hemi_' + str(0 + 1) + '_label_' + str(
                                                             label + 1) + '_method_2.nii.gz')).get_fdata(),
                                   THA[subj]))
        maps_mean.append(np.mean(np.array(temp), axis=0))
    maps_mean = np.array(maps_mean)

    # Set up colors
    cortical_roi_labels = ['V1', 'V2', 'V3', 'hV4', 'OFA', 'FFA', 'aTL-f', 'EBA', 'FBA', 'PPA', 'VWFA', 'MT', 'IPS', 'FEF']

    # Plot
    for view, limx, limy, project_side in zip(['sagittal', 'axial', 'coronal'],
                                              [(107, 88), (60, 85), (63, 82)],
                                              [(60, 85), (88, 108), (60, 85)],
                                              ['left', 'right', 'left']):
        # Create plot
        fig, ax = plt.subplots(figsize=(3, 3))
        # Plot plot
        plot_correlation_contours_with_com(ax, maps_mean, ROI, labels=[0, 1, 2, 3, 4, 5],
                                           colors=ff.cmap_corr(np.linspace(0, 1, 6)), view=view, xlims=limx,
                                           ylims=limy, project=True, project_side_x=project_side,
                                           project_side_y=project_side)
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.gca().set_aspect('equal')
        plt.tight_layout(pad=0)
        plt.savefig(os.path.join('../figures', 'fig_corr_cor_to_sub_mip_group_avg_contours_ventral_stream_' + view + '.png'),
                    dpi=300)
        plt.close()

    a = np.array([[0, 1]])
    plt.figure(figsize=(1, 2.25))
    img = plt.imshow(a, cmap=matplotlib.colors.ListedColormap(list(matplotlib.cm.Set3([0, 1, 2, 3, 4, 5]))))
    plt.gca().set_visible(False)
    cax = plt.axes([0.2, 0.1, 0.2, 0.8])
    cbar = plt.colorbar(cax=cax, ticks=np.linspace(1 / 6 / 2, 1 - 1 / 6 / 2, 6))
    cbar.set_ticklabels([x for x, ii in zip(cortical_roi_labels, range(14)) if ii in [0, 1, 2, 3, 4, 5]])


def plot_fig_corr_cor_to_sub_ventral_stream_avg_group_contours_ap_color_code():
    # Load anatomy and ROI labels
    ROI = ff.load_volume(volume='thalamus')
    THA = ff.load_volume(volume='postthalamus')

    # Load anatomical position data
    anatpos = loadmat(os.path.join(ff.dir_data, 'group', 'cortical_roi_average_positions.mat'))['results']

    # Load in and store average data
    maps_mean = []
    for label in range(14):
        temp = []
        for subj in range(8):
            # Load data
            temp.append(clean_data(nib.load(os.path.join(ff.dir_data, 'subj0' + str(subj + 1), 'mni',
                                                         'corr_cor_to_sub' '_hemi_' + str(0 + 1) + '_label_' + str(
                                                             label + 1) + '_method_2.nii.gz')).get_fdata(),
                                   THA[subj]))
        maps_mean.append(np.mean(np.array(temp), axis=0))
    maps_mean = np.array(maps_mean)

    for anat_axis in range(3):
        # Set up colors
        ap = np.mean(anatpos[anat_axis, 0, :, :], axis=0)
        ap = ap - np.nanmin(ap)
        ap = ap / np.nanmax(ap)
        label_colors = list(matplotlib.cm.plasma(ap))
        cortical_roi_labels = ['V1', 'V2', 'V3', 'hV4', 'OFA', 'FFA', 'aTL-f', 'EBA', 'FBA', 'PPA', 'VWFA', 'MT', 'IPS', 'FEF']

        # Plot
        for view, limx, limy, project_side in zip(['sagittal', 'axial', 'coronal'],
                                                  [(88, 108), (60, 85), (60, 85)],
                                                  [(60, 85), (88, 108), (60, 85)],
                                                  ['left', 'right', 'left']):
            # Create plot
            fig, ax = plt.subplots(figsize=(3, 3))
            # Plot plot
            plot_correlation_contours_with_com(ax, maps_mean, ROI, labels=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                                               colors=label_colors, view=view, xlims=limx, ylims=limy,
                                               contour=False, project=True, max=True, max_mode='com',
                                               project_side_x=project_side, project_side_y=project_side,
                                               annotate=True, annotation_labels=cortical_roi_labels,
                                               annotation_angles=np.linspace(-np.pi/2, 3/4*np.pi, 14))
            plt.gca().get_xaxis().set_visible(False)
            plt.gca().get_yaxis().set_visible(False)
            plt.gca().set_aspect('equal')
            plt.tight_layout(pad=0)
            plt.savefig(os.path.join('../figures', 'fig_corr_cor_to_sub_mip_group_avg_contours_ventral_stream_' +
                                     str(anat_axis) + '_coloring_' + view + '.png'),
                        dpi=400)
            plt.close()


def plot_fig_corr_cor_to_sub_individual_subject_consistency_maps():
    # Load anatomy and ROI labels
    T1 = ff.load_volume(volume='T1')
    ROI = ff.load_volume(volume='thalamus')
    THA = ff.load_volume(volume='postthalamus')
    T1 = np.mean(np.array(T1), axis=0)

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

    # Create plot
    plt.figure(figsize=(8, 2.0))
    view = 'coronal'
    idx_slice = 95
    lims_x = (55, 135)
    lims_y = (35, 115)
    for idx_map, map in enumerate(maps):
        for idx_subj in range(8):
            plt.subplot(2, 8, idx_subj + 1 + idx_map*8)
            ff.plot_slice_with_overlay(T1, map[idx_subj], idx_slice, view=view, mask=map[idx_subj] <= 0, lims_x=lims_x,
                                       lims_y=lims_y, cmap2='hot', clim2=(0, 0.1))
            plt.gca().set_aspect('equal')
            plt.tight_layout(pad=0.1, h_pad=0.1, w_pad=0.1)
        # Plot consistency map
        #plt.subplot(2, 9, 8 + 1 + idx_map * 9)
        #prob_map = np.sum(np.array(map) > 0.01, axis=0)
        #ff.plot_slice_with_overlay(T1, prob_map, idx_slice, view=view, mask=prob_map <= 0, lims_x=lims_x,
        #                           lims_y=lims_y, cmap2='plasma', clim2=(0, 8))
        #plt.gca().set_aspect('equal')
        #plt.tight_layout(pad=0)
    plt.savefig(os.path.join('../figures', 'fig_corr_cor_to_sub_individual_consistency_maps.png'), dpi=500,
                interpolation='none')
    plt.close()


def plot_fig_corr_cor_to_sub_group_consistency_maps():
    # Load anatomy and ROI labels
    T1 = ff.load_volume(volume='T1')
    ROI = ff.load_volume(volume='thalamus')
    THA = ff.load_volume(volume='postthalamus')
    T1 = np.mean(np.array(T1), axis=0)
    # Load in all data and store it in a list in order
    maps = ['contrastNEW', 'bodyauto']
    R2_contrast = np.mean(np.array(ff.load_volume(volume='contrastNEW_R2')), axis=0)
    R2_bodies = np.mean(np.array(ff.load_volume(volume='bodyauto_R2')), axis=0)
    R2_faces = np.mean(np.array(ff.load_volume(volume='faceauto_R2')), axis=0)
    THA_group = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()
    R2_contrast[THA_group == 0] = 0
    R2_bodies[THA_group == 0] = 0

    # Indicate which maps you want to load
    #['V1', 'V2', 'V3', 'hV4', 'OFA', 'FFA', 'aTL-f', 'EBA', 'FBA', 'PPA', 'VWFA', 'MT', 'IPS', 'FEF']
    maps_to_load = [0, 5, 8]

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

    # Create plot
    views = ['coronal', 'sagittal']
    idxs_slice = [95, 110]
    lims_x = [(55, 135), (120, 70)]
    lims_y = [(35, 115), (50, 100)]
    for view, idx_slice, lim_x, lim_y in zip(views, idxs_slice, lims_x, lims_y):
        for idx_map, map in enumerate(maps):
            plt.figure(figsize=(2.0, 2.0))
            # Plot consistency map
            prob_map = np.sum(np.array(map) > 0.01, axis=0)
            ff.plot_slice_with_overlay(T1, prob_map, idx_slice, view=view, mask=prob_map <= 0, lims_x=lim_x,
                                    lims_y=lim_y, cmap2=ff.cmap_cons, clim2=(0, 8))
            # Calculate maxima of R2_contrast and R2_bodies in this view
            def plot_maxima(vol, hemi, idx_c):
                vol = np.copy(vol)
                if hemi == 'lh':
                    vol[92:, :, :] = 0
                else:
                    vol[0:92, :, :] = 0
                if view == 'coronal':
                    idx_max = np.unravel_index(np.nanargmax(vol[:, idx_slice, :]), (182, 182))
                else:
                    idx_max = np.unravel_index(np.nanargmax(vol[idx_slice, :, :]), (218, 182))
                plt.scatter(*idx_max, s=30, color=ff.cmap_feat(idx_c))
            # Plot maxima
            if view == 'coronal':
                plot_maxima(R2_contrast, 'rh', 0)
                plot_maxima(R2_bodies, 'rh', 3)
                plot_maxima(R2_faces, 'rh', 4)
                plot_maxima(R2_contrast, 'lh', 0)
                plot_maxima(R2_bodies, 'lh', 3)
                plot_maxima(R2_faces, 'lh', 4)
            else:
                plot_maxima(R2_contrast, 'rh', 0)
                plot_maxima(R2_bodies, 'rh', 3)
                plot_maxima(R2_faces, 'rh', 4)
            # Plot slice indicator lines
            if view == 'coronal':
                plt.plot([110, 110], lim_y, color='dimgray', linestyle='dashed', linewidth=1.0)
            else:
                plt.plot([95, 95], lim_y, color='dimgray', linestyle='dashed', linewidth=1.0)
            plt.gca().set_aspect('equal')
            plt.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.0)
            plt.savefig(os.path.join('../figures', 'fig_corr_cor_to_sub_group_consistency_maps_' + view + '_' + str(idx_map) + '.png'), dpi=500,
                        interpolation='none')
    plt.close('all')

