import os
os.chdir("code")
from figure_funcs import ff
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
import importlib

# CONTRAST
T1 = ff.load_volume(volume='T1')
R2 = ff.load_volume(volume='contrastNEW_R2')
R2_body = ff.load_volume(volume='bodyauto_R2')
AN = ff.load_volume(volume='contrastNEW_angle')
EC = ff.load_volume(volume='contrastNEW_eccentricity')
SZ = ff.load_volume(volume='contrastNEW_size')
ROI = ff.load_volume(volume='thalamus')
ROI = np.array(ROI)
T1 = np.mean(np.array(T1), axis=0)
R2 = np.median(np.array(R2), axis=0)
R2_body = np.median(np.array(R2_body), axis=0)
# AN_orig: (1) convert AN to radians; (2) turn into complex-valued vector representation;
#          (3) average; (4) extract angle in radians; (5) convert to degrees; (6) take real part of complex 
AN = np.abs(np.angle(np.nanmedian(np.exp(1j * (np.array(AN)*np.pi/180 + np.pi/2)), axis=0)))*180/np.pi
EC = np.median(np.array(EC), axis=0)
SZ = np.median(np.array(SZ), axis=0)
THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()

# # Angle analysis 3D (left hemisphere)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'coronal')
# plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_left_coronal.png'), bbox_inches='tight', pad_inches=0)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'axial')
# plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_left_axial.png'), bbox_inches='tight', pad_inches=0)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'sagittal')
# plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_left_sagittal.png'), bbox_inches='tight', pad_inches=0)

# # Angle analysis 3D (right h)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'coronal')
# plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_right_coronal.png'), bbox_inches='tight', pad_inches=0)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'axial')
# plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_right_axial.png'), bbox_inches='tight', pad_inches=0)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'sagittal')
# plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_right_sagittal.png'), bbox_inches='tight', pad_inches=0)

# Angle analysis summaries
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 20, 'left', 'coronal')
# plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_left_coronal_shuffled_maps.png'), bbox_inches='tight', pad_inches=0)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 500, 'left', 'coronal')
# plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_left_coronal_shuffled_histogram.png'), bbox_inches='tight', pad_inches=0)

# Angle analysis data figures left and right hemi
for wta in [True, False]:
    for r2 in [0.05, 0.1, 0.15, 0.2]:
        if wta:
            importlib.reload(ff); ff.char_prf_angle_3D_fits_vs_data(AN, R2, THA, ROI, R2_body, 100, r2, 'left');
            plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_predicted_maps_lh_r2=%s_wta=true.png' % r2), bbox_inches='tight', pad_inches=0)
            importlib.reload(ff); ff.char_prf_angle_3D_fits_vs_data(AN, R2, THA, ROI, R2_body, 100, r2, 'right');
            plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_predicted_maps_rh_r2=%s_wta=true.png' % r2), bbox_inches='tight', pad_inches=0)
        else:
            importlib.reload(ff); ff.char_prf_angle_3D_fits_vs_data(AN, R2, THA, ROI, None, 100, r2, 'left');
            plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_predicted_maps_lh_r2=%s_wta=false.png' % r2), bbox_inches='tight', pad_inches=0)
            importlib.reload(ff); ff.char_prf_angle_3D_fits_vs_data(AN, R2, THA, ROI, None, 100, r2, 'right');
            plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_predicted_maps_rh_r2=%s_wta=false.png' % r2), bbox_inches='tight', pad_inches=0)

# BODY
T1 = ff.load_volume(volume='T1')
R2 = ff.load_volume(volume='bodyauto_R2')
R2_contrast = ff.load_volume(volume='contrastNEW_R2')
AN = ff.load_volume(volume='bodyauto_angle')
EC = ff.load_volume(volume='bodyauto_eccentricity')
SZ = ff.load_volume(volume='bodyauto_size')
ROI = ff.load_volume(volume='thalamus')
ROI = np.array(ROI)
T1 = np.mean(np.array(T1), axis=0)
R2 = np.median(np.array(R2), axis=0)
R2_contrast = np.median(np.array(R2_contrast), axis=0)
AN = np.abs(np.angle(np.nanmedian(np.exp(1j * (np.array(AN)*np.pi/180 + np.pi/2)), axis=0)))*180/np.pi
EC = np.median(np.array(EC), axis=0)
SZ = np.median(np.array(SZ), axis=0)
THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()

# # Angle analysis in 3D
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'coronal')
# plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_left_coronal.png'), bbox_inches='tight', pad_inches=0)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'axial')
# plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_left_axial.png'), bbox_inches='tight', pad_inches=0)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'sagittal')
# plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_left_sagittal.png'), bbox_inches='tight', pad_inches=0)

# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'coronal')
# plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_right_coronal.png'), bbox_inches='tight', pad_inches=0)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'axial')
# plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_right_axial.png'), bbox_inches='tight', pad_inches=0)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'sagittal')
# plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_right_sagittal.png'), bbox_inches='tight', pad_inches=0)

# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 20, 'left', 'coronal')
# plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_left_coronal_shuffled_maps.png'), bbox_inches='tight', pad_inches=0)
# importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 500, 'left', 'coronal')
# plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_left_coronal_shuffled_histogram.png'), bbox_inches='tight', pad_inches=0)

# Angle analysis data figures left and right hemi
for wta in [True, False]:
    for r2 in [0.05, 0.1, 0.15, 0.2]:
        if wta:
            importlib.reload(ff); ff.char_prf_angle_3D_fits_vs_data(AN, R2, THA, ROI, R2_contrast, 100, r2, 'left');
            plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_predicted_maps_lh_r2=%s_wta=true.png' % r2), bbox_inches='tight', pad_inches=0)
            importlib.reload(ff); ff.char_prf_angle_3D_fits_vs_data(AN, R2, THA, ROI, R2_contrast, 100, r2, 'right');
            plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_predicted_maps_rh_r2=%s_wta=true.png' % r2), bbox_inches='tight', pad_inches=0)
        else:
            importlib.reload(ff); ff.char_prf_angle_3D_fits_vs_data(AN, R2, THA, ROI, None, 100, r2, 'left');
            plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_predicted_maps_lh_r2=%s_wta=false.png' % r2), bbox_inches='tight', pad_inches=0)
            importlib.reload(ff); ff.char_prf_angle_3D_fits_vs_data(AN, R2, THA, ROI, None, 100, r2, 'right');
            plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_predicted_maps_rh_r2=%s_wta=false.png' % r2), bbox_inches='tight', pad_inches=0)

# Correlation analyses
# import os
# os.chdir("code")
# from figure_funcs import ff
# import matplotlib.pyplot as plt
# import numpy as np
# import nibabel as nib
# import importlib

# from figs_gen import fig_corr_cor_to_sub
# importlib.reload(fig_corr_cor_to_sub)
# fig_corr_cor_to_sub.plot_fig_corr_cor_to_sub_ventral_stream_avg_glm_decomposition()
# plt.savefig(os.path.join('../figures', 'fig_corr_cor_to_sub_decomposition.png'), bbox_inches='tight', pad_inches=0)

