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
AN = ff.load_volume(volume='contrastNEW_angle')
EC = ff.load_volume(volume='contrastNEW_eccentricity')
SZ = ff.load_volume(volume='contrastNEW_size')
ROI = ff.load_volume(volume='thalamus')
ROI = np.array(ROI)
T1 = np.mean(np.array(T1), axis=0)
R2 = np.median(np.array(R2), axis=0)
# AN_orig: (1) convert AN to radians; (2) turn into complex-valued vector representation;
#          (3) average; (4) extract angle in radians; (5) convert to degrees; (6) take real part of complex 
AN = np.abs(np.angle(np.nanmedian(np.exp(1j * (np.array(AN)*np.pi/180 + np.pi/2)), axis=0)))*180/np.pi
EC = np.median(np.array(EC), axis=0)
SZ = np.median(np.array(SZ), axis=0)
THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()

# Angle analysis
# importlib.reload(ff); ff.char_prf_angle(AN, R2, THA, np.arange(102, 93, -1), 5, 'left', 'coronal'); plt.show() 
# importlib.reload(ff); ff.char_prf_angle(AN, R2, THA, np.arange(102, 93, -1), 5, 'right', 'coronal'); plt.show() 
# importlib.reload(ff); ff.char_prf_angle(AN, R2, THA, np.arange(62, 85, 1), 5, 'left', 'sagittal'); plt.show() 
# importlib.reload(ff); ff.char_prf_angle(AN, R2, THA, np.arange(182-62, 182-85, -1), 5, 'right', 'sagittal'); plt.show() 
# importlib.reload(ff); ff.char_prf_angle(AN, R2, THA, np.arange(65, 85, 1), 5, 'left', 'axial'); plt.show() 
# importlib.reload(ff); ff.char_prf_angle(AN, R2, THA, np.arange(65, 85, 1), 5, 'right', 'axial'); plt.show() 

# Angle analysis 3D (left hemisphere)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'coronal')
plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_left_coronal.png'), bbox_inches='tight', pad_inches=0)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'axial')
plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_left_axial.png'), bbox_inches='tight', pad_inches=0)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'sagittal')
plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_left_sagittal.png'), bbox_inches='tight', pad_inches=0)

# Angle analysis 3D (right h)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'coronal')
plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_right_coronal.png'), bbox_inches='tight', pad_inches=0)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'axial')
plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_right_axial.png'), bbox_inches='tight', pad_inches=0)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'sagittal')
plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_right_sagittal.png'), bbox_inches='tight', pad_inches=0)

importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 20, 'left', 'coronal')
plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_left_coronal_shuffled_maps.png'), bbox_inches='tight', pad_inches=0)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 500, 'left', 'coronal')
plt.savefig(os.path.join('../figures', 'fig_prf_contrast_char_angle_3D_left_coronal_shuffled_histogram.png'), bbox_inches='tight', pad_inches=0)

# BODY
T1 = ff.load_volume(volume='T1')
R2 = ff.load_volume(volume='bodyauto_R2')
AN = ff.load_volume(volume='bodyauto_angle')
EC = ff.load_volume(volume='bodyauto_eccentricity')
SZ = ff.load_volume(volume='bodyauto_size')
ROI = ff.load_volume(volume='thalamus')
ROI = np.array(ROI)
T1 = np.mean(np.array(T1), axis=0)
R2 = np.median(np.array(R2), axis=0)
AN = np.abs(np.angle(np.nanmedian(np.exp(1j * (np.array(AN)*np.pi/180 + np.pi/2)), axis=0)))*180/np.pi
EC = np.median(np.array(EC), axis=0)
SZ = np.median(np.array(SZ), axis=0)
THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()

# Angle analysis in 3D
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'coronal')
plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_left_coronal.png'), bbox_inches='tight', pad_inches=0)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'axial')
plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_left_axial.png'), bbox_inches='tight', pad_inches=0)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'left', 'sagittal')
plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_left_sagittal.png'), bbox_inches='tight', pad_inches=0)

importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'coronal')
plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_right_coronal.png'), bbox_inches='tight', pad_inches=0)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'axial')
plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_right_axial.png'), bbox_inches='tight', pad_inches=0)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 1, 'right', 'sagittal')
plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_right_sagittal.png'), bbox_inches='tight', pad_inches=0)

importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 20, 'left', 'coronal')
plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_left_coronal_shuffled_maps.png'), bbox_inches='tight', pad_inches=0)
importlib.reload(ff); ff.char_prf_angle_3D(AN, R2, THA, 500, 'left', 'coronal')
plt.savefig(os.path.join('../figures', 'fig_prf_body_char_angle_3D_left_coronal_shuffled_histogram.png'), bbox_inches='tight', pad_inches=0)

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

