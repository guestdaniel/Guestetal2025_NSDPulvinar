# Python script to move existing/pre-prepared data files (e.g., T1 anatomy) to Guestetal2021_data.
import os
import sys
from shutil import copyfile
os.chdir("code")

# Set up paths (Old for surly)
# main_path = '/home/surly-raid4/kendrick-data/nsd/nsddata/ppdata/'
# nsdprf_path = '/home/stone-ext4/generic/Dropbox/nsdfaceprf/freesurfer'
# output_path = '/home/surly-raid3/dguest-data/Guestetal2021_data'

# Set up paths (new for public)
main_path = 'data/nsd/nsddata/ppdata/'
nsdprf_path = 'data/nsdfaceprf/freesurfer'
output_path = 'data/prepared'

# Configure
subjs = ['subj01', 'subj02', 'subj03', 'subj04', 'subj05', 'subj06', 'subj07', 'subj08']

# [0] Make empty directories
for subj in subjs+['group']:
        os.makedirs(os.path.join(output_path, subj), exist_ok=True)
        for space in ['fsaverage', 'mni', 'func1mm']:
                os.makedirs(os.path.join(output_path, subj, space), exist_ok=True)

# [1] Copy anatomy in func1mm and MNI
for subj in subjs:
        for image in ['T1', 'T2', 'TOF', 'SWI', 'EPI']:
                # MNI space
                copyfile(os.path.join(main_path, subj, 'anat', image + '_to_MNI.nii.gz'),
                         os.path.join(output_path, subj, 'mni', image + '.nii.gz'))
        for image in ['T1', 'T2', 'TOF', 'SWI']:
                # func1mm space
                copyfile(os.path.join(main_path, subj, 'func1mm', image + '_to_func1mm.nii.gz'),
                         os.path.join(output_path, subj, 'func1mm', image + '.nii.gz'))

# [2] Copy thalamus ROIs
for subj in subjs:
        for roi in ['thalamus']:
                copyfile(os.path.join(main_path, subj, 'func1mm', 'roi', roi + '.nii.gz'),
                         os.path.join(output_path, subj, 'func1mm', roi + '.nii.gz'))

# [3] Copy func1mm pRF data
for subj in subjs:
        for image in ['prf_angle', 'prf_eccentricity', 'prf_R2', 'prf_size']:
                copyfile(os.path.join(main_path, subj, 'func1mm', image + '.nii.gz'),
                         os.path.join(output_path, subj, 'func1mm', image + '.nii.gz'))

# [4] Copy func1mm nsd-prf data
for subj in subjs:
        for image in ['_angle', '_eccentricity', '_R2', '_size']:
                for prf_type in ['backgroundauto', 'bodyauto', 'contrastNEW', 'faceauto', 'foregroundauto', 'salience', 'wordauto']:
                        copyfile(os.path.join(nsdprf_path, subj, 'label', prf_type + image + '.nii.gz'),
                                 os.path.join(output_path, subj, 'func1mm', prf_type + image + '.nii.gz'))

# [5] Copy Kastner ROIs
copyfile(os.path.join('/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/', 'fsaverage', 'label', 'lh.Kastner2015.mgz'),
         os.path.join(output_path, 'group', 'fsaverage', 'lh.Kastner2015.mgz'))
copyfile(os.path.join('/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/', 'fsaverage', 'label', 'rh.Kastner2015.mgz'),
         os.path.join(output_path, 'group', 'fsaverage', 'rh.Kastner2015.mgz'))

