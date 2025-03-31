python3 copy_existing_data.py
matlab -nodisplay -nosplash -nodesktop -r "run('~/Guestetal2021/preproc/transform_thalamus_anat0pt8_to_MNI.m');exit;" | tail -n +11
matlab -nodisplay -nosplash -nodesktop -r "run('~/Guestetal2021/preproc/transform_prf_func1mm_to_MNI.m');exit;" | tail -n +11
matlab -nodisplay -nosplash -nodesktop -r "run('~/Guestetal2021/preproc/transform_postthalamus_0pt5mm_to_MNI_and_func1mm.m');exit;" | tail -n +11
matlab -nodisplay -nosplash -nodesktop -r "run('~/Guestetal2021/preproc/transform_nsdprf_func1mm_to_MNI.m');exit;" | tail -n +11
matlab -nodisplay -nosplash -nodesktop -r "run('~/Guestetal2019/preproc/transform_corticalroi_subjnat_to_fsaverage.m');exit;" | tail -n +11
