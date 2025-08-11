# Program to sync CMRR data from old name paths to new name paths on desktop
#!/bin/bash

# Move bulk images for supplemental materials surface PRF data 
scp dguest@surly.cmrr.umn.edu:~/Guestetal2021/figs_output_bulk/fsaverage*.png C:\Users\Daniel\Desktop\Guestetal2025_NSDPulvinar\figures_bulk

# Move bulk images from supplemental stability assessment of sub_to_cor results
scp dguest@surly.cmrr.umn.edu:~/Guestetal2025_NSDPulvinar/figures_bulk/*supplemental*.png C:\Users\Daniel\Desktop\Guestetal2025_NSDPulvinar\figures_bulk

# Move NIFTIs containing correlation results to disk
scp dguest@surly.cmrr.umn.edu:/home/surly-raid3/dguest-data/Guestetal2021_data/subj01/mni/corr_cor_to_sub_hemi_*.nii.gz C:\Users\Daniel\Desktop\Guestetal2025_NSDPulvinar\?