# Program to sync CMRR data from old name paths to new name paths on desktop
#!/bin/bash

# Move bulk images for supplemental materials surface PRF data 
scp dguest@surly.cmrr.umn.edu:~/Guestetal2021/figs_output_bulk/fsaverage*.png C:\Users\Daniel\Desktop\Guestetal2025_NSDPulvinar\figures_bulk

# Move bulk images from supplemental stability assessment of sub_to_cor results
scp dguest@surly.cmrr.umn.edu:~/Guestetal2025_NSDPulvinar/figures_bulk/*supplemental*.png C:\Users\Daniel\Desktop\Guestetal2025_NSDPulvinar\figures_bulk

# Move cluster labels manually for a test subject, later we use tar method below
scp dguest@surly.cmrr.umn.edu:/home/surly-raid3/dguest-data/Guestetal2021_data/subj01/mni/beta_clusters*.nii.gz data\prepared\subj01\mni
scp dguest@surly.cmrr.umn.edu:/home/surly-raid3/dguest-data/Guestetal2021_data/subj02/mni/beta_clusters*.nii.gz data\prepared\subj02\mni

# tar all data from `data_dir` and then copy and untar va scp
#  on remote side:
cd /home/surly-raid3/dguest-data
tar -cvf data.tar --exclude='*.mat' Guestetal2021_data
#  on local side:
cd C:\Users\Daniel\Desktop\Guestetal2025_NSDPulvinar
scp dguest@surly.cmrr.umn.edu:/home/surly-raid3/dguest-data/data.tar data
cd data
tar -xf data.tar -C prepared
cd ..