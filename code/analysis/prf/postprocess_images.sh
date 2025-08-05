# Crop two-sided surface plots to nice tidy versions
cd ~/Guestetal2025_NSDPulvinar/figures_bulk
mkdir prf_cropped
magick mogrify -monitor -crop 1050x1050+20+10 -path prf_cropped *fsaverage*.png