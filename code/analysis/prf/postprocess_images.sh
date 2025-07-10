# Crop two-sided surface plots to nice tidy versions
cd ~/Guestetal2025_NSDPulvinar/figs_bulk
mkdir cropped
magick mogrify -crop 1050x800+40+20 -path cropped *.png