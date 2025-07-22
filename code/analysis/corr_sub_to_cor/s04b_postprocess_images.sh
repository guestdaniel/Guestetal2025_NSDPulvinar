# Crop two-sided surface plots to nice tidy versions
cd ~/Guestetal2025_NSDPulvinar/figures_bulk
mkdir cropped
magick mogrify -monitor -crop 1650x1450+20+10 -path cropped *supplemental*.png

# Montage relevant examples
cd cropped
# Contrast method 1
magick montage sub_to_cor_map_supplemental_seed=contrast_Nmax_rh_method=1_threshold=2_label=1_voxel=*.png -tile 4x -geometry +2+2 contrast_montage_method_1.png
# Contrast method 2
magick montage sub_to_cor_map_supplemental_seed=contrast_Nmax_rh_method=2_threshold=2_label=1_voxel=*.png -tile 4x -geometry +2+2 contrast_montage_method_2.png
# Body method 1
magick montage sub_to_cor_map_supplemental_seed=body_Nmax_rh_method=1_threshold=2_label=1_voxel=*.png -tile 4x -geometry +2+2 body_montage_method_1.png
# Body method 2
magick montage sub_to_cor_map_supplemental_seed=body_Nmax_rh_method=2_threshold=2_label=1_voxel=*.png -tile 4x -geometry +2+2 body_montage_method_2.png
