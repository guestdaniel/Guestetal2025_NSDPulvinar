# Rsync
#rsync -rv dguest@stone.cmrr.umn.edu:~dguest/Guestetal2021/../figures_bulk ~/Guestetal2021

# Change directory
cd ~/Guestetal2021/../figures_bulk

# Do cropping and var such stuff
for hemi in lh rh
do
	for vers in 1 2
	do
		for var in 1 2
		do
			for label in 1 2 3
			do
				for crange in 1 2 
				do
				# Do maps that are thresholded by tval
			    convert -crop 1650x1450+10+20 +repage contrast_${hemi}_${vers}_${var}_${label}_${crange}_threshold_by_tval.png contrast_${hemi}_${vers}_${var}_${label}_${crange}_cropped.png
			    convert -crop 1650x1450+10+20 +repage body_${hemi}_${vers}_${var}_${label}_${crange}_threshold_by_tval.png body_${hemi}_${vers}_${var}_${label}_${crange}_cropped.png

			    # Do maps that are thresholded manually
			    convert -crop 1650x1450+10+20 +repage contrast_${hemi}_${vers}_${var}_${label}_${crange}.png contrast_${hemi}_${vers}_${var}_${label}_${crange}_cropped_plain.png
			    convert -crop 1650x1450+10+20 +repage body_${hemi}_${vers}_${var}_${label}_${crange}.png body_${hemi}_${vers}_${var}_${label}_${crange}_cropped_plain.png
			    done
			done
		done
	done
done

# Also crop the plain surface image the same way
convert -crop 1650x1450+10+20 +repage empty_surface.png empty_surface_cropped.png