# Change directory
cd figures

# Make side by side comparisons (contrast peak on top, body peak on bottom)
for hemi in lh rh
do
	for vers in 1 2
	do
		for label in 1 2 3
		do
		  convert -crop 1650x1450+10+20 +repage contrast_${hemi}_${vers}_1_${label}.png contrast_${hemi}_${vers}_1_${label}_cropped.png
		  convert -crop 1650x1450+10+20 +repage body_${hemi}_${vers}_1_${label}.png body_${hemi}_${vers}_1_${label}_cropped.png
		done
	done
done

# Make side by side comparisons (contrast peak on top, body peak on bottom)
for hemi in rh
do
	for vers in 1 2
	do
		for label in 1 2 3
		do
		  convert -crop 1650x1450+10+20 +repage contrast_${hemi}_${vers}_na_${label}_1_bootstrapped_thr_95_7.png contrast_${hemi}_${vers}_na_${label}_1_bootstrapped_thr_95_7_cropped.png
		  convert -crop 1650x1450+10+20 +repage body_${hemi}_${vers}_na_${label}_1_bootstrapped_thr_95_7.png body_${hemi}_${vers}_na_${label}_1_bootstrapped_thr_95_7_cropped.png
		done
	done
done
