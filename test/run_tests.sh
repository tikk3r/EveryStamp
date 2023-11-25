rm -rf tests_output
mkdir tests_output

#
# Download tests
#
printf "Testing Legacy download"
everystamp download --survey legacy --ra 165.345 --dec 65.567 --size 0.1 --mode jpeg --legacy_bands grz --legacy_layer ls-dr9 --legacy_autoscale --download_dir tests_output
everystamp download --survey legacy --ra 202.4841667 --dec 47.2305556 --size 0.1 --mode fits --legacy_bands z --legacy_layer ls-dr9 --legacy_autoscale --download_dir tests_output

printf "\nTesting LoLSS download"
everystamp download --survey lolss --ra 202.4841667 --dec 47.2305556 --size 0.1 --mode fits --download_dir tests_output

printf "\nTesting LoTSS download"
everystamp download --survey lotss --ra 202.4841667 --dec 47.2305556 --size 0.3 --mode fits --download_dir tests_output --lotss_release pdr
everystamp download --survey lotss --ra 202.4841667 --dec 47.2305556 --size 0.3 --mode fits --download_dir tests_output --lotss_release dr1
everystamp download --survey lotss --ra 202.4841667 --dec 47.2305556 --size 0.3 --mode fits --download_dir tests_output --lotss_release dr2

printf "\nTesting Pan-STARRS download"
everystamp download --survey pan-starrs --ra 165.345 --dec 65.567 --size 0.1 --mode jpeg --ps_bands gri --download_dir tests_output

printf "\nTesting TGSS download"
everystamp download --survey tgss --ra 165.345 --dec 65.567 --size 0.1 --download_dir tests_output --mode fits

printf "\nTesting VLASS download"
everystamp download --survey vlass --ra 165.345 --dec 65.567 --size 0.1 --download_dir tests_output --mode fits
everystamp download --survey vlass --ra 180.345 --dec 60.567 --size 0.1 --download_dir tests_output --mode fits --vlass_type se

#
# Plotting tests
#
printf "\nTesting basic plotting"
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits

printf "\nTesting CLAHE"
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --CLAHE --CLAHE-gridsize 21 --CLAHE-cliplim 1

printf "\nTesting gamma"
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --gamma 2.2

#
# Contour plotting tests
#

printf "\nTesting HDR tonemapping"
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap ashikhmin
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap drago
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap duran
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap fattal
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap ferradans
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap ferwerda
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap vanhateren
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap kimkautz
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap lischinski
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap mantiuk06
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap mantiuk08
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap reinhard02
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap reinhard05
everystamp plot --image tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits --hdr-tonemap pattanaik
