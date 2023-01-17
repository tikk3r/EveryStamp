mkdir tests_output
echo Testing Legacy download
everystamp download --survey legacy --ra 165.345 --dec 65.567 --size 0.1 --mode jpeg --legacy_bands grz --legacy_layer ls-dr9 --legacy_autoscale --download_dir tests_output

echo \\nTesting LoLSS download
everystamp download --survey lolss --ra 202.4841667 --dec 47.2305556 --size 0.1 --mode fits --download_dir tests_output

echo \\nTesting LoTSS download
everystamp download --survey lotss --ra 202.4841667 --dec 47.2305556 --size 0.3 --mode fits --download_dir tests_output --lotss_release pdr
everystamp download --survey lotss --ra 202.4841667 --dec 47.2305556 --size 0.3 --mode fits --download_dir tests_output --lotss_release dr1
everystamp download --survey lotss --ra 202.4841667 --dec 47.2305556 --size 0.3 --mode fits --download_dir tests_output --lotss_release dr2

echo \\nTesting Pan-STARRS download
everystamp download --survey pan-starrs --ra 165.345 --dec 65.567 --size 0.1 --mode jpeg --ps_bands gri --download_dir tests_output

echo \\nTesting TGSS download
everystamp download --survey tgss --ra 165.345 --dec 65.567 --size 0.1 --download_dir tests_output --mode fits

echo \\nTesting VLASS download
everystamp download --survey vlass --ra 165.345 --dec 65.567 --size 0.1 --download_dir tests_output --mode fits
