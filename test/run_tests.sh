mkdir tests_output
echo Testing Legacy download
everystamp --survey legacy --ra 165.345 --dec 65.567 --size 0.1 --mode jpeg --legacy_bands grz --legacy_layer ls-dr9 --legacy_autoscale --download_dir tests_output

echo Testing PANSTARRS download
everystamp --survey panstarrs --ra 165.345 --dec 65.567 --size 0.1 --mode jpeg --ps_bands gri --download_dir tests_output