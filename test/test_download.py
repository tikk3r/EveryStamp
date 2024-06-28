from test_basics import run_command

def test_legacy_cli():
    run_command('everystamp download --survey legacy --ra 165.345 --dec 65.567 --size 0.1 --mode jpeg --legacy_bands grz --legacy_layer ls-dr9 --legacy_autoscale')
    run_command('everystamp download --survey legacy --ra 202.4841667 --dec 47.2305556 --size 0.1 --mode fits --legacy_bands z --legacy_layer ls-dr9 --legacy_autoscale')

def test_LoLSS_cli():
    run_command('everystamp download --survey lolss --ra 202.4841667 --dec 47.2305556 --size 0.1 --mode fits --download_dir tests_output')

def test_LoTSS_cli():
    run_command('everystamp download --survey lotss --ra 202.4841667 --dec 47.2305556 --size 0.3 --mode fits --download_dir tests_output --lotss_release pdr')
    run_command('everystamp download --survey lotss --ra 202.4841667 --dec 47.2305556 --size 0.3 --mode fits --download_dir tests_output --lotss_release dr1')
    run_command('everystamp download --survey lotss --ra 202.4841667 --dec 47.2305556 --size 0.3 --mode fits --download_dir tests_output --lotss_release dr2')

def test_PanSTARRS_cli():
    run_command('everystamp download --survey pan-starrs --ra 165.345 --dec 65.567 --size 0.1 --mode jpeg --ps_bands gri --download_dir tests_output')

def test_TGSS_cli():
    run_command('everystamp download --survey tgss --ra 165.345 --dec 65.567 --size 0.1 --download_dir tests_output --mode fits')

def test_VLASS_cli():
    run_command('everystamp download --survey vlass --ra 165.345 --dec 65.567 --size 0.1 --download_dir tests_output --mode fits --vlass_type ql')
    #run_command('everystamp download --survey vlass --ra 180.345 --dec 60.567 --size 0.1 --download_dir tests_output --mode fits --vlass_type se')
