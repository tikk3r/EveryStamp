import subprocess

def run_command(cmd: str) -> int:
    """_summary_

    Parameters
    ----------
    cmd : str
        Command to execute.

    Raises
    ------
    RuntimeError
        If execution failed.

    Returns
    -------
    success : bool
        True if the command exited successfully (return code 0)
    """
    cmd_list = cmd.split(' ')
    output = subprocess.run(cmd_list)
    success = True if output.returncode == 0 else False
    if not success:
        raise RuntimeError('Command "' + cmd + '" failed.')
    return success

def test_help():
    run_command('everystamp -h')

def test_help_download():
    run_command('everystamp download -h')

def test_help_plot():
    run_command('everystamp plot -h')

def test_help_cutout():
    run_command('everystamp cutout -h')