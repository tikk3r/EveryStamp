![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Ftikk3r%2FEveryStamp%2Fmaster%2Fpyproject.toml)
![PyPI - Version](https://img.shields.io/pypi/v/EveryStamp)

# EveryStamp
EveryStamp aims to be an over-arching package to easily obtain postage stamps from a variety of surveys.

# Installation
EveryStamp can be installed through pip from PyPI via

    pip install EveryStamp

or directly from the git repository via

    pip install git+https://github.com/tikk3r/EveryStamp.git

# Usage
Four main types of functionality are currently supported: downloading cutouts, simple plotting, making cutouts from local FITS files and more advanced image composites:

    everystamp -h
    usage: everystamp [-h] {download,plot,cutout} ...

    EveryStamp 1.6.0 by Frits Sweijen

    subcommands:
      Description of sub commands.

      {download,plot,cutout}
        download            Download a cutout from a specified survey.
        plot                Plot a user-supplied FITS image.
        cutout              Cut a user-supplied FITS image to size.
        composite           Create composite images with more advanced blending modes.

See the help of each sub command for a complete overview of available options. For more detailed examples, see the [documentation](https://tikk3r.github.io/EveryStamp/).

# Acknowledgements
If you use EveryStamp in your work, please consider acknowledging this package through

    This work made use of EveryStamp\footnote{https://tikk3r.github.io/EveryStamp/}.

and acknowledging any surveys you used by their preferred acknowledgements.
