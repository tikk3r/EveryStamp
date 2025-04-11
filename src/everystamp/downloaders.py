import glob
import logging
import os
import subprocess
import sys
from collections.abc import Iterable
from urllib.request import urlopen

import astropy.units as u
import numpy as np
import pyvo
import requests
import tqdm
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astroquery.hips2fits import hips2fits
from astroquery.skyview import SkyView

logging.basicConfig(
    format="[%(name)s] %(asctime)s - %(levelname)s: %(message)s", level=logging.INFO
)
logger = logging.getLogger("EveryStamp:Downloader")


def flatten(xs):
    """Flatten a nested list, as per the example on https://stackoverflow.com/a/40857703.

    Args:
        xs : list
            Nested list to flatten.

    Returns:
        x : list
            Flattened list.
    """
    for x in xs:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            yield from flatten(x)
        else:
            yield x


class FileDownloader(object):
    """From https://medium.com/better-programming/python-progress-bars-with-tqdm-by-example-ce98dbbc9697
    Copyright 2019 tiptapcode Authors. All Rights Reserved.
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
         http://www.apache.org/licenses/LICENSE-2.0
    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
    """

    def get_url_filename(self, url):
        """
        Discover file name from HTTP URL, If none is discovered derive name from http redirect HTTP content header Location
        :param url: Url link to file to download
        :type url: str
        :return: Base filename
        :rtype: str
        """
        try:
            filename = os.path.basename(url)
            basename, ext = os.path.splitext(filename)
            if ext:
                return filename
            header = requests.head(url, allow_redirects=False).headers
            return (
                os.path.basename(header.get("Location"))
                if "Location" in header
                else filename
            )
        except requests.exceptions.HTTPError as errh:
            print("Http Error:", errh)
            raise errh
        except requests.exceptions.ConnectionError as errc:
            print("Error Connecting:", errc)
            raise errc
        except requests.exceptions.Timeout as errt:
            print("Timeout Error:", errt)
            raise errt
        except requests.exceptions.RequestException as err:
            print("OOps: Something Else", err)
            raise err

    def download_file(self, url, filename=None, target_dir=None):
        """
        Stream downloads files via HTTP
        :param url: Url link to file to download
        :type url: str
        :param filename: filename overrides filename defined in Url param
        :type filename: str
        :param target_dir: target destination directory to download file to
        :type target_dir: str
        :return: Absolute path to target destination where file has been downloaded to
        :rtype: str
        """
        if target_dir and not os.path.isdir(target_dir):
            raise ValueError("Invalid target_dir={} specified".format(target_dir))
        local_filename = self.get_url_filename(url) if not filename else filename

        req = requests.get(url, stream=True, verify=True)
        req.raise_for_status()
        try:
            file_size = int(req.headers["Content-Length"])
        except KeyError:
            if req.headers["Transfer-Encoding"] == "chunked":
                file_size = 0
        chunk_size = 1024  # 1 MB
        num_bars = int(file_size / chunk_size)

        base_path = os.path.abspath(os.path.dirname(__file__))
        target_dest_dir = (
            os.path.join(base_path, local_filename)
            if not target_dir
            else os.path.join(target_dir, local_filename)
        )
        with open(target_dest_dir, "wb") as fp:
            for chunk in tqdm.tqdm(
                req.iter_content(chunk_size=chunk_size),
                total=num_bars,
                unit="KB",
                leave=True,
                file=sys.stdout,
            ):
                fp.write(chunk)

        return target_dest_dir


class LoTSSDownloader(FileDownloader):
    """Downloader sub-class for the LOFAR Two-metre Sky Survey."""

    supported_keywords = ["ra", "dec", "size_arcmin"]
    logger = logging.getLogger("EveryStamp:LoTSSDownloader")

    def __init__(self):
        self.url = "https://lofar-surveys.org/{release:s}-cutout.fits?pos={coord_str:s}&size={size_arcmin:f}"

    def format_url(self, ra: str, dec: str, release: str, size: float = 1.0) -> str:
        """

        Args:
            ra: right ascension in HH:MM:SS format.
            dec: declination in HH:MM:SS format.
            size: size of the cutout area in degrees.

        Returns:
            url: the formatted URL.
        """
        size_arcmin = size * 60
        coord_str = (
            SkyCoord(ra, dec, unit="deg", frame="icrs")
            .to_string("hmsdms")
            .replace("h", ":")
            .replace("m", ":")
            .replace("d", ":")
            .replace("s", "")
        )
        url = self.url.format(
            coord_str=coord_str, size_arcmin=size_arcmin, release=release
        )
        return url

    def download(self, **kwargs):
        if kwargs["mode"] != "fits":
            raise ValueError("LoTSSDownloader only supports FITS downloads.")
        furl = self.format_url(
            ra=kwargs["ra"],
            dec=kwargs["dec"],
            size=kwargs["size"],
            release=kwargs["release"],
        )
        logger.info(furl)
        fname = "LoTSS-{release:s}_{ra:f}_{dec:f}_{size:.3f}.{mode:s}".format(
            ra=kwargs["ra"],
            dec=kwargs["dec"],
            mode="fits",
            size=kwargs["size"],
            release=kwargs["release"].upper(),
        )
        if not kwargs["ddir"]:
            self.logger.info(
                "Download directory not specified, downloading to %s instead",
                os.getcwd(),
            )
            ddir = os.getcwd()
        else:
            ddir = kwargs["ddir"]
        self.download_file(furl, filename=fname, target_dir=ddir)


class LegacyDownloader(FileDownloader):
    """Downloader sub-class for the DESI Legacy Imaging Surveys."""

    supported_keywords = ["ra", "dec", "mode", "layer", "pixscale", "bands", "size_pix"]
    logger = logging.getLogger("EveryStamp:LegacyDownloader")

    def __init__(self):
        self.url = "https://www.legacysurvey.org/viewer/{mode:s}-cutout/?ra={ra:f}&dec={dec:f}&layer={layer:s}&pixscale={pixscale:.3f}&bands={bands:s}&size={size_pix:d}"

    def format_url(
        self,
        ra=None,
        dec=None,
        size=None,
        bands="grz",
        mode="jpeg",
        layer="ls-dr9",
        pixscale=0.262,
        autoscale=False,
        **kwargs,
    ):
        """Returns a properly formatted URL that can be used to obtain a cutout from Legacy."""
        size_pix = int(size * 3600 / pixscale)
        dlpixscale = pixscale
        dlsize_pix = size_pix
        if (size_pix > 3000) and autoscale:
            # Jump to the next available pixel size by scaling from the (approximate) native pixel scale.
            new_pixscale = pixscale
            new_size_pix = int(size * 3600 / new_pixscale)
            while new_size_pix > 3000:
                new_pixscale += 0.262
                new_size_pix = int(size * 3600 / new_pixscale)
            self.logger.warn(
                "Image size of {:.2f} deg with pixel scale {:.3f} exceeds server limit of 3000 pixels! Automatically adjusting pixel scale to {:.3f} giving {:d} pixels.".format(
                    size, pixscale, new_pixscale, new_size_pix
                ),
                stacklevel=2,
            )
            dlpixscale = new_pixscale
            dlsize_pix = new_size_pix
        elif size_pix > 3000:
            self.logger.warn(
                "Image size of {:.2f} deg with pixel scale {:.3f} exceeds server limit of 3000 pixels! Image will be truncated! Use --legacy_autoscale or pass autoscale=True to automatically switch pixel scales.".format(
                    size, pixscale
                ),
                stacklevel=2,
            )
        return self.url.format(
            ra=ra,
            dec=dec,
            size_pix=dlsize_pix,
            bands=bands,
            mode=mode,
            layer=layer,
            pixscale=dlpixscale,
        )

    def download(self, **kwargs):
        if kwargs["mode"] == "both":
            furl = self.format_url(**kwargs)
            self.logger.info("Downloading cutout from %s", furl)
            if not kwargs["ddir"]:
                self.logger.info(
                    "Download directory not specified, downloading to %s instead",
                    os.getcwd(),
                )
                ddir = os.getcwd()
            else:
                ddir = kwargs["ddir"]
            fname = "legacystamps_{ra:f}_{dec:f}_{layer:s}.{mode:s}".format(
                ra=kwargs["ra"],
                dec=kwargs["dec"],
                layer=kwargs["layer"],
                mode="jpeg",
            )
            self.download_file(furl, filename=fname, target_dir=ddir)

            # Download FITS
            self.logger.info("Downloading cutout from %s", furl)
            furl = furl.replace("jpeg", "fits")
            fname = fname.replace("jpeg", "fits")
            self.download_file(furl, filename=fname, target_dir=ddir)
        else:
            furl = self.format_url(**kwargs)
            self.logger.info("Downloading cutout from %s", furl)
            if not kwargs["ddir"]:
                self.logger.info(
                    "Download directory not specified, downloading to %s instead",
                    os.getcwd(),
                )
                ddir = os.getcwd()
            else:
                ddir = kwargs["ddir"]
            fname = "legacystamps_{ra:f}_{dec:f}_{layer:s}.{mode:s}".format(
                ra=kwargs["ra"],
                dec=kwargs["dec"],
                layer=kwargs["layer"],
                mode=kwargs["mode"],
            )
        try:
            self.download_file(furl, filename=fname, target_dir=ddir)
        except requests.exceptions.HTTPError:
            self.logger.warning(f"Failed to download {fname}")


class PanSTARRSDownloader:
    """Downloader sub-class for the PanSTARRS survey."""

    from panstamps.downloader import downloader as psdownloader

    def __init__(self):
        self.logger = logging.getLogger("EveryStamp:Pan-STARRSDownloader")

    def download(self, ra, dec, size, mode="jpeg", ddir="", bands="gri"):
        if mode == "jpeg":
            get_jpeg = True
            get_fits = False
        elif mode == "fits":
            get_jpeg = False
            get_fits = True
        elif mode == "both":
            get_jpeg = True
            get_fits = True
        arcsecsize = size * 3600
        self.logger.info("Downloading cutout from PANSTARRS")
        d = self.psdownloader(
            ra=ra,
            dec=dec,
            downloadDirectory=ddir or os.getcwd(),
            fits=get_fits,
            jpeg=get_jpeg,
            color=True,
            singleFilters=True,
            filterSet=bands,
            imageType="stack",
            arcsecSize=arcsecsize,
            log=self.logger,
        )
        fitspath, jpegpath, colorpaths = d.get()
        # Rename the output slightly so the user can find it easier.
        for colorpath in colorpaths:
            os.rename(
                colorpath,
                os.path.join(
                    os.path.dirname(colorpath),
                    "panstamps_" + os.path.basename(colorpath),
                ),
            )


class VLASSDownloader(FileDownloader):
    """Downloader sub-class for the VLASS survey.

    Based on the original code by Anna Ho (https://github.com/annayqho/Query_VLASS) and edits by R. Timmerman.
    """

    def __init__(self, datatype: str = "ql"):
        self.summary_url = "https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php"
        self.logger = logging.getLogger("EveryStamp:VLASSDownloader")
        self.datatype = datatype.lower()
        if self.datatype != "ql" and self.datatype != "se":
            raise ValueError(
                "Unknown VLASS image type specificed. Choose either ql or se for Quick Look or Single Epoch."
            )
        if self.datatype == "ql":
            self.pixel_scale = 2.777777777778e-4 * 3600  # arcsec / pixel
        elif self.datatype == "se":
            self.pixel_scale = 1.666666666667e-04 * 3600

    def get_tiles(self, summary_file="VLASS_dyn_summary.php"):
        """
        Read tiles from tile catalog. If file missing, try wget https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php

        Parameters
        ----------
        summary_file : str
            Location where the VLASS summary file is located.

        Returns
        -------
        tab : astropy Table
            Table containing tile information.
        """
        if not os.path.isfile(summary_file):
            self.logger.warn(f"Could not find VLASS summary file {summary_file}!")
            self.logger.info("Attempting to download VLASS summary file")
            subprocess.run(
                [
                    "wget",
                    "https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php",
                    "-O",
                    summary_file,
                ]
            )

        # Put it in a more managable format by replacing consecutive white space with commas.
        # Assumes no more than 1 space in valid entries.
        subprocess.run(["sed", "-i", "-e", r"s/ \{2,\}/,/g", summary_file], check=True)

        tab = Table.read(summary_file, data_start=3, format="ascii.csv")
        tab.rename_column("Observing", "Epoch")
        tab.rename_column("Observation", "Date")
        return tab

    def search_tiles(self, tiles, c):
        """Search the tile catalog for tiles containing the input coordinate

        Parameters
        ----------
        tiles : astropy Table
            Tile catalogue as obtained from get_tiles()
        c : SkyCoord
            Location to check for coverage in one of the observed tiles.

        Returns
        -------
        tile name : str
            Best matching tile.
        observing epoch : str
            Epoch the best matching tile was observed in.
        observing date : str
            Date the best matching tile was observed at.
        """
        ra_h = c.ra.hour
        dec_d = c.dec.deg

        has_dec = np.logical_and(dec_d >= tiles["Dec min"], dec_d < tiles["Dec max"])
        has_ra = np.logical_and(ra_h >= tiles["RA min"], ra_h < tiles["RA max"])
        in_tile = np.logical_and(has_ra, has_dec)
        name = tiles["Tile"][in_tile]
        epoch = tiles["Epoch"][in_tile]
        date = tiles["Date"][in_tile]
        if len(name) == 0:
            raise IndexError("Zero VLASS tiles available for the given coordinate")
        c_grid = SkyCoord(
            7.5 * (tiles["RA min"][in_tile] + tiles["RA max"][in_tile]),
            0.5 * (tiles["Dec min"][in_tile] + tiles["Dec max"][in_tile]),
            unit="deg",
            frame="icrs",
        )
        dist = c_grid.separation(c)
        best_idx = np.argmin(dist)
        return name[best_idx], epoch[best_idx], date[best_idx]

    def get_subtiles(self, tilename, epoch, consider_QA_rejected):
        """For a given tile name, get the subtile filenames in the VLASS directory

        Parse those filenames and return a list of subtile RA and Dec.
        RA and Dec returned as a SkyCoord object

        Parameters
        ----------
        tilename : str
            Name of the tile to extract a subtile from.
        epoch : str
            The epoch the tile was observed in.
        consider_QA_rejected : bool
            Also consider tiles that did not pass the quality assurance checks.

        Returns
        -------
        fname : str
            Name of the subtile.
        c : astropy SkyCoord
            Coordinate of the subtile.
        """

        # Obtain the HTML for the given tile
        if "1." not in epoch:
            if self.datatype == "se":
                URLPATH = f"https://archive-new.nrao.edu/vlass/se_continuum_imaging/{epoch}/{tilename}"
            elif self.datatype == "ql":
                URLPATH = (
                    f"https://archive-new.nrao.edu/vlass/quicklook/{epoch}/{tilename}"
                )
            self.logger.info(f"Downloading from {URLPATH}")
            urlpath = urlopen(f"{URLPATH}")
        else:
            self.logger.info(
                f"Downloading from https://archive-new.nrao.edu/vlass/quicklook/{epoch}v2/{tilename}"
            )
            urlpath = urlopen(
                f"https://archive-new.nrao.edu/vlass/quicklook/{epoch}v2/{tilename}"
            )
        string = urlpath.read().decode("utf-8").split("\n")

        if self.datatype == "ql" and consider_QA_rejected:
            # Obtain the HTML for the QA Rejected
            urlpath_rejected = urlopen(
                f"https://archive-new.nrao.edu/vlass/quicklook/{epoch}v2/QA_REJECTED"
            )
            string += urlpath_rejected.read().decode("utf-8").split("\n")

        # Select only the subtile parts
        vals = np.array(
            [
                val.strip()
                for val in string
                if ("href" in val.strip()) and (tilename in val.strip())
            ]
        )

        # Select the coordinate part. You want the 'VLASS1.1.ql.T25t12.J150000+603000.10.2048.v1/' bit
        fname = np.array([val.split('"')[7] for val in vals])

        # Split out the actual coordinate string
        pos_raw = np.array(
            [val.split(".")[4] for val in fname if val.startswith("VLASS")]
        )

        if "-" in pos_raw[0]:
            # dec < 0
            ra_raw = np.array([val.split("-")[0] for val in pos_raw])
            dec_raw = np.array([val.split("-")[1] for val in pos_raw])
        else:
            # dec > 0
            ra_raw = np.array([val.split("+")[0] for val in pos_raw])
            dec_raw = np.array([val.split("+")[1] for val in pos_raw])
        ra = []
        dec = []
        for ii, val in enumerate(ra_raw):
            if val[1:3] == "24":
                rah = "00"
            else:
                rah = val[1:3]
            ra.append(f"{rah}h{val[3:5]}m{val[5:]}s")
            dec.append(f"{dec_raw[ii][:2]}d{dec_raw[ii][2:4]}m{dec_raw[ii][4:]}s")
        ra = np.array(ra)
        dec = np.array(dec)
        c = SkyCoord(
            ra, dec, frame="icrs"
        )  # .directional_offset_by(45*u.deg, 0.75*u.deg)
        return fname, c

    def get_cutout(self, imname, c, crop_scale):
        """Get a smaller cutout from the subtile.

        Parameters
        ----------
        imname : str
            Name of the image to make a cutout from.
        c : astropy SkyCoord
            Coordinate around which to make a cutout.
        crop_scale : int
            Size of the cutout in pixels.

        Returns
        -------
        output_fits : str
            Name of the output FITS file.
        """
        # Define output name
        output_fits = os.path.join(
            "/".join(imname.split("/")[:-1]),
            "VLASS_{:.6f}_{:.6f}_".format(c.ra.value, c.dec.value)
            + imname.split("/")[-1].rstrip(".fits")
            + "_poststamp.fits",
        )

        # Get header info
        hdu_list = fits.open(imname)
        header = hdu_list[0].header
        data = hdu_list[0].data[0, 0, :, :]

        # Obtain header and drop useless axes
        wcs = WCS(header)
        wcs = wcs.dropaxis(2).dropaxis(2)

        pixel_coords = skycoord_to_pixel(
            SkyCoord(c.ra.deg, c.dec.deg, unit="deg", frame="icrs"), wcs
        )

        if (
            pixel_coords[0] < 0
            or pixel_coords[1] < 0
            or pixel_coords[0] > data.shape[0]
            or pixel_coords[1] > data.shape[1]
        ):
            subprocess.call(f"rm -f {imname}", shell=True)
            raise Exception(
                "Requested coordinate not within the available subtiles. Consider running with consider_QA_rejected=True to also search additional subtiles which failed initial QA checks"
            )

        # Produce a cutout
        cutout = Cutout2D(data, c, (crop_scale, crop_scale), wcs=wcs)

        # Update the HDU
        hdu_list[0].data = cutout.data
        new_header = cutout.wcs.to_header()
        hdu_list[0].header.update(new_header)
        hdu_list[0].header.set("NAXIS", 4)
        hdu_list[0].header.insert("NAXIS2", ("NAXIS3", 1), after=True)
        hdu_list[0].header.insert("NAXIS3", ("NAXIS4", 1), after=True)
        hdu_list[0].header.remove("WCSAXES", ignore_missing=True)
        hdu_list[0].header.remove("MJDREF", ignore_missing=True)
        hdu_list[0].header.remove("MJD-OBS", ignore_missing=True)

        # Write the new fits
        hdu_list.writeto(output_fits, overwrite=True)

        # Cleanup
        subprocess.call(f"rm -f {imname}", shell=True)

        return output_fits

    def search_vlass(
        self,
        c,
        crop=False,
        crop_scale=256,
        consider_QA_rejected=False,
        ddir=os.getcwd(),
    ):
        """
        Searches the VLASS catalog for a source

        Parameters
        ----------
        c : astropy SkyCoord
            Coordinate to search for in tiles.
        crop : bool
            Make a cropped cutout of the area of interest.
        crop_scale : int
            Crop the cutout to this amount of pixels centred around c.
        consider_QA_rejected : bool
            Also consider tiles that failed the Quality Assurance checks.
        ddir : str
            Location to download the cutout to.

        Returns
        -------
        imname : str
            Name of the output image.
        """
        # Find the VLASS tile
        tiles = self.get_tiles()
        tilename, epoch, obsdate = self.search_tiles(tiles, c)
        if self.datatype == "se":
            # No other epoch available for now, so hardcode until the code is updated.
            epoch = "VLASS2.1"

        subtiles, c_tiles = self.get_subtiles(tilename, epoch, consider_QA_rejected)
        dist = c.separation(c_tiles)
        subtile = subtiles[np.argmin(dist)]

        imname = f"{subtile[:-1]}.I.iter1.image.pbcor.tt0.subim.fits"
        if len(glob.glob(imname)) == 0:
            if self.datatype == "ql":
                if "1." not in epoch:
                    url_get = f"https://archive-new.nrao.edu/vlass/quicklook/{epoch}/{tilename}/{subtile}"
                else:
                    url_get = f"https://archive-new.nrao.edu/vlass/quicklook/{epoch}v2/{tilename}/{subtile}"
            elif self.datatype == "se":
                # Dirty workaround until proper implementation.
                epoch = "VLASS2.1"
                imname = imname.replace("iter1", "iter3")
                url_get = f"https://archive-new.nrao.edu/vlass/se_continuum_imaging/{epoch}/{tilename}/{subtile}"
            fname = f"{url_get}{imname}"
            self.logger.info("Downloading to " + ddir)
            self.download_file(fname, target_dir=ddir)
            if self.datatype == "ql" and consider_QA_rejected:
                if "1." not in epoch:
                    url_get = f"https://archive-new.nrao.edu/vlass/quicklook/{epoch}/QA_REJECTED/{subtile}"
                else:
                    url_get = f"https://archive-new.nrao.edu/vlass/quicklook/{epoch}v2/QA_REJECTED/{subtile}"
                fname = f"{url_get}{imname}"
                self.download_file(fname, target_dir=ddir)
        if crop:
            out = self.get_cutout(os.path.join(ddir, imname), c, crop_scale=crop_scale)
            return out
        else:
            return imname

    def download(
        self,
        ra=0.0,
        dec=0.0,
        size=0.1,
        crop=True,
        consider_QA_rejected=False,
        ddir=os.getcwd(),
    ):
        """Download a cutout from the VLASS survey.

        Parameters
        ----------
        ra : float
            Right ascension of the coordinate of interest.
        dec : float
            Declination of the coordinate of interest.
        size : float
            Size of the area of interest in degrees.
        crop : bool
            Crop the image to the area of interest.
        consider_QA_rejected : bool
            Also consider tiles that failed the Quality Assurance checks.
        ddir : str
            Location to download the cutout to.
        """
        c = SkyCoord(ra, dec, unit="deg")

        crop_scale = size * 3600 / self.pixel_scale
        self.logger.info("Downloading cutout from VLASS")
        self.search_vlass(
            c,
            crop=crop,
            crop_scale=crop_scale,
            consider_QA_rejected=consider_QA_rejected,
            ddir=ddir,
        )


class VODownloader:
    """Downloader sub-class for surveys offeret through a VO."""

    def __init__(self, url, name=""):
        if not url:
            raise ValueError("VO url cannot be empty.")
        else:
            self.url = url
        if not name:
            self.name = self.url
            self.logger = logging.getLogger("EveryStamp:VODownloader")
        else:
            self.name = name
            self.logger = logging.getLogger(
                "EveryStamp:VODownloader[{:s}]".format(self.name)
            )

    def download(self, ra=0.0, dec=0.0, size=0.1, ddir=os.getcwd(), suffix=""):
        """Download a cutout from a VO survey.

        Parameters
        ----------
        ra : float
            Right ascension of the coordinate of interest in degrees.
        dec : float
            Declination of the coordinate of interest in degrees    .
        size : float
            Size of the area of interest in degrees.
        ddir : str
            Location to download the cutout to.
        """
        c = SkyCoord(ra, dec, unit="deg")

        vo = pyvo.sia.SIAService(self.url)
        query = vo.search(c, size=size)
        if not query:
            raise ValueError("Requested coordinates not covered by the specified VO!")
        im = query.getrecord(0)
        if im.format == "image/fits":
            self.logger.info("Downloading cutout from {:s}".format(self.name))
            im.cachedataset(
                filename=os.path.join(
                    ddir,
                    "{:s}_{:.4f}_{:.4f}_{:.3f}.fits".format(self.name, ra, dec, size),
                )
            )


class HiPSDownloader:
    """Sub-class to download a file from a HiPS image using hips2fits."""

    def __init__(self, hips, name=""):
        if not hips:
            raise ValueError("HiPS name cannot be empty.")
        else:
            self.hips = hips
        if not name:
            self.name = self.hips.replace("/", "_")
        else:
            self.name = name
        self.logger = logging.getLogger(
            "EveryStamp:HiPSDownloader[{:s}]".format(self.name)
        )
        self.logger.warning(
            "downloading from a HiPS survey. Scientific accuracy of resulting FITS images is not guaranteed!"
        )

    def download(
        self, ra=0.0, dec=0.0, size=0.1, ddir=os.getcwd(), pixsize=1.0, mode="jpg"
    ):
        """Download a cutout from a HiPS survey.

        Parameters
        ----------
        ra : float
            Right ascension of the coordinate of interest in degrees.
        dec : float
            Declination of the coordinate of interest in degrees    .
        size : float
            Size of the area of interest in degrees.
        ddir : str
            Location to download the cutout to.
        pixsize : float
            Pixel scale of the survey in arcsec. Default is 1.0 arcsec per pixel.
        mode : str
            What image format to download. Can be jpg or fits, default is jpg.
        """
        imsize = int(size / (pixsize / 3600))
        img = hips2fits.query(
            hips=self.hips,
            format=mode,
            width=imsize,
            height=imsize,
            projection="SIN",
            fov=size * u.deg,
            ra=ra * u.deg,
            dec=dec * u.deg,
        )
        if not os.path.exists(ddir):
            self.logger.info("Download directory does not exist, creating it")
            os.mkdir(ddir)
        if mode == "jpg":
            from PIL import Image

            imdata = Image.fromarray(img)
            imdata.save(
                os.path.join(
                    ddir,
                    "{:s}_{:.4f}_{:.4f}_{:.5f}.jpeg".format(self.name, ra, dec, size),
                )
            )
        elif mode == "fits":
            img.writeto(
                os.path.join(
                    ddir,
                    "{:s}_{:.4f}_{:.4f}_{:.5f}.fits".format(self.name, ra, dec, size),
                )
            )


class SkyViewDownloader:
    """Downloader sub-class for surveys offeret through a VO."""

    def __init__(self, survey):
        if not survey:
            raise ValueError("SkyView survey cannot be empty.")
        else:
            self.survey = survey
            self.logger = logging.getLogger(
                "EveryStamp:SkyViewDownloader[{:s}]".format(self.survey)
            )

    def download(
        self, ra=0.0, dec=0.0, size=0.1, pixsize=1, ddir=os.getcwd(), suffix=""
    ):
        """Download a cutout from a SkyView survey.

        Parameters
        ----------
        ra : float
            Right ascension of the coordinate of interest in degrees.
        dec : float
            Declination of the coordinate of interest in degrees    .
        size : float
            Size of the area of interest in degrees.
        pixsize
            Pixel size to use for the cutout.
        ddir : str
            Location to download the cutout to.
        """
        c = SkyCoord(ra, dec, unit="deg")

        sv = SkyView()
        pixsize_deg = pixsize / 3600.0
        pixels = size / pixsize_deg
        self.logger.info(
            f"Requesting a {size} x {size} deg cutout of {pixels} x {pixels} pixels."
        )
        hdul = sv.get_images(c, self.survey, radius=size * u.deg, pixels=int(pixels))
        if not hdul:
            raise ValueError(
                "SkyView did not return a result. If you requested a large cutout, try increasing the pixel size through --skyview_pixsize or reducing the cutout area."
            )
        hdul[0].writeto(
            os.path.join(
                ddir,
                "{:s}_{:.4f}_{:.4f}_{:.3f}.fits".format(
                    self.survey.replace(" ", "_"), ra, dec, size
                ),
            )
        )
