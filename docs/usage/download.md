---
title: Downloading cutouts
layout: default
nav_order: 1
parent: General usage
---

Downloading cutouts
{: .no_toc}
This page gives a brief overview of how to download cutouts with EveryStamp and which surveys are supported.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

The download subcommand allows you to download FITS or JPEG cutouts from various surveys. Available surveys include those offered by SkyView and a number of custom additions not (yet) present in SkyView.

## Basic cutout download
A basic cutout download requires five arguments and can be obtained via

> everystamp download --survey \<survey> --ra \<ra> --dec \<dec> --size \<size> --mode \<mode>

This will download a square cutout centred on the given right ascension and declination of the given size. The mode can be `fits`, `jpeg` or `both` to download the respecitve formats. Note that not all surveys support JPEG cutouts. Various surveys may also have further arguments specific to them. See `everystamp download -h` for a complete summary of available parameters.

## Optical surveys
### DESI Legacy Imaging Surveys
Cutouts of the [DESI Legacy Imaging Surveys](https://www.legacysurvey.org/) are provided through a wrapper aroud their own cutout service. Examples of how to query that service can be found at https://www.legacysurvey.org/viewer/urls.

### Pan-STARRS
Cutouts of the Pan-STARRS survey are provided through the [panstamps](https://panstamps.readthedocs.io/en/master/) package.

## Radio surveys
The table below summarises available radio surveys.

| Survey | Argument name | Frequency/bandwidth | Angular resolution | Description |
|--------|---------------|---------------------|--------------------|-------------|
| TGSS | tgss | 150 MHz | 25'' | TFIR GMRT Sky Survey |
| LoLSS | lolls | 54 MHz | 45'' (PDR), 15'' (DR1) | LoFAR LBA Sky Survey |
| LoTSS | lotss | 120-168 MHz | 25'' (PDR), 6'' (DR1, DR2) | LoFAR Two-metre Sky Survey |
| VLASS | vlass | 2-4 GHz | 2.5'' | VLA Sky Survey |

### LoFAR LBA Sky Survey
Two versions of the LoFAR LBA Sky Survey are available:

1. Preliminary data release (PDR): region around the HETDEX Spring Field with an angular resolution of 45'' and a median noise level of 5 mJy/beam.
2. Data Release 1 (DR1): great improvements over PDR with an angular resolution of 15'' and a median noise level of 1.5 mJy/beam.

### LoFAR Two-metre Sky Survey
Three versions of the LoFAR Two-metre Sky Survey (LoTSS) are available:

1. Preliminary data release (PDR): a preliminary data release covering the HETDEX Spring Field at an angular resolution of 25'' and a typical noise level of <0.5 mJy/beam.
2. Data Release 1 (DR1): the first data release, covering the HETDEX Spring Field with a greatly improved angular resolution of 6'' and a median noise level of ~71 uJy/beam.
3. Data Release 2 (DR2): the second data release, covering 27% of the northern sky at an angular resolution of 6'' and a median noise level of ~81 uJy/beam.

Cutous from the PDR or DR1 releases are made through the [PDR image cutout](https://vo.astron.nl/lofartier1/q_img/cutout/form) and [DR1 image cutout](https://vo.astron.nl/hetdex/lotss-dr1-img/cutout/form) VO services hosted by ASTRON. For DR2 such a service is not yet available. Instead, CDS' [hips2fits service](https://alasky.u-strasbg.fr/hips-image-services/hips2fits) is used to generate FITS cutouts from the [6'' DR2 HiPS](https://lofar-surveys.org/public_hips/LoTSS_DR2_high_hips/).

For detailed information about LoTSS, see https://lofar-surveys.org/.

### TIFR GMRT Sky Survey
[TIFR GMRT Sky Survey (TGSS) Alternative Data Release 1 (ADR1)](https://tgssadr.strw.leidenuniv.nl/doku.php) cutouts are provided through a wrapper around the [VO image cutout service](https://vo.astron.nl/tgssadr/q_fits/cutout/form) hosted at ASTRON.

### VLA Sky Survey
Cutouts from the VLA Sky Survey (VLASS) are provided by code based on scripts written by [annayqho](https://github.com/annayqho/Query_VLASS) and [RedTimm](https://github.com/RedTimm/VLASS-poststamp). It queries the VLASS observation summary hosted at https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php to determine the epoch and tile to download a cutout from.
