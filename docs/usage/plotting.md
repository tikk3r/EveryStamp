---
title: Plotting FITS images
layout: default
nav_order: 2
parent: General usage
---

<script type="text/javascript" id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>

# Plotting FITS images
{: .no_toc}
This page gives a brief overview of how to make basic plots of FITS files.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## Basic sky plot
The simplest plotting command is

> everystamp plot --image \<fits image>

This will create a simple RA-DEC sky plot of the provided image using a linear stretch.


# Image manipulation
EveryStamp offers image manipulation that can be used to visualise an image to your liking. The options outlined below are applied in the following order, if used:

1. HDR tonemapping
2. CLAHE
3. Gamma correction

## Gamma correction
[Gamma correction](https://en.wikipedia.org/wiki/Gamma_correction) can be applied as a final step by adding `--gamma <gamma>` to the plotting arguments. As implemented here, this will gamma stretch an image by

$$
I_\mathrm{out} = I_\mathrm{in}^{1 / \gamma}
$$

A value of $\gamma = 2$ results in a square root stretch of the input image.

## Contrast limited adaptive histogram equalisation
[Contrast limited adaptive histogram equalisation](https://en.wikipedia.org/wiki/Adaptive_histogram_equalization) (CLAHE) tries to improve contrast in the image by applying a local histogram equalisation while simultaneously limiting the contrast enhancement in an attempt to reduce the amplification of noise. It is applied as a second-to-last step, before any gamma correction, and can be enabled by using `--CLAHE`. The `--CLAHE-gridsize` and `--CLAHE-cliplim` control the size of the local region and the amount of contrast enhancement, respectively.

## HDR tonemapping
When dealing with a high dynamic range (HDR) image, simple tonemapping operators may not be enough. If LuminanceHDR is installed (see [installation]({{ site.baseurl }}{% link docs/installation.md %})), several HDR tonemapping algorithms are available through the `--hdr-tonemap` argument.