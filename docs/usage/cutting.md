---
title: Cutouts from local images
layout: default
nav_order: 3
parent: General usage
---

<script type="text/javascript" id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>

# Making cutouts from local FITS images
{: .no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## Single cutouts
The simplest cutout command is

> everystamp cutout --image \<fits image> --ra \<ra> --dec \<dec> --size \<size>

where the position and size are in degrees. The option `--cutout-mode` option can be set to either `astropy` to use its Cutout2D, or `fast` to try and simply slice a subarray out of the main image. The latter can be faster for large images.

## Batch cutouts from a catalogue
If you have a catalogue with RA and DEC columns, EveryStamp can make cutouts of all those sources automatically by specifcying the catalogue through `--from_catalogue`.

> everystamp cutout --image \<fits image> --ra \<ra> --dec \<dec> --size \<size> --from_catalogue \<mycatalogue>
