---
title: Blending multiple images
layout: default
nav_order: 3
parent: General usage
---

<script type="text/javascript" id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

See `everystamp composite --help` for the full list of arguments.

The `everystamp composite` module allows the composition of multiple images. The functionality offered at the moment is to compose images using (blending modes)[https://en.wikipedia.org/wiki/Blend_modes]. Several foreground images can be composited on top of a background image. Each image can be manipulated via:

* `--blend-modes`: sets the blending mode this image will use to merge with the previous layers.
* `--blend-opacities`: the opacity for the specific layer
* `--blend-cmaps`: the colour map to use for a layer

Some experimental presets are defined under the `--preset` argument. Currently they aim for either optical and LOFAR, or optical, x-ray and LOFAR.
