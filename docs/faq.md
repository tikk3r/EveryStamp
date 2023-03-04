---
title: FAQ
layout: default
nav_order: 99
---

# Frequently asked questions
{: .no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## I see warnings about a build mismatch when starting a container.
If you see the following message when entering a container (bar specific values)

> WARNING - software has been built with -march=znver2 but current machine reports -march=bdver4.
> If you encounter strange behaviour or Illegal instruction warnings, consider building a container with the appropriate architecture set.
>
> WARNING - software has been built with -mtune=znver2 but current machine -mtune=bdver4.
> If you encounter strange behaviour or Illegal instruction warnings, consider building a container with the appropriate architecture set.

it means the container was built for a different CPU architecture than the machine you intend to run on. If the container was built with `-march=x86-64 -mtune=generic`, this warning should be safe to ignore. If any other mismatch is reported, successful software runs will depend on how far apart the target architecture and your architecture are in terms of supported features. Running older containers on newer hardware will generally work and even a newer generation does not necessarily guarantee problems. Your mileage may vary and you have three options: use at your own risk, obtain a container built for your architecture or obtain a generic container.

## My software crashes with "Illegal instruction"
The software you are trying to use was built for a different architecture than your machine. If you see the warning of the previous question, then this is probably the result of attempting to use an incompatible container. If you think that this is a bug and that it should work, feel free to open a [ticket](https://github.com/tikk3r/lofar-grid-hpccloud/issues).

## Directories have disappeared or software can't find files that exist.
If directories or files that exist on the host system cannot be found inside the container, double check if all the required directories have been passed along to `--bind/-B` before entering the container.