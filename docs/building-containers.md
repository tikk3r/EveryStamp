---
title: Building containers
layout: default
nav_order: 2
---

# Building containers
{: .no_toc}
## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

Rebuilding containers is necessary if you want to change or add software, or optimise for a specific machine.

---

The recipes offer customisation options through environment variables near the top. Environment variables that a user may want to modify are summarised in the table below. Other variables should not requrie modification unless there are specific reasons for doing so.

| Variable | Values | Purpose |
|----------|--------|---------|
|`HAS_CUDA`|true/false| Install Nvidia CUDA and link IDG to it if set to `true`. |
|`HAS_MKL`| true/false| Install the Intel Math Kernel Library en link IDG to it if set to `true`.|
|`MARCH`|[GCC compatible architecture](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html)| Passed to GCC's `-march` option to generate code optimised for the specified CPU architecture.|
|`MTUNE`|[GCC compatible architecture](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html)| Passed to GCC's `-mtune` option to optimise code for the specified CPU architecture.|
|`NOAVX512`|true/false| Disables the use of AVX512 instructions if set to `true`. Recommended when building for older or generic hardware on a modern machine, or when building AMD-specific containers on an Intel machine.|
|`CPPSTD`|[C++ standard](https://gcc.gnu.org/projects/cxx-status.html)| Passed to GCC's `-std` option when building C++ code, specifying which C++ standard to use.|

{: .warning}
> In 2021 Singularity was renamed to Apptainer. There is little if any change in functionality, but if apptainer-related commands or environment variables mentioned below do not work for you, replace instances of `apptainer` with `singularity` and instances of `APPTAINER_` with `SINGULARITY_`. See [the official announcement](https://apptainer.org/news/community-announcement-20211130/) for more information.

### Building with sudo rights
If have root permissions on the machine you are building on, a container can be built with
```bash
sudo apptainer build <container name> <recipe file>
```

### Building without sudo rights
If you do not have root permissions on the machine you are building on, the Singularity/Apptainer installation has to have the fakeroot feature enabled. Once that is available, a container can be built by running
```bash
apptainer build --fakeroot <container name> <recipe file>
```

### Building with small temporary directories
If the machine you are running at has little space available on e.g. `/var` or `/tmp` you may run into trouble during the build. To avoid this, the build's temporary and cache directories can be set as follows:
```bash
sudo APPTAINER_TMPDIR=<custom tempdir> APPTAINER_CACHEDIR=<custom cachedir> apptainer build <container name> <recipe file>
```

## Building a generic container
A generic container is intended to be compatible with a wide variety of systems. This entails:

1. Setting `MARCH=x86-64` and `MTUNE=generic` to disable CPU architecture optmisation by setting.
2. Setting `NOAVX512=true` to disable AVX512 instructions.
3. Setting `HAS_MKL=false` to not use Intel MKL.

This should produce a container capable of running on a reasonably wide range of CPUs launched in the past decade or so.

## Building an optimised container
Enabling compiler optimisations can result is a substantial performance increase over a generic container. This entails:

1. Setting `MARCH` and `MTUNE` to values suitable for your CPU.
2. Setting `NOAVX512=false` if your CPU supports AVX512 instructions.
3. Setting `HAS_MKL=true` if you have an Intel CPU or if you have an AMD CPU using the AMD-specific recipe that uses the suite of AMD Optimised CPU Libraries (AOCL).

You can find the recommended march, mtune and AVX512 settings by running [`obtain_march_mtune.sh`](https://github.com/tikk3r/lofar-grid-hpccloud/blob/fedora-py3/obtain_march_mtune.sh). It is important to use a recent GCC such that it can recognise your CPU properly. The container will currently use GCC 11.