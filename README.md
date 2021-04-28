# pigx_sarscov2_ww

PiGx SARS-CoV2 wastewater sequencing pipeline

## Hacking

We recommend using Guix to hack on this pipeline.  To enter a
reproducible environment with a known-good version of Guix use this
(slow) command:

```sh
USE_GUIX_INFERIOR=t guix environment -m manifest.scm
```

To use your current Guix channels instead of the fixed set of
channels, just omit the `USE_GUIX_INFERIOR` shell variable:

```sh
guix environment -m manifest.scm
```

Inside the environment you can then perform the usual build steps:

```sh
./bootstrap.sh # to generate the "configure" script
./configure
make
```
