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

To fetch code that is common to all PiGx pipelines run this:

```sh
git submodule update --init
```

Inside the environment you can then perform the usual build steps:

```sh
./bootstrap.sh # to generate the "configure" script
./configure
make
make check
```

Please, also see https://github.com/BIMSBbioinfo/pigx_sarscov2_ww/blob/main/documentation/user_installation_doc.md for setting up the pipeline and the needed databases.
