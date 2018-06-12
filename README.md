# RETRO
( **R**adio n**E**u**TR**ino simulati**O**n )

## Description

Coming soon ...

## Building the binaries

Clone the present repository and run the [build](build.sh) script, e.g:

```bash
git clone https://github.com/grand-mother/retro
cd retro
./build.sh
```

Note that you might need to install some external dependencies first, e.g.
`gfortran` for TAUOLA. In addition, for a detailed topography you'll need to
fetch the [SRTMGL1](https://lpdaac.usgs.gov/node/527) tiles corresponding to
your location.

## Documentation

In order to run a simulation you'll need a configuration card,
e.g. [share/cards/ulastai.json](share/cards/ulastai.json). Then run:

```bash
bin/retro share/cards/ulastai.json
```

## License

The RETRO suite is under the **GNU LGPLv3** license. See the provided
[LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files.
