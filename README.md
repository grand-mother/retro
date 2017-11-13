# RETRO
( **R**adio n**E**u**TR**ino simulati**O**n )

## Description

Coming soon ...

## Installation

Clone the present repository and run the [installation](install.sh) script, e.g:

```bash
git clone https://github.com/grand-mother/retro
cd retro
./install.sh
```

Note that you might need to install some external dependencies first, e.g.
`libtiff` for reading the [ASTER-GDEM2](https://asterweb.jpl.nasa.gov/gdem.asp)
topography data with TURTLE.

## Documentation

In order to run a simulation you'll need a configuration card,
e.g. [share/cards/ulastai.json](share/cards/ulastai.json). Then run:

```bash
. setup.sh
retro-run share/cards/ulastai.json
```

Note that you'll need to _initialise you environment only once_, i.e. `. setup.sh`.

## License

The RETRO suite is under the **GNU LGPLv3** license. See the provided
[LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files.
