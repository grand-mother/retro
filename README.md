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
`libtiff` for reading the [ASTER-GDEM2](https://asterweb.jpl.nasa.gov/gdem.asp)
topography data with TURTLE.

## Documentation

In order to run a simulation you'll need a configuration card,
e.g. [share/cards/ulastai.json](share/cards/ulastai.json). Then run:

```bash
bin/retro share/cards/ulastai.json
```

## License

The RETRO suite is under the **GNU LGPLv3** license. See the provided
[LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files.
