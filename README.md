This is the implementation of 

[Erik-Oliver Blass, Florian Kerschbaum, “BOREALIS: Building Block for Sealed Bid Auctions on Blockchains”, Proceedings of ACM Asia Conference on Computer and Communications Security (AsiaCCS), Taipei, Taiwan, 2020](https://eprint.iacr.org/2019/276.pdf).


The code in this repository is a fork of https://github.com/baz1/GS-NIZK

See below for installation requirements. If your platform meets all requirements, just run `make` and then `lib/gsnizk/gsnizk`. The relevant proofs are in file `lib/gsnizk/tests.cpp`. 


Below is the original README.
***
# GS-NIZK

Efficient implementation of the [Non-Interactive Zero-Knowledge proofs
invented by Groth and Sahai][1]

## Usage

Even though the source code has been designed to be cross-platform compatible,
the Makefile script has only been written for a Linux-64 based platform.
For now, the GS-NIZK project itself is managed with Qt:
you will find the project file `gsnizk.pro` under `lib/gsnizk/`.
Once the settings at the beginning of that file have been set,
you can compile the project with the `Makefile`, as described below
(be sure to be consistent with your choices).
It will need the tool `qmake` as well as `g++` to compile the main project.

### Compiling the library for MIRACL

```
$ make miracl-build
```

### Compiling the library for PBC

First of all, make sure that the following packages are installed
(else, run the corresponding installation line)
```
# apt-get install libgmp3-dev
# apt-get install flex
# apt-get install bison
```
Then, you can start compiling the PBC library with:
```
$ make pbc-build
```
This compiles the library and creates a static library file.
If you want to install PBC on your computer, run this instead:
```
# make pbc-build-install
```

### Creating the documentation

```
$ make doc
```
This will create a folder `doc` containing the documentation for
GS-NIZK's library. Note that you need to have Doxygen (`doxygen`) installed
to do this, as well as the `dot` tool.

### Doing all of this at once

```
$ make
```
(see notes for the PBC build and the documentation)

[1]: https://eprint.iacr.org/2007/155
