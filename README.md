<img align="right" width=20% src="logo.gif" />

# quantum_slit
![License](https://img.shields.io/github/license/martapisci/quantum_slit)
![Language](https://img.shields.io/badge/language-c%2B%2B-blue)
![Size](https://img.shields.io/github/repo-size/martapisci/quantum_slit)

Basic implementation of a solver for the two-dimensional time-dependent Schrödinger equation for a single particle.

Here you can find the code needed for studying the time evolution of the wave function passing through different types of slits in a 2D frame.

This repo is mantained by group 100 of the course FYS3150/4150 at UiO.


## Table of contents

- [Requirements](#requirements)
- [Building](#building)
- [Filesystem](#filesystem)
- [Tests](#tests)
- [Aknowledgements](#aknowledgements)
- [License](#license)

## Requirements

In order to be able to run everything succesfully you will need:
* A `C++11` compiler
* `armadillo`, an efficient linear algebra library
* `cmake`, for building (optional, recommended)

## Building

Clone this repo with

```bash
git clone https://github.com/martapisci/quantum_slit
```

or

```bash
git clone git@github.com:martapisci/quantum_slit
```

### g++

You can compile and link with `g++`.

```bash
g++ main.cpp src/*.cpp -I include -larmadillo -o quantum
```
You might need to add the  `-std=gnu++11` if you are a Mac user.
You can then run the executable with

```bash
./quantum
```

### CMAKE

Alternatively you can build everything with `cmake`. First you need to create a `build/` directory and move into it
```bash
mkdir build; cd build
```
Now you just need to run
```bash
cmake ..
```
If it outputs correctly you should be able to build the project by runnning
```bash
make
```
from the same `build/` directory. 

## Filesystem
The repo is organized as follows:
```
penning_trap
│
└───build/
│   └───data/
│   └───plots/
│
└───include/
│
└───plots/
│
└───src/
│
└───test/
```
Inside the `build/` directory you can find the executables and two subfolders: `build/data/` for storing the resulting data and `build/plots/` for storing the plots of the same resulting data.
In the `include/` directory are stored all the header filese and in the `src/` directory are stored the source files.
In `plots/` you can find also the python scripts that make the graphs and store them in `build/plots/`. For example, to make the plots of `twobody.cpp` simply run
```bash
python3 probability.py
```
from `src/`.
The `test/` directory is meant for the testing of the built code.



## Tests
CMake provides an easy command for code testing, i.e. from `build/` you can run
```bash
ctest
```
In this way all tests are run.


## License

The code here presented is released under version 3 of the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html).


## Aknowledgements
<p align="center">
<button  style="border: transparent; background-color: transparent;">
    <img align="left" width=10% src="https://avatars.githubusercontent.com/u/51904841?v=4"> 
</button>
<button style="border: transparent; background-color: transparent;">
    <img align="left" width=10% src="https://avatars.githubusercontent.com/u/112166702?v="> 
</button>
<button style="border: transparent; background-color: transparent;">
    <img align="left" width=10% src="https://avatars.githubusercontent.com/u/79975678?s=400&u=6770b5f0354ed29bf9a54e7f27a8250bb812c279&v=4"> 
</button>
<button style="border: transparent; background-color: transparent;">
    <img align="left" width=10% src="https://avatars.githubusercontent.com/u/112163092?v=4">
</button>
</p>


