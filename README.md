<img align="right" width=20% src="logo.gif" />

[comment]: <![](./title.svg)>
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
- [HowToRun](#howtorun)
- [Tests](#tests)
- [Demo](#demo)
- [Aknowledgements](#aknowledgements)
- [License](#license)

## Requirements

In order to be able to run everything succesfully you will need:
* A `C++11` compiler
* `armadillo`, an efficient linear algebra library
* `cmake`, for building (optional, recommended)

This project also uses the header-only [`progressbar`](https://github.com/gipert/progressbar) library, which is already included in this repository in `include/progressbar.hpp`.

## Building

Clone this repo with

```bash
git clone https://github.com/martapisci/quantum_slit
```

or

```bash
git clone git@github.com:martapisci/quantum_slit
```

### CMAKE

In order to build everything you need `cmake`. First you need to create a `build/` directory and move into it
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
quantum_slit
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

## HowToRun
There are two possible ways of running the code. You can either modify the `parameters.txt` file and then compile with `cd build; make`  and run `./probability`, OR you can just run the shell-script `./run.sh` form the parent directory. Remember to make it executable, i.e.
```bash
chmod +rx run.sh
```
When running the shell-script you give the parameters as command line arguments (we modify file `parameters.txt` under the hood for you :wink:):
```bash
./run.sh <M> <dT> <T> <xc> <yc> <sigmax> <sigmay> <px> <py> <v0> <nslit>
```
An example would be
```bash
./run.sh 200 0.000025 0.008 0.25 0.5 0.05 0.1 200 0 1e10 3
```
## Tests
CMake provides an easy command for code testing, i.e. from `build/` you can run
```bash
ctest
```
In this way all tests are run.

## Demo
Here you can find a little demonstration of what you can do with our beautiful little program. First of all run
```bash
./run.sh 200 0.000025 0.003 0.25 0.5 0.05 0.2 200 0 0. 0
./run.sh 200 0.000025 0.003 0.25 0.5 0.05 0.2 200 0 1e10 1
./run.sh 200 0.000025 0.003 0.25 0.5 0.05 0.2 200 0 1e10 2
./run.sh 200 0.000025 0.003 0.25 0.5 0.05 0.2 200 0 1e10 3
```
After all of them finished you can move to `plots/` and make the animations with
```bash
python3 animation.py
```
This is going to be the result:
<p align="center">
<button  style="border: transparent; background-color: transparent;">
    <img align="left" width=10% src="demo/no_slit"> 
</button>
<button  style="border: transparent; background-color: transparent;">
    <img align="left" width=10% src="demo/single_slit"> 
</button>
<button style="border: transparent; background-color: transparent;">
    <img align="left" width=10% src="demo/double_slit"> 
</button>
<button style="border: transparent; background-color: transparent;">
    <img align="left" width=10% src="demo/triple_slit"> 
</button>
</p>



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


