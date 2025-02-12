# BlipSim
Simulation code for Blip


## Installation instructions

### Prerequisites
- ROOT
- GEANT4.11 w/ Multithreading, GDML, and QT installation ON.
- cmake

#### Recommended
- ccmake (install using `sudo apt-get install cmake-curses-gui`

### Installation
```bash
git clone https://github.com/infophysics/BlipSim.git
cd BlipSim
git submodule init
git submodule update --recursive
mkdir build; cd build/
ccmake ../
make -j16
