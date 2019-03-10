# Narrow Band FLIP using APIC Particles

This repo is based in Mantaflow Project.

## Installation

1. Create a `build` folder and enter with

```sh
mkdir build && cd build
```

2. Run cmake

```sh
cmake .. -DOPENMP=ON -DGUI=ON
```

3. Compile

```sh
make
```

## Run the project

There's two files in the folder `scenes`. The `flip_nbofficial.py` runs the initial paper presented in SIGGRAPH 2016, the other one contains the implementaion using APIC particles.

```sh
./manta ../scenes/flip_josue.py
```
