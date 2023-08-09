# Adjoint Accelerator Optimization (PyAAO)

## Install 

Use pip to install the required python packages:

```
python3 -m pip install -r requirements.txt
```

C/C++ compiled libraries (.dll and .so)  already exist under the bindings/ folder, so no need to compile. If you want to recompile:

1. Windows: open the bindings/src/AdjointFTR folder in visual studio as a cmake project and build the project, this will produce a new .dll

2. Linux: Make sure you have cmake and gcc installed, then follow these steps to generate a .so:
```
cd bindings/src/AdjointFTR
mkdir build
cd build
cmake ..
make -j4
```

## Running

Use the pyaao/ folder for example run and optimizations:

```
python3 pyaao/bbcMatchRun.py
```

See the systems/ folder for example lattices and optimization parameter setups

