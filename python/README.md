# Packed Compressed Sparse Row for Python

This README will guide you through how to use the PCSR implementation using Python.

## Compiling PCSR.cpp into the .so file

Run the following to get the shared object file, that can be pip installed to be used inside any python program

```bash
c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) PCSR.cpp -o PCSR$(python3-config --extension-suffix)
```