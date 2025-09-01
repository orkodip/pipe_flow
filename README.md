# pipe_flow
This repository contain the codes (written in C++) to compute the fluid flow inside a circular pipe. Finite Volume Method is used to discretize the governing equations in axisymmetric cylindrical coordinates.

# Software requirements
This solver needs:

- gcc

# How to install the required packages (on a Linux system)

To install gcc (compiler for C++ code)

```bash
sudo apt install build-essential
```

# How to compile and run the code

To compile the code

```bash
g++ driver.cpp -o output
```
To run this code

```bash
./output
```
