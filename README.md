# FixedPoint support class

[![Ubuntu](https://github.com/Galfurian/fixedpoint/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/Galfurian/fixedpoint/actions/workflows/ubuntu.yml)
[![Windows](https://github.com/Galfurian/fixedpoint/actions/workflows/windows.yml/badge.svg)](https://github.com/Galfurian/fixedpoint/actions/workflows/windows.yml)
[![MacOS](https://github.com/Galfurian/fixedpoint/actions/workflows/macos.yml/badge.svg)](https://github.com/Galfurian/fixedpoint/actions/workflows/macos.yml)
[![Documentation](https://github.com/Galfurian/fixedpoint/actions/workflows/documentation.yml/badge.svg)](https://github.com/Galfurian/fixedpoint/actions/workflows/documentation.yml)

**Author**: Enrico Fraccaroli

**Date**: July 24 2019

**Language**: C/C++

## Introduction

This repository contains a easy to use class which represents fixed-point 
values and allows to perform operations between fixed-point values.

## Requirements

A C or C++ compiler.

## Compile

To compile, rely on CMake, and peform the following steps:

```C++
mkdir build
cd build
cmake ..
make
```

## Execute

After compling, you can run the example:

```bash
./test_fixed_point
``` 