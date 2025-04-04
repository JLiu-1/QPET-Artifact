# QPET integration for SZ3/QoZ/HPEZ (pointwise and regional QoI)

## Introduction

This branch includes the codes for SZ3-QPET, QoZ-QPET, and HPEZ-QPET on pointwise and regional QoIs (single-snapshot compression).

For vector QoI preservation on cross-snapshot compression, check the codes at the szfamily_qpet_vec branch.

## Dependencies

Please Install the following dependencies before compiling HPEZ:

* cmake>=3.13
* gcc>=6.0
* SymEngine (https://github.com/symengine/symengine)
* GMP (https://gmplib.org/) (The dependency of SymEngine)
* Zstd >= 1.3.5 (https://facebook.github.io/zstd/). It is not mandatory to manually install as Zstandard v1.4.5 is included and will be used if libzstd can not be found by pkg-config.

## Installation

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include. A Cmake version >= 3.13.0 is needed. 
Before you proceed to the following evaluations, you may add the installation path of HPEZ to your system path so that you can directly run the **hpez** command in your machine for further evaluation. 
Otherwise, regard the **hpez** command name in the following as **[INSTALL_DIR]/bin/hpez** (the path of your installed executable).

## HPEZ (base) Compression/Decompression Examples

This section is irrelevant to QoI-preserving compression but provides fundamental instructions for how to run several error-bounded lossy compressors in the SZ family. 

The **hpez** command integrates 5 different compression levels by the argument -q, corresponding to 3 compressors (they share the same argument list):

* hpez -q 0: **SZ3.1** compression.
* hpez -q 1: **QoZ 1.1** compression.
* hpez -q 2/3/4: different optimization levels of **HPEZ** compression (level 3 recommended, which is the default).

In the terminal, just run **hpez** or **hpez -h** to see the detailed usage. For validation tasks, the most convenient way is to run the compression, decompression, and data validation in a single command:

**hpez -q [level] -f/-d -a -[Dim_num] [fastest_dim_size] [second_fastest_dim_size] [slowest_dim_size] -i [input_file_name] -o [output_file_name] -M REL [error_bound (e.g. 1e-3)]**

-f: single-precision floating point data, -d: double=precision floating point data. -a: Validate decompression data quality. -M REL: value-range-based error bound. The input error bound will be multiplied by the data range.

The input file should be a binary file of data array. A 100x200x300 3D data array should use the dimensional arguments as **-3 300 200 100.**

More examples are shown in the output of **hpez -h**.

## QPET QoI-preserving HPEZ Compression/Decompression Examples

The QPET functionalities in the **hpez** command can just be activated by specifying additional arguments (including a config file).

Command Line arguments: **-m REL [relative QoI error tolerance]** or **-m ABS [absolute QoI error tolerance]**

Config file: check the  **qoi_configs** folder to find templates, instructions, and examples.

Example of QPET-Integrated HPEZ compression: **[HPEZ base command] -m REL 1E-3 -c qoi.config**

## QoI validation tool

This branch contains another executable, **[INSTALL_DIR]/bin/qoi_val**. **qoi_val** can evaluate the QoI errors between any original data and decompressed data (they need to share the same data type and shape).

Usage: **qoi_val -f/-d -3 dim3 dim2 dim1 -i [original_file] -o [decompressed_file] -c qoi.config**

The QoI to be evaluated should be described in the qoi.config file.

## Test Dataset

4 evaluated datasets in the paper (Miranda, NYX, Scale, Hurricane) can be accessed at [SDRBench](https://sdrbench.github.io/). For Miranda, we converted it to float32 before the evaluation (the original data is double). For Hurricane, we didn't use the logarithmic fields.

The other 2 datasets (RTM, SegSalt) are not publicly accessible due to their commercial source.
