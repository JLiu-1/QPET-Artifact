# QPET integration for SZ3/QoZ/HPEZ (3D vector QoI)

## Introduction

This branch includes the codes for SZ3-QPET, QoZ-QPET, and HPEZ-QPET on 3D vector QoIs (three-snapshot compression).

For other QoIs, check the codes at the szfamily_qpet branch.

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
Before you proceed to the following evaluations, you may add the installation path to your system path so that you can directly run the **hpez_vec** command in your machine for further evaluation. 
Otherwise, regard the **hpez** command name in the following as **[INSTALL_DIR]/bin/hpez_vec** (the path of your installed executable).


## QPET QoI-preserving HPEZ-Vec Compression/Decompression Examples

The **hpez_vec** command can do both normal compression and QPET-integrated compression you can still apply **hpez_vec -q 0/1/3** to perform SZ3/QoZ/HPEZ vector compression. The arguments have similar structure with the **hpez** command, but some of them now need 3 values;  

Base command: **hpez_vec -q [level] -f/-d -a -[Dim_num] [fastest_dim_size] [second_fastest_dim_size] [slowest_dim_size] -i [input_file_name1] [input_file_name2] [input_file_name3] -o [output_file_name1] [output_file_name2] [output_file_name3] [output_file_name1] -M REL [error_bound (e.g. 1e-3)]** (the error bound is shared among three snapshots).

QoI Config file: check the  **qoi_configs** folder to find templates, instructions, and examples.

QPET-Integrated HPEZ compression: **[HPEZ base command] -m REL [QoI error threshold] -c qoi.config**

## Vector QoI validation tool

This branch contains another executable, **[INSTALL_DIR]/bin/qoi_val_vec**. **qoi_val_vec** can evaluate the vector QoI errors between any original data and decompressed data (they need to share the same data type and shape).

Usage: **qoi_val -f/-d -3 dim3 dim2 dim1 -i [original_file1] [original_file2] [original_file3] -o [decompressed_file1] [decompressed_file2] [decompressed_file3] -c qoi.config**

The QoI to be evaluated should be described in the qoi.config file.


## Test Dataset

The evaluated datasets in the paper (Miranda, NYX, Scale, Hurricane) can be accessed at [SDRBench](https://sdrbench.github.io/). For Miranda, we converted it to float32 before the evaluation (the original data is double).

Use the velocity fields (vx, vy, vz or u, v, w) to evaluate vector compression.


