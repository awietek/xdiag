---
title: FileH5
---

A file handler for [hdf5](https://www.hdfgroup.org/solutions/hdf5/) files. The proper tool to write results of XDiag simulations to disk. Only provided for the C++ version as Julia already provides good tools handling hdf5 files.

**Sources** [file_h5.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/io/file_h5.hpp), [file_h5.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/io/file_h5.cpp)

---

## Constructors

=== "C++"	
	```c++
	FileH5(std::string filename, std::string iomode = "r");
	```

	
| Name     | Description                      | Default |
|:---------|:---------------------------------|---------|
| filename | filename of the hdf5 file        |         |
| iomode   | whether to read or write to file | `r`     |

There are four possible values of `iomode`

1. `r`: read-only mode
2. `w`: secure write mode, new file is created if it does not exist
3. `w!`: forced write mode, existing file will be overwritten
4. `a`: append mode, new datasets in an existing file can be created

---

## Methods


#### getindex / operator[]

Returns a handler to a value to be read or written from or to the hdf5 file.

=== "C++"	
	```c++
	hdf5::FileH5Handler operator[](std::string key);
	```
	

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:FileH5"
	```
