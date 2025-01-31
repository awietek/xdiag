---
title: FileToml
---

A file handler for TOML files. Most funtionality is only provided for the C++ version as Julia already provides good tools handling TOML files.

**Sources** [file_toml.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/io/file_toml.hpp), [file_toml.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/io/file_toml.cpp)

---

## Constructors

=== "C++"	
	```c++
	FileToml(const char *filename);
	FileToml(std::string filename);
	```
=== "Julia"
	```julia
    FileToml(filename::String)
	```
	
| Name     | Description               | Default |
|:---------|:--------------------------|---------|
| filename | filename of the TOML file |         |

---

## Methods

#### defined

Returns whether or not the TOML file has a certain key defined.

=== "C++"	
	```c++
	bool defined(FileToml const &fl, std::string key);
	```
		
---

#### getindex / operator[] (C++ only)

Returns a handler to a value to be read from the TOML file.

=== "C++"	
	```c++
	io::FileTomlHandler operator[](std::string key);
	```
	

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:FileToml"
	```
