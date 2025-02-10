---
title: read_opsum
---

Reads an [OpSum](../operators/opsum.md) object from a TOML file.

**Sources** [read.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/io/read.hpp), [read.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/io/read.cpp)

---

## Definition


=== "C++"
	```c++
	OpSum read_opsum(FileToml file, std::string tag);
	```
	
=== "Julia"
	```julia
	read_opsum(file::FileToml, tag::String)::OpSum
	```

---

## Parameters

| Name | Description                                                                           |   |
|:-----|:--------------------------------------------------------------------------------------|---|
| file | [FileToml](file_toml.md) object from which the [OpSum](../operators/opsum.md) is read |   |
| tag  | tag which holds the information in the TOML file about the OpSum                      |   |

---

## Data format

An [OpSum](../operators/opsum.md) can be defined in a TOML file as a simple list. The entries of the list are themselves also lists, which contain two or more entries:

1. The first entry must be either a string or a real / complex number denoting the coupling constant of the term.
2. The second entry must be a string and denotes the [operator type](../operators/operator_types.md)
3. The (optional) following entries are integer numbers which denote the sites of the [Op](../operators/op.md)

Here is typical example of a OpSum specification in a TOML file:

```toml
Interactions = [
  ['J1', 'SdotS', 0, 1],
  ['J1', 'SdotS', 1, 2],
  ['J1', 'SdotS', 2, 0],
  ['Jchi', 'ScalarChirality', 0, 1, 2],
]
```

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:read_opsum"
	```

=== "Julia"
	```julia
	--8<-- "examples/usage_examples/main.jl:read_opsum"
	```

