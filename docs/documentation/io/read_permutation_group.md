---
title: read_permutation_group
---

Reads an [PermutationGroup](../symmetries/permutation_group.md) object from a TOML file.

**Sources** [read.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/io/read.hpp), [read.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/io/read.cpp)

---

## Definition


=== "C++"
	```c++
	OpSum read_permutation_group(FileToml file, std::string tag);
	```
	
=== "Julia"
	```julia
	read_permutation_group(file::FileToml, tag::String)::PermutationGroup
	```

---

## Parameters

| Name | Description                                                                                                   |   |
|:-----|:--------------------------------------------------------------------------------------------------------------|---|
| file | [FileToml](file_toml.md) object from which the [PermutationGroup](../symmetries/permutation_group.md) is read |   |
| tag  | tag which holds the information in the TOML file about the PermutatioGroup                                    |   |

---

## Data format

A  [PermutationGroup](../symmetries/permutation_group.md) can be defined in a TOML file as an integer matrix, where the rows are the integers of the permutation.

A typical specification of a $C_4$ PermutationGroup is shown here:

```toml
Symmetries = [
  [0, 1, 2, 3],
  [1, 2, 3, 0],
  [2, 3, 0, 1],
  [3, 0, 1, 2]
]
```

!!! warning "1-indexing in Julia / 0-indexing in C++"

	To enumerate the sites of a Permutation, we start counting at 1 in Julia and 0 in C++.
	

---

## Usage Example

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:read_permutation_group"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:read_permutation_group"
	```

