---
title: read_representation
---

Reads an [Representation](../symmetries/representation.md) object from a TOML file.

**Sources** [read.hpp](https://github.com/awietek/xdiag/blob/main/xdiag/io/read.hpp), [read.cpp](https://github.com/awietek/xdiag/blob/main/xdiag/io/read.cpp)

---

## Definition

=== "C++"
	```c++
	OpSum read_representation(FileToml file, std::string irrep_tag, 
	                          std::string group_tag = "Symmetries");
	```
	
=== "Julia"
	```julia
	read_representation(file::FileToml, irrep_tag::String, 
	                    group_tag::String = "Symmetries")::Representation
	```

---

## Parameters

| Name      | Description                                                                                                       |              |
|:----------|:------------------------------------------------------------------------------------------------------------------|--------------|
| file      | [FileToml](file_toml.md) object from which the [Representation](../symmetries/representation.md) is read          |              |
| irrep_tag | tag which holds the information about the [Representation](../symmetries/representation.md) in the TOML file      |              |
| group_tag | tag which holds the information about the [PermutationGroup](../symmetries/permutation_group.md) in the TOML file | "Symmetries" |

---

## Data format

A  [Representation](../symmetries/representation.md) can be defined in a TOML file by specifying up two or three things:

1. The [PermutationGroup](../symmetries/permutation_group.md) as an integer matrix
2. The `characters` of the representation 
3. (optional) the `allowed_symmetries` of the representation, i.e. a list of the number of symmetries used in the irrep. By default all symmetries of the group are used in the representation.

A typical specification of several Representations is shown here:

```toml
	--8<-- "misc/data/irreps.toml"
```

!!! warning "1-indexing in Julia / 0-indexing in C++"

	To enumerate the sites of a Permutation, we start counting at 1 in Julia and 0 in C++.
	

---

## Usage Example

The example reads the representation defined in the `irreps.toml` file, whose contents are shown up in the section [Data format](#data-format).

=== "C++"
	```c++
	--8<-- "examples/usage_examples/main.cpp:read_representation"
	```

=== "Julia"
	```c++
	--8<-- "examples/usage_examples/main.jl:read_representation"
	```

