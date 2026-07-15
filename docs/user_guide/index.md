---
title: User Guide
---

# User Guide

This user guide is a hands-on, step-by-step introduction to XDiag that takes you from a fresh installation all the way to running your own exact diagonalization calculations. Each section introduces one core concept and demonstrates it with parallel **Julia** and **C++** code examples, so you can follow along in whichever language you prefer.

The guide is organized to mirror a typical workflow. We start with [installation](01-installation.md) and [setting up an application](02-application.md), then cover the fundamental building blocks of a quantum many-body problem — [Hilbert spaces](03-hilbert-spaces.md), [operators](04-operators.md), and [states](05-states.md). Building on these, we turn to the numerical linear algebra at the heart of XDiag: [dense](06-dense.md) and [iterative](07-diagonalization.md) diagonalization, [time evolution](08-time-evolution.md), and [measurements](09-measurements.md). The remaining sections deal with [input/output](10-io.md), the use of [symmetries](11-symmetries.md), and advanced topics such as [sparse matrices](12-sparse.md) and [distributed-memory computing](13-distributed.md).

For a complete reference of every function and object, see the [Documentation](../documentation/index.md); for fully worked applications, see the [Examples](../examples.md).

## Getting started
<div class="grid cards" markdown>

-   :material-download: **Installation**

    ---

    Installation and compilation of the XDiag core library

    [:octicons-arrow-right-24: Read more](01-installation.md)

-   :material-cog: **Application code**

    ---

    How to set up your application code using Julia or C++

    [:octicons-arrow-right-24: Read more](02-application.md)

</div>


## Basic quantum mechanics

<div class="grid cards" markdown>

-   :material-react: **Hilbert spaces**

    ---

    How to define a Hilbert space of fermions, bosons, or spins. 

    [:octicons-arrow-right-24: Read more](03-hilbert-spaces.md)

-   :material-react: **Operators**

    ---

    How to define operators and perform algebraic operations with them.

    [:octicons-arrow-right-24: Read more](04-operators.md)
	
-   :material-react: **States**

    ---

    How to define a quantum state and act on it with operators

    [:octicons-arrow-right-24: Read more](05-states.md)
		
</div>

## Linear Algebra

<div class="grid cards" markdown>

-   :material-matrix: **Dense matrices**

    ---

    How to compute a dense matrix representation of an operator and perform a full diagonalization 

    [:octicons-arrow-right-24: Read more](06-dense.md)

-   :material-matrix: **Diagonalization**

    ---

    How to perform a diagonalization using iterative numerical algorithms

    [:octicons-arrow-right-24: Read more](07-diagonalization.md)
	
-   :material-matrix: **Time evolution**

    ---

    How to evolve a state in real or imaginary time

    [:octicons-arrow-right-24: Read more](08-time-evolution.md)
	
	
-   :material-matrix: **Measurements**

    ---

    How to compute expectation values of observables 

    [:octicons-arrow-right-24: Read more](09-measurements.md)
	
</div>


## Input / Output

<div class="grid cards" markdown>

-   :material-file: **Input / Output**

    ---

    Reading and writing to TOML and hdf5 files. 

    [:octicons-arrow-right-24: Read more](10-io.md)


</div>

## Symmetries

<div class="grid cards" markdown>

-   :material-octahedron: **Symmetries**

    ---

    How to work with symmetries and symmetry-adapted bases

    [:octicons-arrow-right-24: Read more](11-symmetries.md)


</div>

## Advanced topics

<div class="grid cards" markdown>

-   :material-star: **Sparse matrices**

    ---

    How to create sparse matrices and use them in iterative algorithms

    [:octicons-arrow-right-24: Read more](12-sparse.md)
	
	
-   :material-star: **Distributed computing**

    ---

    How to use the distributed XDiag library for distributed memory computing 

    [:octicons-arrow-right-24: Read more](13-distributed.md)
	
</div>

