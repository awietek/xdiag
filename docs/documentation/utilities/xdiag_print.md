# Debug printing

For quick debugging in C++, XDiag features a simple macro which outputs the name and content of a variable calles `XDIAG_PRINT(x)`. For example

```c++
Spinhalf block(16, 8);
XDIAG_PRINT(block);
```

will write an output similar to

```text
block:
  n_sites  : 16
  n_up     : 8
  dimension: 12,870
  ID       : 0xa9127434d66b9878
```

The `XDIAG_PRINT(x)` macro can be used on any XDiag object and several other standard C++ objects as well.
