# Timing

In standard C++ measuring time is a bit awkward. To quickly monitor the CPU time spent by XDiag by simple functions.

## Simple timing using tic() / toc()
Similar as in Matlab one can use `tic()` and `toc()` to measure the time spent between two points in the code. 

```c++
tic();
double e0 = eigval0(bonds, block);
toc();
```

`toc()` will output the time spent since the last time `tic()` has been called.

## Detailed timing

To get the present time, simply call 

```c++
auto time = rightnow();
```

A timing (in second) between two time points can be written to output using

```c++
timing(begin, end);
```

This can even be accompanied by a message about what is being timed and a verbosity level (see [Logging](logging.md)) can also be set. The full call signature is

```c++
timing(begin, end, message, level);
```

| Name    | Description                                | Default |
|:--------|:-------------------------------------------|---------|
| begin   | starting time computed using `rightnow()`  |         |
| end     | end time computed using `rightnow()`       |         |
| message | message string to be prepended to timing   | ""      |
| level   | verbosity level at which timing is printed | 0       |
