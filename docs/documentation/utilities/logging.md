# Logging

## Setting the verbosity

Algorithms implemented in XDiag do not output anything during their execution by default. However, it is typically useful to get some information on how the code is performing and even intermediary results at runtime. For this, the verbosity of the internal XDiag logging can be set using the function `set_verbosity`, which is defined as ([logger.hpp](https://github.com/awietek/xdiag/blob/master/xdiag/utils/logger.hpp))

=== "C++"
	```c++ 
    void set_verbosity(int64_t level);
	```

=== "Julia"

	```julia 
	set_verbosity(level::Integer);
	```
	
There are several levels of verbosity, defining how much information is shown.

| level | outputed information |
|:------|:---------------------|
| 0     | no information       |
| 1     | some information     |
| 2     | detailed information |

For example, when computing a ground state energy using the [eigval0](algorithms/eigval0.md) function, we can set a higher verbosity level using

=== "C++"
	```c++ 
    set_verbosity(2);
	double e0 = eigval0(bonds, block);
	```

=== "Julia"

	```julia 
    set_verbosity(2);
    e0 = eigval0(bonds, block);
	```
This will print detailed information, which can look like this

```text
Lanczos iteration 1
MVM: 0.00289 secs
alpha: -0.2756971549998545
beta: 1.7639347562074059
eigs: -0.2756971549998545
Lanczos iteration 2
MVM: 0.00244 secs
alpha: -0.7116140394927443
beta: 2.3044797637130743
eigs: -2.2710052270892791 1.2836940325966804
Lanczos iteration 3
MVM: 0.00210 secs
alpha: -1.2772539678430306
beta: 2.6627870395174456
eigs: -3.7522788386927637 -0.6474957945455240 2.1352094709026579
```

## Log mechanism (C++ only)
Producing nicely formatted output is unfortunately a bit cumbersome in standard C++. For this, the `Log` mechanism in XDiag can help. To simply write out a line of information you can call,

```c++
Log("hello from the logger");
```

By default, a new line is added. It is also possible to set verbosity by handing the level as the first argument,

```c++
Log(2, "hello from the logger only if global verbosity is set to >= 2");
```

This message will only appear if the global verbosity level is set to a value $\geq 2$. Finally, XDiag also supports formatted output by using the [fmtlib](library). For example, numbers can be formated this way

```c++
Log("pi is around {:.4f} and the answer is {}", 3.141592, 42);
```
