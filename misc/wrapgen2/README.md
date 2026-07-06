* Compile the commands using 

```bash
cmake -B build -S . -D XDIAG_DISABLE_HDF5=On -D CMAKE_EXPORT_COMPILE_COMMANDS=On -D CMAKE_CXX_COMPILER=clang++
```

Comment: clang++ is needed for OpenMP to be detected on OSX
This generates the file `compile_commands.json` in the `build` directory.;

* Compile wrapper

```bash
cmake -S . -B build -D XDIAG_JULIA_WRAPPER=On -D CMAKE_PREFIX_PATH=/Users/awietek/.julia/artifacts/44722d868723dcb2366a9e90f02b34885cd80aba -D CMAKE_CXX_COMPILER=clang++

```
