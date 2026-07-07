* Compile the commands using 

```bash
cmake -B build -S . -D XDIAG_DISABLE_HDF5=On -D CMAKE_EXPORT_COMPILE_COMMANDS=On -D CMAKE_CXX_COMPILER=clang++
```

Comment: clang++ is needed for OpenMP to be detected on OSX
This generates the file `compile_commands.json` in the `build` directory.;

* then first dump the API
```bash
./dump_api.py
```

* then create wrapper
```bash
./create_wrapper_cpp.sh
```

* then install wrapper
```bash
./install_wrapper_cpp.sh
```

