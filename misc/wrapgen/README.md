* Steps:

1. Remove previous build and install directories
2. Run `configure.sh`: This will set up the build of the wrapper and extract `compile_commands.json`
3. Run `generate.sh`: This will create the actual wrapper code
4. Run `build.sh`: This will now compile the C++ part of the wrapper
5. Run `install.sh`: This will copy the wrapper code where it belongs and installs the libxdiagjl library
