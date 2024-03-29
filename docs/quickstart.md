---
title: Quick start
---

## Hello World

Let us set up our first program using the `hydra` library. For this we need to create two files. The first is the actual `C++` code.

```C++
#include <hydra/all.h>

using namespace hydra;

int main() try {
  say_hello();
} catch (std::exception const &e) {
  traceback(e);
}
```

The function `say_hello()` prints out a welcome message, which also contains information which exact hydra version is used. What is maybe a bit unfamiliar is the `try / catch` block. Hydra implements a traceback mechanism for runtime errors, which is activated by this idiom. While not stricly necessary here, it is a good practice to make use of this.

Now that the application program is written, we next need to set up the compilation instructions using [CMake](https://cmake.org/). To do so we create a second file called `CMakeLists.txt` in the same directory.

```cmake
cmake_minimum_required(VERSION 3.19)

project(
  hello_world
  LANGUAGES CXX
)

find_package(hydra REQUIRED HINTS "/path/to/hydra/install")
add_executable(main main.cpp)
target_link_libraries(main PUBLIC hydra::hydra)
```

You should replace `"/path/to/hydra/install"` with the appropriate directory where your hydra library is installed after compilation. This exact `CMakeLists.txt` file can be used to compile any hydra application.

We then compile the application code,

```bash
cmake -S . -B build
cmake --build build
```

and finally run our first `hydra` application.

```bash
./build/main
```