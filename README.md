# 2D/3D Elliptic Multiscale Problem with MsFEM in Deal.II (Shared Memory Parallel)

This repository is for teaching purposes only.

It demonstrates how to compue multiscale basis functions in a shared memory parallel fashion in both 2D and 3D.

## Building the executable

To build the executable and Eclipse project files you must clone the repository:

```
git clone https://github.com/konsim83/deal.ii_SS19_demo_elliptic_multiscale_fem.git elliptic_msfem
```
We want an out-of-source-build with build files in a folder parallel to the code:

```
mkdir elliptic_msfem_build
cd elliptic_msfem_build
```
Then create the build files with `cmake`:

```
cmake -DDEAL_II_DIR=/path/to/dealii -G"Eclipse CDT4 - Unix Makefiles" ../elliptic_msfem
```
You can now import an existing project in Eclipse. To generate the executable in debug mode type

```
make debug
make
```
If you want to produce a faster reslease version type

```
make release
make
```

## Building the Documentation

You will need `doxygen`, `mathjax` and some other packages such as `GraphViz` installed.

To build the documentation with `doxygen` enter the code folder

```
cd elliptic_msfem/doc
```
and type

```
doxygen Doxyfile
```
This will generate a html documentation of classes in the `diffusion_equation/documentation/html` directory.
To open it open the `index.html` in a web browser.
