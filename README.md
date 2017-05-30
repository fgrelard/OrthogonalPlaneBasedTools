# MyDGtalContrib

Description
===========
Algorithms implemented thanks to the [DGtal Library](http://dgtal.org/)

On tubular volumes:

* Orthogonal plane estimation
* Pruning algorithm using orthogonal planes
* Recentering procedure
* Curve-skeleton extraction



Quick Build Instructions
========================
The main instructions on linux/unix-based systems are the following:

```shell
git clone https://github.com/fgrelard/MyDGtalContrib.git
cd MyDGtalContrib ; mkdir build ; cd build
cmake ..
make
```

Minimum system requirements: C++11 enabled compiler, [cmake](http://cmake.org), [DGtal](http://dgtal.org/), [boost](http://boost.org) (>= 1.46), [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (>=3.2.0)

In order to compile TESTS and EXAMPLES, the following additional libraries are necessary:
* [QGLViewer](http://libqglviewer.com/) (>=2.5.0)
* Optionally, [ITK](https://itk.org/)

DGtal needs to be compiled with these libraries as well (checkout WITH_QGLVIEWER, WITH_EIGEN, and WITH_ITK options with CMake)
