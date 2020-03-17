# lattice-c
[Lattice](https://github.com/rjb3977/Lattice), but in C. And way better. Requires [gmp](https://gmplib.org/) to link, and a recent C compiler to compile. I personally use [clang-10](https://clang.llvm.org/).

This implementation provides a ~28x improvement in speed on an 12-core, 24-thread AMD Ryzen 3900X. There are likely still improvements to be made to the Simplex implementation.
