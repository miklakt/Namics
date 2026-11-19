from setuptools import Extension, setup

setup(ext_modules=[Extension("namics._namics_kernels", ["cpp/propagate.cpp"],
                             language="c++", extra_compile_args=["-std=c++17", "-O3"])])
