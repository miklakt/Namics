"""Load a C ABI kernel that operates on caller-owned NumPy buffers."""
from __future__ import annotations

import ctypes
from functools import lru_cache
from importlib.util import find_spec
import numpy as np



@lru_cache(maxsize=1)
def load_kernel():
    spec = find_spec("namics._namics_kernels")
    if spec is None or spec.origin is None:
        raise RuntimeError("C++ kernel not built; install the package or run make in python-next")
    path = spec.origin
    lib = ctypes.CDLL(str(path))
    array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
    fn = lib.namics_linear_density
    indices = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS")
    fn.argtypes = [array, indices, indices, array, ctypes.c_int, ctypes.c_int,
                   ctypes.c_int, ctypes.c_int, ctypes.c_double,
                   array, array, array, array]
    fn.restype = ctypes.c_int
    return fn
