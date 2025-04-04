from setuptools import setup, Extension
import numpy

module = Extension(
    "generic",
    sources=["generic.c"],
    include_dirs=[numpy.get_include()]  # <- This line is key
)

setup(
    name="generic",
    version="1.0",
    ext_modules=[module],
)
