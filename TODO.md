# to do list

## most urgent task

Generate python interface for Omnes function via pybind11.
There are two short projects showing how to combine cmake and setuptools to
be able to install the code nicely via pip, namely
1. project 1 at `https://www.benjack.io/2018/02/02/python-cpp-revisited.html`;
2. project 2 at `https://github.com/pybind/cmake_example`.
The `setup.py` of project 1 is basically a slightly modified version of the one
of project 2 (at least before project 1 introduces tests, both scripts are
essentially the same). However, be aware that project 1 uses deprecated
pybind11 features (this does not affect `setup.py`). The `setup.py` seems to
be sufficiently general.

Proceed along the following lines:
1. adapt directory structure from project 1, because it seems to work fine
   and is simple
2. add pybind11 as a submodule (hence, one needs to use `git clone --recursive`)
3. include `setup.py` from projects 1/2 (first: without tests)
4. write interface for Omnes function and check if it compiles properly
5. include python tests (also via adaption of `setup.py` as far as necessary),
   see if this requires a modified directory structure
6. make travis-ci work again (including both cpp and python tests)

## further tasks

* get documentation going, incorporate doxygen and python docstrings
* add install instructions including a list of software dependencies
