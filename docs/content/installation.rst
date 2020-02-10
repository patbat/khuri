Installation
============

Requirements
------------

To use this software, you need to have

* `cmake <https://cmake.org/>`_ (version 3.5 or newer)
* a compiler that supports c++17
* the `GNU scientific library <https://www.gnu.org/software/gsl/>`_ (version
  2.1 or newer)
* the `Eigen library <http://eigen.tuxfamily.org>`_
* python (version 3.7.1 or newer) as well as
  standard python packages (numpy, scipy, setuptools, pytest), as
  provided e.g. via a standard installation of
  `anaconda <https://www.anaconda.com/distribution/>`_

Installation of khuri
---------------------

To install this package, enter the following commands into a shell (without
switching directories in between)::

    git clone --recursive https://github.com/patbat/khuri.git
    pip install ./khuri

To check whether the install was successful, try the following::

   cd khuri/khuri/tests/
   pytest -v

If everything works fine, several tests should run successfully.

To build the documentation for offline access (optional), go to ``khuri/docs``.
Subsequently enter::

   make html

Then open the file ``khuri/docs/_build/html/index.html`` in a web browser.
