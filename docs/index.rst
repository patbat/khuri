.. khuri documentation master file, created by
   sphinx-quickstart on Thu Sep 19 21:06:37 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====
khuri
=====

.. toctree::
   :maxdepth: 2
   :caption: Contents:

The aim of this project is to develop a general solver for Khuri-Treiman
equations including facilities that are of general use when dealing with
dispersion relations, e.g. Omnes function(s).


Installation
============

Requirements
------------

To use this software, you need to have

* `cmake <https://cmake.org/>`_ (version 3.5 or newer)
* a compiler that supports c++14
* the `GNU scientific library <https://www.gnu.org/software/gsl/>`_ (version
  2.1 or newer)
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

(You may remove the source code directory afterwards).

First steps
-----------

To check whether the install was successful, try the following inside a
python script (or interactive session) to produce a plot of the Omnes function
of the Madrid p-wave:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt

    from khuri.tests.test_omnes import omnes_function

    o = omnes_function()

    energies = np.linspace(0, 1200, 200)
    omnes_values = o(energies**2)

    plt.title('The Omnes function of the Madrid p-wave')
    plt.plot(energies, np.real(omnes_values), label='Re')
    plt.plot(energies, np.imag(omnes_values), label='Im')
    plt.xlabel("E/MeV")
    plt.legend()
    plt.show()
