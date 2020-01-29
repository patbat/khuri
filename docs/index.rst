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

.. contents::
   :local:

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

Conventions
===========

Throughout the code, the following conventions apply:

* parameters like thresholds etc. are given in units of Mandelstam variables
* whenever an explicit choice of units is necessary (e.g. for the default
  parameters of phase parametrizations by the Madrid group), (powers of) GeV
  are used

First steps
===========

Try the following inside a python script (or interactive session) to produce a
plot of the Omnes function of the Madrid p-wave:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt

    from khuri import madrid_global, phases, omnes


    THRESHOLD = (2.0 * madrid_global.PION_MASS)**2


    @phases.asymptotic1(matching_point=1.12**2)
    def phase(s):
        return madrid_global.p_wave_phase(s)


    omnes_function = omnes.generate_omnes(phase, threshold=THRESHOLD,
                                          constant=np.pi, cut=1e10)

    energies = np.linspace(0, 1.2, 200)
    omnes_values = omnes_function(energies**2)

    plt.title('The Omnes function of the Madrid p-wave')
    plt.plot(energies, np.real(omnes_values), label='Re')
    plt.plot(energies, np.imag(omnes_values), label='Im')
    plt.xlabel('E/GeV')
    plt.legend()
    plt.show()
