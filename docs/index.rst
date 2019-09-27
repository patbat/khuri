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

* cmake (version 3.5 or higher)
* a compiler that supports c++14
* python (version 3.7.1 or higher) as well as
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
python script (or interactive session):

.. code-block:: python

   from khuri.madrid import phase

   # the madrid phase (p-wave) evaluated at the mass of the rho meson
   print(phase(760.0**2))
