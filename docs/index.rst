.. khuri documentation master file, created by
   sphinx-quickstart on Thu Sep 19 21:06:37 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

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

To install this package, enter the following commands (without switching
directories in between)::

    git clone --recursive https://github.com/patbat/khuri.git
    pip install ./khuri

To check whether the install was successful, try the following inside a
python script (or interactive session):

.. code-block:: python

   from khuri.madrid import phases

   # the madrid phase evaluated at the mass of the rho meson
   print(phases(760.0**2))
