.. RAFIKI documentation master file, created by
   sphinx-quickstart on Tue Jul 11 10:36:06 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RAFIKI's documentation!
==================================
.. image:: images/rafiki.png
   :width: 400

RAFIKI (Realizing AGN Feedback in Keyframe Images) is an open-source python package developed to analyze
simulated microwave and X-ray data. In its current form, it can be used to generate and stack tSZ images
from raw hdf5 files from SIMBA or any other GIZMO-based code.

.. note:: 
	This project is under active development. The current version is only functional for tSZ analysis

Features
--------
        -Create stacked images of Compton-y parameter around galaxies with properties of your choosings
       
	-Generate radial profiles of the Compton-y parameter and moment analyses of asymmetries
        
	-Calculate the total thermal energy in an aperature of your choosing and plot against stellar or
        halo mass

Use
---
At it's base, RAFIKI's functionality is contained in three steps:
        #. Generate maps of tSZ signal from raw hdf5 files
        
	#. Extract data around galaxies of interest in your map
        
	#. Stack and visualize your data by comparing a range of models and redshifts

Each of these steps is run individually with several user input parameters.
See the :doc:`usage` section for further information on requirements and :doc:`running` section for information on how each of these steps is executed. 



Support
-------
If you have questions, contact Skylar Grayson sigrayso@asu.edu



Indices and tables
==================
.. toctree::
   usage
   running
   outputs
   reference

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
