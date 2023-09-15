Code Reference
==============

Here you will find reference to the specific functions used in the RAFIKI pacakge. These are found in two python files: Generating_SZ_Data.py which contains the methods for 
cutting particles, and creating a projected map of the tSZ signal in your box, and functions.py which contains all the functions relevant for analysis. Below we summarize 
the role of each function as well as its inputs and outputs. 

..autofunction:: Generating_SZ_Data.determining_frb_size

.. autofunction:: Generating_SZ_Data.generating_sz_data

..autofunction:: Generating_SZ_Data.determinging_caesar_conversion

.. autofunction:: functions.sorting

.. autofunction:: functions.cropping_sz_x

.. autofunction:: functions.cropping_sz_y

.. autofunction:: functions.cropping_sz_z	

.. autofunction:: functions.azimuthalAverage 

.. autofunction:: functions.make_radial_profiles

The next set of functions are used when generating radial profiles of the moments

.. autofunction:: functions.topolar

.. autofunction:: functions.make_moment_profiles


The next function is used when calculating thermal energy

.. autofunction:: functions.thermal_energy


The next set of functions are used when calculating error via the bootstrapping method

.. autofunction:: functions.gen_random_indices

.. autofunction:: functions.single_catalog_bootstrap


The last function is used for plotting and visualizing data

.. autofunction:: functions.make_thermal_energy_plot

