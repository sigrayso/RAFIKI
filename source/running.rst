Running RAFIKI
**************

RAFIKI allows for a range of analyses and is thus run in a series of steps requiring a minimum amount of user input. RAFIKI can be run in command line via its .py 
files or Jupyter notebook. In order to identify galaxies and conduct analysis for different populations, it is necessary that you have a CAESAR output file for each snapshot. Some 
simulations, such as Simba provide these Caesar files along with the snapshots. You can also run CAESAR yourself, following the instructions `here 
<https://caesar.readthedocs.io/en/latest/>`_
 

Once you have your snapshot output (in hdf5 format) and CAESAR galaxy/halo catalog, the process of running RAFIKI is done in two steps:

Generating and Analyzing tSZ Map (see analysis_example.ipynb)

Stacking and Visualizing Data (see stacking_example.ipynb)

Here we walk through the general steps of these processes. For more detailed information on the code structure and functions used, see 

Generating and Analyzing tSZ Map
--------------------------------
There are three steps in the process of analyzing tSZ data, as walked through in analysis_example.ipynb. Those are: 

1. Creating a projected map of the Compton-y parameter
2. Extracting galaxy sample information 
3. Extracting data around galaxies to generate results

There are three main results RAFIKI is equipped to produce: 

1. Radial profiles of the Compton-y signal
2. Radial profiles of moments of asymmetry
3. Thermal energy measurements 

The example notebook shows how to conduct all three analyzes, which generates a series of .csv files which are described in :doc:`outputs` 

There are several key parameters that need to be adjusted, as described briefly in the notebook and more extensively here

**Step One: Creating a Projected Map**

The first stage of analysis is generating a map of the Compton-y parameter projected along the snapshot box. This converts the particle data to a two-dimensional 
pixelated map of the Compton-y parameter, defined as :math:`y = \sigma_T \int dl \ n_e \frac{k(T_e-T_{CMB})}{m_e c^2}` This is done via the generating_sz_data function 
and requires the following inputs: 

.. code-block:: python

	filename = 'snap_m50n512_105.hdf5'
	projection_direction = 'x' 
	output_name = 'noagn_x_szy.npy' 
	z = 0.9927 
	theta = 2.5
	comov =3289.7*1000   #Calculate with https://www.astro.ucla.edu/~wright/CosmoCalc.html
	frb=determining_frb_size(50, z, comov, theta) 

- **filename** - Path and name of snapshot for analysis. Must be .hdf5 format
- **projection_direction** - ‘x’, ‘y’, or ‘z’. It is reccomended that each snapshot is projected in all three directions. Each projection is then treated as a separate 
map with a separate galaxy sample. This is not required, as the statistical legitimacy of this process depends on the radius around the galaxies you are interested in. While the 
tSZ is not expected to be different for different projections while looking within the radius of the galaxy, at larger radii there is no significant spherical symmetry.
- **output_name** - Name of the .npy file containing the projected Compton-y map. This name MUST END with _szy.npy and it is reccomended that the direction of projection 
is included as shown above.
- **z** - redshift corresponding to the snapshot. Deriving the Compton-y parameter from pressure fields as done here is redshift dependent.
- **theta** - Angular resolution you want your frb to have
- **comov** - Comoving distance corresponding to the redshift of your snapshot. Suggested that you calculate this `here <https://www.astro.ucla.edu/~wright/CosmoCalc.html>`_
- **frb** - The number of pixels on each side for a generated fixed resolution buffer. This calculated for you with the ``Generating_SZ_Data.determining_frb_size()`` function. 


**Step Two: Extracting Galaxy Information**

Here, we load the CAESAR file and extract key characteristics of the galaxies, such as stellar mass, halo mass, and age. Aside from the path to the file, there is only one key component that must be changed here, and that is a conversion factor between CAESAR's units (comoving kpc) and the pixel size of the frb you chose above. This conversion factor is calculated using the ``Generating_SZ_Data.determining_caesar_conversion()`` function. 


**Step Three: Data Analysis**

Now we get to the good part. The final inputs set what results we want from RAFIKI, with the following inputs:

.. code-block:: python

	sz_dat_path = 'noagn_'
	analysis = 'rpmpte'
	low_sm = 10**(11)
	width = 300
	beam_kernel = 1.68
	label = 'example'
	aperature = 50

- **sz_dat_path**: where is the .npy file generated above? Here we put the path to the file as well as the name, but without the suffix including the projection 
direction. This is to allow for the loop shown in the notebook to consider all three projections as one big data set.
 
- **analysis**: What analysis would you like to run? - RAFIKI is equipped to calculate three things from the tSZ data. The first is simply radial profiles of the 
Compton-y signal around each galaxy in your sample. If you would like this data, input rp. It can also generate radial profiles for the moments of symmetry for each galaxy. By default it will 
output moments 0, 1, and 2. If you would like this data, input mp. Finally, the tSZ values can be used to calculate the total thermal energy within a given radius around 
your galaxies. If you choose this option, you will be asked later on for the radius of the aperture you would like, in terms of pixels. If you would like this data, input 
te. You may choose any combination of these, and in the example above all three analyses will be conducted.

- **low_sm**: Minimum stellar mass for your sample - In order to speed up the analysis process, if you know you are only interested in large galaxies you may make a lower 
cut with this variable. Particularly at high redshifts, where only the most massive galaxies are visible, it is suggested to make a cut similar to the 1e11 shown above. This value must 
have units of solar masses.

- **width**: Pixel length of box you want to cut around each galaxy - Particularly relevant for the radial profile analyses, this input determines how far out from each 
galaxy you want to study. Note this is the length of the side of a stamp around the galaxy, not a radius of a circular aperature. Also note this is in terms of pixels, so make sure 
you have a good understanding of the relationship between the pixel size and physical or angular units.

- **beam_kernel**: Standard deviation of Gaussian kernel - A key element of RAFIKI is the ability to compare against observational data by convolving with a variety of 
beam sizes. In its current form, RAFIKI only convolves with a Gaussian beam, and this input determines the standard deviation of that beam (in terms of pixels).

- **label**: Output Label - String added to each output file. The outputs are described more in detail in OUTPUTS.

- **aperature**: Size of region for thermal energy measurements - Thermal energy is calculated in a circular aperature around each galaxy. Here, you set the radius of 
that aperature in pixels. 


Stacking and Visualizing Data
-----------------------------
The file stacking_example.ipynb shows how to take the csv files generated above and stack the data around a galaxy sample of your choice, creating plots. Once again, you 
can chose if you want radial profiles, moments, or thermal energy data. One key place for user input comes in specifying the galaxy sample. Where above, we set a lower 
limit on the stellar mass to reduce computational time, we now set full constraints. 

stacking_example.ipynb generates a sample with stellar masses larger than 1e11, ages above 1 Gyr, and SSFR less than 0.01, a way of selecting quiescent galaxies. 

This notebook both generates plots and csv files containing the plot data. Errors are calculated using a bootstrapping method. 




