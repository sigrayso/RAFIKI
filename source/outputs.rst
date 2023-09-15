RAFIKI Outputs
**************

At various stages in running RAFIKI, output files will be generated. Nearly all will be csv files and their data is explained here ([label] indicates a user input label 
chosen when running. For more details, see :doc:`running`).

**[label]_szy.npy** - Output of Generating_SZ_Data.py. Image data of the projected map of the Compton-y parameter. 

**galaxy_sample_[label].csv** - Contains properties of your galaxy sample. Data is output in five columns: 
Stellar Mass, Halo Mass, R200, Age, SFR. Retaining the order of this csv is very important, as it is the same as 
the files below. 

**radial_[label].csv** - Radial profile of Compton-y signal around each galaxy

**m_0_[label].csv** - Radial profile of the zeroth moment around each galaxy

**m_1_[label].csv** - Radial profile of the first moment around each galaxy

**m_2_[label].csv** - Radial profile of the second moment around each galaxy

**thermal_energy_[label].csv** - Thermal energy within the radius specified when running around each galaxy

You should never need to directly deal with any of these outputs, as they will go directly into the stacking and 
visualization process. The outputs of that process will be more useful and are as follows:

**stacked_radial_profiles_[label].csv** - Radial profile of stacked galaxies. Columns are radius (arcmin), Compton-y signal, and y-error

**stacked_moment_profiles_[label.csv** - Radial profile of the moments of stacked galaxies. Columns are radius (arcmin), ratio of M1/M0, the y-error of that, ratio of 
M2/M0, and y-error of that

**thermal_energy_plots_[label].csv** - Thermal energy against mass plot data. Columns are stellar mass bins, thermal energy in each bin, error, halo mass bins, thermal 
energy in each bin, error. 


