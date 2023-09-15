'''
This script generates a projected map of the Compton-y signal for use in tSZ analysis. 
Inputs are set in lines 17-22 and explained more in detail in the documentation.
'''

import yt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import gridspec
import numpy as np
import caesar
import csv
import pandas as pd



'''
You shouldn't need to adjust below this line
-------------------------------------------------------------
'''

def generating_sz_data(filename, projection_direction, output_name, frb, redshift):
    '''
    Cuts ISM and current wind particles, makes a new field of gas pressure, projects this field in the specified direction, and saves a file containing the projected map of the tSZ signal

    :param filename: path to snapshot
    :type filename: str
    :param projection_direction: 'x', 'y', or 'z', the direction you want to project your data in
    :type projection_direction: str
    :param output_name: name of the .npy file that will store the projected SZ data
    :type output_name: str
    :param frb: number of pixels in your fixed resolution buffer. Suggested to correspond to twice the resolution of your observational comp
    :type frb: int
    :param redshift: redshift of the snapshot
    :type redshift: int
    :return: None. Saves a file with the name given in the inputs containing the projected SZ data
    
    '''



    obj = yt.load(filename)
    #Creating derived field of gas pressure
    def _pressure(field, data):
        return (
            data["PartType0", "density"]
            * data["PartType0", "Temperature"]
        )

    obj.add_field(
        name=("PartType0", "pressure"),
        function=_pressure,
        sampling_type="local",
        units="K*code_mass/code_length**3",
    )

    #Cutting ISM and wind particles
    def cuts(pfilter, data):
        filter = np.logical_and(data[(pfilter.filtered_type, "H_nuclei_density")]  < 0.1, 
                                data[(pfilter.filtered_type, "DelayTime")]  <= 0 )
        return filter
    yt.add_particle_filter(
        "szcuts", function=cuts, filtered_type="PartType0",requires=("H_nuclei_density","DelayTime"))
    obj.add_particle_filter("szcuts")

    #Generating projection plot, can change direction of projection in second input
    prj = yt.ProjectionPlot(obj, projection_direction, ('szcuts', 'pressure'))

    #Generating an frb to set pixel size of your map. Resolution currently set to be 1/2 the smallest beam width
    prj.set_buff_size((frb,frb)) 
    data = prj.frb[('szcuts', 'pressure')]

    sz_dat3a = np.array(data)*1.11*10**(-32)*(1+redshift)**3 
    sz_dat = np.array(sz_dat3a)
    #Saves 2D array of projected SZ-y data as a .npy file
    with open(output_name, 'wb') as f:
        np.save(f, np.array(sz_dat))