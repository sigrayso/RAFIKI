import yt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import gridspec
import numpy as np
import caesar
import pandas as pd
import csv
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from astropy.convolution import convolve, convolve_fft
import scipy.ndimage
from random import randint
import math
from scipy.ndimage.interpolation import geometric_transform
#from Galaxy_data import *
from scipy.stats import bootstrap

def sorting(attribute,minimum):
    '''
    Identifies galaxies by a specified attribute that are above the minimum value for that attribute. Used to cut low-mass galaxies to speed up the analysis process
        
    :param attribute: list of galaxy stellar masses (can be altered to other galaxy properties if desired)
    :type attribute: list[float]
    :param mininimum: minimum value you want to allow in the sample
    :type minimum: float
    :return: gals_index-index of the galaxies in Caesar catalog that are above the set minimum
    :rtype: list[int]

    '''
    gals_index = []
    for i in range(len(attribute)):
        if attribute[i]>minimum:
            gals_index.append(i)

    #GPT: gals_index = [i for i, mass in enumerate(attribute) if mass > minimum]

    return gals_index   

def cropping_sz_x(sz_dat, xs, ys, zs, sample, width, frb):
    ''' 
    Generates a nested 2D array of SZ data in the regions around our sampled galaxies using a data box projected in the x-direction
        
    :param sz_dat: 2D array of projected Compton-y map
    :type sz_dat: list[float]
    :param xs: x-coordinates of locations of galaxies in your box, adjusted to match the pixel frb
    :type xs: list[float]
    :param ys: y-coordinates of locations of galaxies in your box, adjusted to match the pixel frb
    :type ys: list[float]
    :param zs: z-coordinates of locations of galaxies in your box, adjuste to match the pixel frb
    :type zs: list[float]
    :param sample: list of index locations of the galaxies we want to sample (output of ``sorting``)
    :type sample: list[int]
    :param width: the pixel size we want to crop around the galaxies
    :type width: int
    :param frb: pixel size of the SZ data used
    :type frb: int
    :return: lists-Nested array of SZ data for the regions around our sample gaalxies
    :rtype: list[float]
        
    '''

    lists = []
    q = 0
    for i in sample:
        indices_z = range(int(zs[i])-int(width/2),int(zs[i])+int(width/2))
        indices_y = range(int(ys[i])-int(width/2),int(ys[i])+int(width/2))
        #print(indices_z)
        #print(indices_y)
        if width/2<zs[i] < frb-(width/2) and width/2<ys[i] < frb-(width/2): 
            lists.append([])
            lists[q].append(sz_dat[int(zs[i]-width/2):int(zs[i]+width/2),int(ys[i]-width/2):int(ys[i]+width/2)])
            q+=1
        else: #Wrapping data around periodic box if needed
            lists.append([])
            z = sz_dat.take(indices_z, axis = 0, mode='wrap').take(indices_y, axis = 1, mode='wrap')
            lists[q].append(z)
            q+=1
    return lists

def cropping_sz_y(sz_dat, xs, ys, zs,sample, width, frb):
    ''' 
    Generates a nested 2D array of SZ data in the regions around our sampled galaxies using a data box projected in the y-direction
        
    :param sz_dat: 2D array of projected Compton-y map
    :type sz_dat: list[float]
    :param xs: x-coordinates of locations of galaxies in your box, adjusted to match the pixel frb
    :type xs: list[float]
    :param ys: y-coordinates of locations of galaxies in your box, adjusted to match the pixel frb
    :type ys: list[float]
    :param zs: z-coordinates of locations of galaxies in your box, adjuste to match the pixel frb
    :type zs: list[float]
    :param sample: list of index locations of the galaxies we want to sample (output of ``sorting``)
    :type sample: list[int]
    :param width: the pixel size we want to crop around the galaxies
    :type width: int
    :param frb: pixel size of the SZ data used
    :type frb: int
    :return: lists-Nested array of SZ data for the regions around our sample gaalxies
    :rtype: list[float]
        
    '''
    lists = []
    q = 0
    for i in sample:
        indices_z = range(int(zs[i])-int(width/2),int(zs[i])+int(width/2))
        indices_x = range(int(xs[i])-int(width/2),int(xs[i])+int(width/2))

        if width/2<zs[i] < frb-(width/2) and width/2<xs[i] < frb-(width/2): 

            lists.append([])
            lists[q].append(sz_dat[int(xs[i]-width/2):int(xs[i]+width/2),int(zs[i]-width/2):int(zs[i]+width/2)])
            q+=1
        else: #Wrapping data around periodic box if needed
            lists.append([])
            z = sz_dat.take(indices_x, axis = 0, mode='wrap').take(indices_z, axis = 1, mode='wrap')
            lists[q].append(z)
            q+=1
    return lists

def cropping_sz_z(sz_dat, xs, ys, zs, sample, width, frb):
    ''' 
    Generates a nested 2D array of SZ data in the regions around our sampled galaxies using a data box projected in the z-direction
        
    :param sz_dat: 2D array of projected Compton-y map
    :type sz_dat: list[float]
    :param xs: x-coordinates of locations of galaxies in your box, adjusted to match the pixel frb
    :type xs: list[float]
    :param ys: y-coordinates of locations of galaxies in your box, adjusted to match the pixel frb
    :type ys: list[float]
    :param zs: z-coordinates of locations of galaxies in your box, adjuste to match the pixel frb
    :type zs: list[float]
    :param sample: list of index locations of the galaxies we want to sample (output of ``sorting``)
    :type sample: list[int]
    :param width: the pixel size we want to crop around the galaxies
    :type width: int
    :param frb: pixel size of the SZ data used
    :type frb: int
    :return: lists-Nested array of SZ data for the regions around our sample gaalxies
    :rtype: list[float]
        
    '''
    lists = []
    q = 0
    for i in sample:
        indices_y = range(int(ys[i])-int(width/2),int(ys[i])+int(width/2))
        indices_x = range(int(xs[i])-int(width/2),int(xs[i])+int(width/2))

        if width/2<ys[i] < frb-(width/2) and width/2<xs[i] < frb-(width/2): 

            lists.append([])
            lists[q].append(sz_dat[int(ys[i]-width/2):int(ys[i]+width/2),int(xs[i]-width/2):int(xs[i]+width/2)])
            q+=1
        else: #Wrapping data around periodic box if needed
            lists.append([])
            z = sz_dat.take(indices_y, axis = 0, mode='wrap').take(indices_x, axis = 1, mode='wrap')
            lists[q].append(z)
            q+=1
    return lists

def azimuthalAverage(image, center=None):
    """
    Calculates the azimuthally averaged radial profile. Taken from with some alterations https://github.com/mkolopanis/python/blob/master/radialProfile.py
        
    :param image: 2D array of projected Compton-y map around our galaxy of interest
    :type image: list[float]
    :param center: The [x,y] pixel coordinates used as the center. Default is None, which useds the center of the image
    :type center: list[float] or None
    :return: radial_prof-Radial profile of the azimuthally averaged signal
    :rtype: list[float]
    
    """
    # Calculate the indices from the image

    y, x = np.indices(np.shape(image))

    if not center:
        center = np.array([(x.max()-x.min())/2.0+1, (x.max()-x.min())/2.0+1]) #+1 added because rounding down below
    
    
    r = np.hypot(x - center[0], y - center[1]) #np.hypot gives the hypotenuse of a triangle with the given legs
    
    # Get sorted radii
    ind = np.argsort(r.flat)  

    r_sorted = r.flat[ind]  #Sorted list of the radii of the pixels, converted to integers below
    i_sorted = image.flat[ind]  #Sorting the image pixels by the radii
     
    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]# Assumes all radii represented
    rind3 = np.where(deltar)[0] # location of changed radius
    rind = np.insert(rind3, 0, 0)
    nr = rind[1:] - rind[:-1] # number of radius bin
    
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]
    
    radial_prof = tbin / nr
    
    #Below is my addition for calculating error, breaks up pixels into separate sub-arrays by radius range, averages. 
    split = [i_sorted[rind[r]:rind[r+1]] for r in range(len(rind)-1)]
 
    errors = [np.std(split[i])/np.sqrt(len(split[i]))for i in range(len(split))]
    return radial_prof




def make_radial_profiles(stamps, kernel, label):
    '''
    Convolves data and generates radial profiles for all galaxies in your sample
        
    :param stamps: A nested array of the SZ data around galaxies. Output of ``cropping_sz`` function
    :type stamps: list[float]
    :param kernel: The standard deviation of the Gaussian kernel in units of pixels
    :type kernel: float
    :param label: Label you want in the file name containing radial profiles
    :type label: str
    :return: Makes a csv file entitled radial_label.csv containing radial profiles around each galaxy in your sample
    '''

    gauss_kernel = Gaussian2DKernel(kernel) 
    a = []
    for i in stamps: 
        dd =convolve_fft(i[0], gauss_kernel)     
        aa=azimuthalAverage(dd, center = None)     
        a.append(aa)


    with open('radial_%s.csv'%label,'a', newline = '') as f: 
        w = csv.writer(f)
        for z in a:
            c = w.writerow(z)

    return 



"""FUNCTIONS FOR FINDING MOMENTS """

def topolar(img, order=1):
    """
    Transforms an image into polar coordinates
        
    :param img: 2D data of the image wanting to transform
    :type img: list[float]
    :param order: The spline interpolation order, default 1
    :type order: int
    :return: Polar-Nested array of image by polar coordinates, (rads, angs)- Values of the radii and angles corresponding to the data in polar
    :rtype: list[float]
    """
    # max_radius is the length of the diagonal 
    # from a corner to the mid-point of img.
    max_radius = 0.5*np.linalg.norm( img.shape )

    def transform(coords):
        # Put coord[1] in the interval, [-pi, pi]
        theta = 2*np.pi*coords[1] / (img.shape[1] - 1.)

        # Then map it to the interval [0, max_radius].
        #radius = float(img.shape[0]-coords[0]) / img.shape[0] * max_radius
        radius = max_radius * coords[0] / img.shape[0]

        i = 0.5*img.shape[0] - radius*np.sin(theta)
        j = radius*np.cos(theta) + 0.5*img.shape[1]
        return i,j

    polar = geometric_transform(img, transform, order=order)

    rads = max_radius * np.linspace(0,1,img.shape[0])
    angs = np.linspace(0, 2*np.pi, img.shape[1])

    return polar, (rads, angs)

def make_moment_profiles(stamps, kernel, label):
    '''
     Generates radial profiles of transformed maps for m=0, 1, and 2
        
    :param stamps: A nested array of the SZ data around galaxies. Output of ``cropping_sz`` function
    :type stamps: list[float]
    :param kernel: The standard deviation of the Gaussian kernel in units of pixels
    :type kernel: float
    :param label: Label you want in the file names containing radial profiles
    :type label: str
    :return: Saves three csv files containing the moments radial profiles for m=0, 1, and 2
    '''
    data0 = []
    data1 = []
    data2 = []
    n=0
    gauss_kernel = Gaussian2DKernel(kernel) 
    for i in stamps: 
        d =convolve_fft(i[0], gauss_kernel)  
      
        pol, (rads,angs) = topolar(d)
        reals0 = []
        imaginaries0 = []
        reals1 = []
        imaginaries1 = []
        reals2 = []
        imaginaries2 = []
        for k in list(pol):
            
            real0 = [k[j]*math.cos(0*angs[j]) for j in range(len(k))]
            imag0 = [k[j]*math.sin(0*angs[j]) for j in range(len(k))]
            reals0.append(np.sum(real0)/(len(real0)))
            imaginaries0.append(np.sum(imag0)/(len(imag0)))
            real1 = [k[j]*math.cos(1*angs[j]) for j in range(len(k))]
            imag1 = [k[j]*math.sin(1*angs[j]) for j in range(len(k))]
            reals1.append(np.sum(real1)/(len(real1)))
            imaginaries1.append(np.sum(imag1)/(len(imag1)))
            real2 = [k[j]*math.cos(2*angs[j]) for j in range(len(k))]
            imag2 = [k[j]*math.sin(2*angs[j]) for j in range(len(k))]
            reals2.append(np.sum(real2)/(len(real2)))
            imaginaries2.append(np.sum(imag2)/(len(imag2)))
        
        
        amplitude0 = [np.sqrt(reals0[h]**2+imaginaries0[h]**2) for h in range(len(reals0))]
        amplitude1 = [np.sqrt(reals1[h]**2+imaginaries1[h]**2) for h in range(len(reals1))]
        amplitude2 = [np.sqrt(reals2[h]**2+imaginaries2[h]**2) for h in range(len(reals2))]
        data0.append(amplitude0)
        data1.append(amplitude1)
        data2.append(amplitude2)
        #data.append(profile)
        n+=1
    '''with open('radii.csv','a', newline = '') as f:
        w = csv.writer(f)
            
        c = w.writerow(rads)'''
    with open('m_0_%s.csv'% label,'a', newline = '') as f: 
        w = csv.writer(f)
        for z in data0:
            c = w.writerow(z)
    with open('m_1_%s.csv'% label,'a', newline = '') as f: 
        w = csv.writer(f)
        for z in data1:
            c = w.writerow(z)
    with open('m_2_%s.csv'% label,'a', newline = '') as f: 
        w = csv.writer(f)
        for z in data2:
            c = w.writerow(z)


"""THERMAL ENERGY"""

def thermal_energy(stamps, kernel, label, aperature, comov, z, theta):
    """
    Calculates the thermal energy in a given aperature from maps of SZ-y data
        
    :param stamps: A nested array of the SZ data around galaxies. Output of ``cropping_sz`` function
    :type stamps: list[float]
    :param kernel: The standard deviation of the Gaussian kernel in units of pixels
    :type kernel: float
    :param label: Label you want in the file name containing thermal energy data
    :type label: str
    :param aperature: pixel size of the radius of the aperature within which you want to caculate thermal energy
    :type aperature: float
    :param comov: comoving distance in Gpc
    :type comov: float
    :param z: Redshift of snapshot
    :type z: float
    :param theta: pixel size in arcseconds
    :type theta: float
    :return: Saves a csv file containing the thermal energy around each of the galaxies in your sample

    """

    gauss_kernel = Gaussian2DKernel(kernel) 
    energies = []
    for i in stamps:
        image= convolve_fft(i[0], gauss_kernel)
        y, x = np.indices(image.shape)
        
        center = np.array([(x.max()-x.min())/2.0+1, (x.max()-x.min())/2.0+1]) #+1 added because rounding down below
        r = np.hypot(x - center[0], y - center[1]) #np.hypot gives the hypotenuse of a triangle with the given legs
        # Get sorted radii
        ind = np.argsort(r.flat)  
        r_sorted = r.flat[ind]  #Sorted list of the radii of the pixels, converted to integers below
        i_sorted = image.flat[ind]  #Sorting the image pixels by the radii 
        # Get the integer part of the radii (bin size = 1)
        r_int = r_sorted.astype(int)

        new_r_int=[]
        for i in r_int:
            if i<= aperature:
                new_r_int.append(i)

        l = len(new_r_int)
        i_sorted_new = i_sorted[0:l]
        total = np.sum(i_sorted_new)
        energies.append(total)
    
    therm = []
    for i in energies:
        therm.append(2.9 * (comov/(1+z))**2 * (i*(theta/60)**2)/(10**(-6)))


    with open('thermal_energy_%s.csv'% label,'a', newline = '') as f: 
        w = csv.writer(f)
        c = w.writerow(therm)


#Bootstrapping tools

def gen_random_indices(index_set, gen_size):
    """
    Generates a list of indicies by random sampling with replacement
        
    :param index_set: List of values to sample from
    :type index_set: list[float]
    :param gen_size: Length of final resampled list you want
    :type gen_size: int
    :return: a list of length gen_size randomly chosen from index_set
    """
    return np.random.choice(index_set, size=gen_size, replace=True)

def single_catalog_bootstrap(data, boot_size, loop_size):
    """
    Calculates the means of a list of catalogs, useful when determining things like correlation matrices

    :param data: List of pandas tables
    :type data: list[float]
    :param boot_size: Sample size
    :type boot_size: int
    :param loop_size: How many times to do the bootstrapping
    :type boot_size: int
    """
    if type(data)!=list:
        print("Data must be provided as a list... exiting...")
        return None
    dlen = len(data)
    dlen2 = len(data[0])    
    print(dlen, len(data[0]))
    
    result = []
    for l in range(loop_size):
        indices = gen_random_indices(np.arange(dlen2), boot_size)
        result.append([np.mean(np.take(data[d], indices)) for d in range(dlen)])
    result = np.array(result)
    return result


def make_thermal_energy_plot(data, thresholds_stellar, thresholds_halo):
    '''
    Generates data and error on thermal energy values for stacked galaxies

    :param data: Nested list containing stellar masses of galaxies, halo masses of galaxies, and thermal energy around each galaxy
    :type data: list[float]
    :param thresholds_stellar: Bin limits for stacking by stellar mass
    :type thresholds_stellar: list[float]
    :param thresholds_halo: Bin limits for stacking by halo mass
    :type thresholds_halo: list[float]
    :return: y-Thermal energy of stacked galaxies in each stellar mass bin, err-Error on thermal energy of stacked galaxies in each stellar mass bin, y2$
    :rtype: list[float]

    '''
    bins = [[] for _ in range(len(thresholds_stellar) - 1)]
    therms = [[] for _ in range(len(thresholds_stellar) - 1)]


    for i, s in enumerate(data[0]):
        for j, threshold in enumerate(thresholds_stellar[:-1]):
          
            if threshold < s < thresholds_stellar[j + 1]:
                bins[j].append(i)
                therms[j].append(data[2][i])
                break

    bins2 = [[] for _ in range(len(thresholds_halo) - 1)]
    therms2 = [[] for _ in range(len(thresholds_halo) - 1)]

    for i, s in enumerate(data[1]):
        for j, threshold in enumerate(thresholds_halo[:-1]):
            if threshold < s < thresholds_halo[j + 1]:
                bins2[j].append(i)
                therms2[j].append(data[2][i])
                break
    # Combine results for therm and mass
    comb = therms
    comb_mass = [data[0][b] for b in bins]
    comb2 = therms2
    comb_mass2 = [data[1][b] for b in bins2]


    y = []
    err = []
    for j in range(len(thresholds_stellar)-1):
        i = comb[j]
        if len(i) != 0:
            y.append(np.mean(i))
            ii = (i,)
            bootstrap_ci = bootstrap(ii, np.mean, confidence_level=0.95,
                         random_state=1, method='percentile')
            err.append(bootstrap_ci.standard_error)

    y2 = []
    err2 = []
    for j in range(len(thresholds_halo)-1):
        i = comb2[j]
        if len(i) != 0:
            y2.append(np.mean(i))
            ii = (i,)
            bootstrap_ci = bootstrap(ii, np.mean, confidence_level=0.95,
                         random_state=1, method='percentile')
            err2.append(bootstrap_ci.standard_error)

    return(y,err,y2,err2)
