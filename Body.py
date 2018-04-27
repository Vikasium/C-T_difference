# -*- coding: utf-8 -*-
"""
Created on Mon May 29 14:44:09 2017

Beans...... Star, Planet, Telescope_band
@author: vikash
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const

# Classes
class Star(object):
    """ Data from the Phoenix Stellar Spectrum files """
    def __init__(self, name):      
        
        self.name = name
        self.temperature = None
        self.radius = None
        self.radius_from_temperature = None
        self.mass = None
        self.RA = None
        self.dec = None
        self.metallicity = None 
        self.magnitude = []             # 4 component. V, J, H, K,
        self.B_minus_V = None
        self.distance = None
                                
        self.directory = None           # location of stellar spectra files
        self.filename = None
        self.variable_names = None
        self.wavelength_array = None
        self.flux_array = None
        self.BB_flux_array = None
        
        self.pol_linear = None
                    
    def read_property(self, file):
        with open(file,'r') as f:
             heading = f.readline().rstrip('\n').split(',')                              # Read the first line for headings
             
             for l in f.read().splitlines():
                 line = l.split(',')
                 nm = line[heading.index('STAR_NAME')].replace(" ","")
                 if nm == self.name:                                   # get the index of a parameter using heading
                     self.mass = float(line[heading.index('MSTAR')]) #* const.M_sun.value ..... coz this way the errors will propagate (MSing i = 1.19 +- 0.03)
                     self.radius = float(line[heading.index('RSTAR')]) #* const.R_sun.value
                     self.temperature = float(line[heading.index('TEFF')])
                     #self.log_G = line[heading.index('LOGG')]
                     self.metallicity = float(line[heading.index('FE')])
                     self.B_minus_V = float(line[heading.index('BMV')])
                     self.magnitude.append(float(line[heading.index('V')]))
                     self.magnitude.append(float(line[heading.index('J')]))
                     self.magnitude.append(float(line[heading.index('H')]))
                     self.magnitude.append(float(line[heading.index('KS')]))
                     self.distance = float(line[heading.index('DIST')])             # In parsecs
                     self.RA = line[heading.index('RA_STRING')]
                     self.dec = line[heading.index('DEC_STRING')]
             if self.mass == None:
                print('Star name- {} is not present (or Mass value is blank)'.format(self.name))
                exit(0)                 
    def set_file_name(self):
        temp = np.round(self.temperature/100.0)                                       # to find the corresponding file with round-off temperature 
        for file in os.listdir(self.directory):
            T = float(file.partition('.')[0].partition('e')[-1])                      # Temperature is in the short_name
            if T==temp:
                self.filename = file
        if self.filename == None:
            print('Did not find the stellar spectra file')
            exit(0)
    def read_spectra(self, DF = -8.0):
        infile = self.directory + '/' + self.filename
        data = np.genfromtxt(infile, delimiter=',', dtype='str')
        self.variable_names = data[0,:]
        float_data = data[1:,:].astype(float)
        self.wavelength_array = float_data[:,0]
        self.flux_array = pow(10,(float_data[:,1] + DF))
        self.BB_flux_array = pow(10,(float_data[:,2] + DF))                                        # 10**(F_lam + DF) to convert to Ergs/sec/cm**2/A
      
    def set_radius_from_temperature(self, p1):
        if(self.temperature != None):
            self.radius_from_temp = np.polyval(p1,self.temperature)                           # Solar radius
            print(self.temperature, self.radius_from_temp)
        
    def plot_spectra(self, x, y, ax1):
        #fig = plt.figure()    
        #ax1 = fig.add_subplot(1, 1, 1)
        ax1.set_xlabel(self.variable_names[0])
        ax1.set_ylabel(self.variable_names[1] + '(Ergs/sec/cm**2/A)')
        ax1.set_xlim(3300,11000)
        #ax1.set_ylim(ymin=0,ymax=5e6)
        ax1.plot(x, y)
        
#==============================================================================
#         ax2 = fig.add_subplot(2, 1, 2)
#         ax2.set_xlabel(self.variable_names[0])
#         ax2.set_ylabel(self.variable_names[2] + '(Ergs/sec/cm**2/A)')
#         ax2.set_xlim(3300,11000)
#         ax2.set_ylim(ymin=0,ymax=5e6)
#         ax2.plot(self.wavelength_array, self.BB_flux_array, 'r')
#==============================================================================
        
class Planet(object):
    """ Data from the Exotransmit Spectra files """
    def __init__(self,name):    
        self.name = name
        self.mass = None                # M_sin(i)
        self.temperature = None
        self.radius = None
        self.surface_G = None
        self.orbital_period = None 
        self.transit_depth = None
        self.transit_duration = None
        self.a_by_R = None
        
        self.directory = None                           # location of planet spectra_files
        self.filename = None
        self.wavelength_array = None
        self.transit_depth_array = None
        self.interpolated_depths = None
        
    def read_property(self, file):
        with open(file,'r') as f:
            heading = f.readline().rstrip('\n').split(',')                            # Read the first line for headings
             
            for l in f.read().splitlines():
                line = l.split(',')
                nm = line[heading.index('STAR_NAME')].replace(" ","")+'b'
                
                if nm == self.name:                                   # get the index of a parameter using heading
                    self.mass = float(line[heading.index('MSINI')]) #* const.M_jup.value
                    self.radius = float(line[heading.index('R')]) #* const.R_jup.value
                    self.surface_G = float(line[heading.index('GRAVITY')])
                    self.orbital_period = float(line[heading.index('PER')])
             
                    self.transit_depth = float(line[heading.index('DEPTH')])
                    self.transit_duration = float(line[heading.index('TRANS_DUR')])             # In parsecs
                    self.a_by_R = float(line[heading.index('A_BY_R')])
            if self.mass == None:
                print('Planet name is not present (or Mass value is blank). Property is not read')
                exit(0)       
    def read_spectra(self):
        infile = self.directory + '/' + self.filename
        #row1, row2 = np.loadtxt(infile, delimiter='\t', dtype='S', unpack=True)
        row1, row2 = np.loadtxt(infile, delimiter='        ', dtype='S', unpack=True)
        self.wavelength_array = row1[2:].astype(np.float) * 1e4                    # To convert into angstrom
        self.transit_depth_array = row2[2:].astype(np.float)                        # R_p/R*
        
    def set_interpolated_transit_depths(self, wavelengths):
        self.interpolated_depths = np.interp(wavelengths,self.wavelength_array,self.transit_depth_array)
        
    def set_planet_equilibrium_temperature(self, T_eff, R_star, f=0.25, A_bond=0.3):
        a = self.a_by_R * R_star
        self.temperature = T_eff * np.sqrt(R_star/a) * np.power(f * (1 - A_bond) , 0.25)
         
    def get_other_planet_spectra(self,this_planets_star_radius, planet_mass, planet_eq_temp, star_radius):
        # Assuming scale height ratio equal to 1
        mass_ratio = planet_mass/self.mass       # Both in units of M_jup
        temp_ratio = planet_eq_temp/self.temperature
        radius_ratio = this_planets_star_radius/star_radius      # both in units of R_sun
        
        interpolation_factor = (mass_ratio/temp_ratio) * np.power(radius_ratio,2)
        planet_spectra = interpolation_factor * self.interpolated_depths
        return planet_spectra
        
    def plot_spectra(self, x, y, ax):
        ax.set_title('Transit_Spectra of ' + self.name)
        ax.set_xlabel('Wavelength(A)')
        ax.set_ylabel('Transit Depth')
        ax.set_xlim(3300,11000)
        ax.semilogx(x, y)
        #ax.plot(wavelengths, self.interpolated_depths)
        #ax.plot(wavelengths, self.transit_depth_array)
        
class Telescope(object):
    """ Cheops and TESS response function from their respective throughputs """
    
    variable_names = ''
    def __init__(self,name):
    
        self.name = name
        self.wavelength_array = None
        self.response_array = None
        self.interpolated_responses = None
        self.variable_names = None
        
        self.directory = None
        self.filename = None
        self.radius_ratio_error = 0.0012766378/0.1104833233         #0.0024744097/0.0839589465 #0.0001733285 / 0.1615604542  #  Relative error. 0.1839086320
        
    def read_response(self):
        infile = self.directory + '/' + self.filename
        data = np.genfromtxt(infile, delimiter=' ', dtype='str')
        float_data = data[1:,:]
        
        self.variable_names = data[0,:]
        self.wavelength_array = float_data[:,0].astype(float) *10                    # To convert into angstrom
        self.response_array = float_data[:,1].astype(float)    
        
    def set_interpolated_responses(self, wavelengths):
        self.interpolated_responses = np.interp(wavelengths,self.wavelength_array,self.response_array, left=0, right=0)
        
    def plot_response(self,x, y, ax):
        #fig = plt.figure()
        #ax = fig.add_subplot(1,1,1)
        ax.set_title(self.name + ' ' + self.variable_names[1])
        ax.set_xlabel(self.variable_names[0].partition('(')[0] + '(A)')
        ax.set_ylabel('response')
        ax.set_xlim(3300,11000)
        ax.plot(x, y)        
        