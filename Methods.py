# -*- coding: utf-8 -*-
"""
Created on Mon May 29 14:37:04 2017

Methods called in the "main2" and not related to the "Body" objects.
 
Methods - Formulae - Utility
@author: vikash
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
import importlib.util
spec = importlib.util.spec_from_file_location("Body.py", "/vik/WORK/TESS-CHEOPS/Code/Body.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

file = '/vikash/WORK/TESS-CHEOPS/DVS_stars_planets.csv'
ylim = 1e5
# Functions
def read_spectras(star,planet,directory):                        #Band_1 = Cheops. Band_2 = TESS
    #star.read_property(file)
    #planet.read_property(file)  
    #planet.set_planet_equilibrium_temperature(star.temperature,star.radius)
    
    star.directory = directory + '/Stellar_Spectra/interpolated'
    #planet.directory = directory + '/Planet_Spectra/D_Sing' + '/' + planet.name
    planet.directory = directory + '/Planet_Spectra/GJ1214' 

    star.set_file_name()                                    # Using Temperature
    #planet.filename = 'Spectra'
    planet.filename = 'spectrum_GJ1214b.txt'

    # Reading Spectra and Responses                       
    star.read_spectra()                                           # ; star.plot_spectra()               # Plot here for its entire range 
    planet.read_spectra()
    
    #planet.radius = np.mean(planet.transit_depth_array) * star.radius
            
def adjust_planetary_spectra(star, planet):
    planet.transit_depth_array * (star.radius / star.radius_from_temp)   
    
def get_composite_stellar_flux(activity, star, star_spot):   
    temp_diff = activity[0]
    filling_factor = activity[1]
    star_spot.temperature = star.temperature - temp_diff
    #star_spot.directory = directory + '/Stellar_Spectra/interpolated'
    star_spot.set_file_name()
    star_spot.read_spectra()
    
    A = (100 - filling_factor) / 100
    B = filling_factor / 100
    composite_stellar_flux = A * star.flux_array + B * star_spot.flux_array
    #star.flux_array = composite_stellar_flux
    
#==============================================================================
#     fig = plt.figure(2)
#     ax1 = fig.add_subplot(4,1,1)
#     ax1.set_xlabel('wavelength')
#     ax1.set_ylabel('Flux')
#     ax1.set_xlim(3000,11500)
#     ax1.plot(star.wavelength_array, star.flux_array)
#     
#     difference = star.flux_array - star_spot.flux_array
#     ax2 = fig.add_subplot(4,1,2)
#     ax2.set_xlabel('wavelength')
#     ax2.set_ylabel('Flux')
#     ax2.set_xlim(3000,11500)
#     ax2.plot(star.wavelength_array, difference)  
#     
#     ax3 = fig.add_subplot(4,1,3)
#     ax3.set_xlabel('wavelength')
#     ax3.set_ylabel('Flux')
#     ax3.set_xlim(3000,11500)
#     ax3.plot(star.wavelength_array, composite_stellar_flux)
#     
#     ax4 = fig.add_subplot(4,1,4)
#     ax4.set_xlabel('wavelength')
#     ax4.set_ylabel('Flux')
#     ax4.set_xlim(3000,11500)
#     ax4.plot(star.wavelength_array, star.flux_array)
#     ax4.plot(star.wavelength_array, star_spot.flux_array)
#     ax4.plot(star.wavelength_array, composite_stellar_flux)
#==============================================================================
    
    return composite_stellar_flux
          
          
def get_activity_contribution(activity, star, star_spot, planet, band_1, band_2):
    
    temp_diff = activity[0]
    filling_factor = activity[1]
    star_spot.temperature = star.temperature - temp_diff
    #star_spot.directory = directory + '/Stellar_Spectra/interpolated'
    star_spot.set_file_name()
    star_spot.read_spectra()    
    
    surface_flux = star.flux_array
    spot_flux = star_spot.flux_array
    wavelength_range = star.wavelength_array
    
# Photons
    #surface_photons = flux_to_photons(surface_flux, wavelength_range)
    #spot_photons = flux_to_photons(spot_flux, wavelength_range)
    
    transit_depth_times_stellar_flux = surface_flux * np.square(planet.interpolated_depths)       # Fraction of blocked lights * stellar flux = Blocked/Scattered stellar Flux
    active_part_times_stellar_flux = surface_flux * (1 - (1 - spot_flux / surface_flux) * filling_factor)    
    
    transit_photons = flux_to_photons(transit_depth_times_stellar_flux, wavelength_range)
    active_photons = flux_to_photons(active_part_times_stellar_flux, wavelength_range)
    
# Responses
    band_1_active_photons = active_photons * band_1.interpolated_responses                      # Stellar flux in the pass bands
    band_2_active_photons = active_photons * band_2.interpolated_responses
    
    band_1_blocked_photons = transit_photons * band_1.interpolated_responses                    
    band_2_blocked_photons = transit_photons * band_2.interpolated_responses
    
    band_1_depth = band_1_blocked_photons.sum()/band_1_active_photons.sum()
    band_2_depth = band_2_blocked_photons.sum()/band_2_active_photons.sum()

    avg_depth = (band_1_depth + band_2_depth)/2.0
    diff = band_1_depth - band_2_depth
    percent_diff = (diff / avg_depth) * 1e2
    
    error_in_percent_difference = calculate_error(band_1)

    output = [band_1_depth, band_2_depth, avg_depth, diff, percent_diff, error_in_percent_difference]
    return output
    
def read_band_responses(band_1, band_2, directory):
    # Reading Telescope band Responses
    band_1.directory = directory + '/Telescope_Response'
    band_2.directory = directory + '/Telescope_Response'
    
    band_1.filename = 'CHEOPSresponse.dat' #'V_band.dat'#
    band_2.filename = 'TessResponse.dat' #'I_band.dat'#'
    
    band_1.read_response()
    band_2.read_response()
           
def setting_interpolated_ranges(star,planet,band_1,band_2):
    # Interpolation
    wavelength_range = star.wavelength_array                                            # range of wavelength data points for 3300-11027.16 A (From 2000A)
    stellar_flux = star.flux_array                                                      # Theoretical flux in Cheops window   

    star.wavelength_array = wavelength_range
    star.flux_array = stellar_flux
    planet.set_interpolated_transit_depths(wavelength_range)
    band_1.set_interpolated_responses(wavelength_range)
    band_2.set_interpolated_responses(wavelength_range) 
    #planet.interpolated_depths = transit_depth(wavelength_range,0.0241221)
    #planet.interpolated_depths = planet.interpolated_depths*transit_depth(wavelength_range,0.0241221)
    
    #print('Successfully read Spectra and Responses and the interpolation range is {} - {}'.format(wavelength_range[0],wavelength_range[-1]))
    
def flux_to_photons(flux,wavelength): 
    h = const.h.value * 1e+7
    c = const.c.value * 1e+2
    energy = (h*c)/wavelength 
    noofphotons = flux/energy
    
    return noofphotons

def calculate_percentage_difference(star,planet,band_1,band_2):
    output = None
    stellar_flux = star.flux_array
    wavelength_range = star.wavelength_array
    
# Photons
    stellar_photons = flux_to_photons(stellar_flux, wavelength_range)
    transit_depth_times_stellar_flux = stellar_flux * np.square(planet.interpolated_depths)       # Fraction of blocked lights * stellar flux = Blocked/Scattered stellar Flux
    transit_photons = flux_to_photons(transit_depth_times_stellar_flux, wavelength_range)

# Responses
    band_1_product = stellar_photons * band_1.interpolated_responses                      # Stellar flux in the pass bands
    band_2_product = stellar_photons * band_2.interpolated_responses
    
    band_1_blocked_photons = transit_photons * band_1.interpolated_responses                    
    band_2_blocked_photons = transit_photons * band_2.interpolated_responses

    band_1_depth = band_1_blocked_photons.sum()/band_1_product.sum()
    #band_1_depth = ((band_1_blocked_photons+1)/(band_1_photon_number+1)).sum()             # +1 to avoid 0/0 scenarios
    band_2_depth = band_2_blocked_photons.sum()/band_2_product.sum()
    #band_2_depth = ((band_2_blocked_photons+1)/(band_2_photon_number+1)).sum()
    
    avg_depth = (band_1_depth + band_2_depth)/2.0
    diff = band_1_depth - band_2_depth
    percent_diff = (diff / avg_depth) * 1e2
    
    error_in_percent_difference = calculate_error(band_1)

    output = [band_1_depth, band_2_depth, avg_depth, diff, percent_diff, error_in_percent_difference]
    return output
    
def calculate_error(band_1):              # sigma_band = uncertainity in measuring the transit depth.
#==============================================================================
#     sigma_of_percent_difference = None    
#     sigma_diff = np.sqrt(np.square(sigma_band_1) + np.square(sigma_band_2))
#     #sigma_avg = 0.5 * np.sqrt(np.square(sigma_band_1) + np.square(sigma_band_2))
#     sigma_avg = sigma_band_1
#     
#     sigma_of_percent_difference = 100 * (diff_in_depth/avg_depth) * np.sqrt(np.square(sigma_diff/diff_in_depth) + np.square(sigma_avg/avg_depth))
#==============================================================================
    sigma_of_percent_difference = 4 * band_1.radius_ratio_error * 100
    
    return sigma_of_percent_difference
      
def percentage_difference_with_activity(activity, star, star_spot, planet,band_1,band_2):

    hot_stellar_flux = star.flux_array
    wavelength_range = star.wavelength_array  

    composite_flux = get_composite_stellar_flux(activity, star, star_spot)
    composite_photons =  flux_to_photons(composite_flux,wavelength_range)
       
    transit_depth_times_stellar_flux = hot_stellar_flux * np.square(planet.interpolated_depths)    # Fraction of blocked lights * stellar flux = Blocked/Scattered stellar Flux
    transit_photons = flux_to_photons(transit_depth_times_stellar_flux,wavelength_range)
    
    band_1_stellar_photons = composite_photons * band_1.interpolated_responses                      # Stellar flux in the pass bands
    band_2_stellar_photons = composite_photons * band_2.interpolated_responses
          
    band_1_blocked_photons = transit_photons * band_1.interpolated_responses                    
    band_2_blocked_photons= transit_photons * band_2.interpolated_responses
    
    band_1_depth = band_1_blocked_photons.sum()/band_1_stellar_photons.sum()
    band_2_depth = band_2_blocked_photons.sum()/band_2_stellar_photons.sum()
    
    avg_depth = (band_1_depth + band_2_depth)/2.0
    diff = band_1_depth - band_2_depth
    percent_diff = (diff / avg_depth) * 1e2

#==============================================================================
#     fig = plt.figure(1)
#     ax2 = fig.add_subplot(2,1,1)
#     ax2.set_xlabel('wavelength')
#     ax2.set_ylabel('Cheops_product')
#     ax2.set_xlim(3000,11500)
#     ax2.plot(star.wavelength_array, band_1_stellar_photons, 'b')
#     ax2.plot(star.wavelength_array, band_2_stellar_photons, 'r')
#==============================================================================
    
#==============================================================================
#     ax3 = fig.add_subplot(2,1,1)
#     ax3.set_xlabel('wavelength')
#     ax3.set_ylabel('Tess_product')
#     ax3.set_xlim(3000,11500)
#     ax3.plot(star.wavelength_array, band_2_product)
#     
#==============================================================================
#==============================================================================
#     ax1 = fig.add_subplot(2,1,2)
#     ax1.set_xlabel('wavelength')
#     ax1.set_ylabel('Cheops_blocked_product')
#     ax1.set_xlim(3000,11500)
#     ax1.plot(star.wavelength_array, band_1_blocked_photons, 'b')
#     ax1.plot(star.wavelength_array, band_2_blocked_photons, 'r')
#==============================================================================
    
#==============================================================================
#     ax4 = fig.add_subplot(2,1,2)
#     ax4.set_xlabel('wavelength')
#     ax4.set_ylabel('Cheops_blocked_product')
#     ax4.set_xlim(3000,11500)
#==============================================================================
    
    #sigma_band_1 = np.sqrt(2) * band_1.radius_ratio_error * band_1_depth
    #sigma_band_2 = np.sqrt(2) * band_2.radius_ratio_error * band_2_depth
    error_in_percent_defference = calculate_error(band_1)

    output = [band_1_depth, band_2_depth, avg_depth, diff, percent_diff, error_in_percent_defference]
    return output    
    
    
def plot_spectras(fig, star, planet, band_1, band_2):
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    star.plot_spectra(star.wavelength_array, star.flux_array, ax1)
    planet.plot_spectra(star.wavelength_array, planet.interpolated_depths, ax2)
    band_1.plot_response(star.wavelength_array, band_1.interpolated_responses, ax3)
    band_2.plot_response(star.wavelength_array, band_2.interpolated_responses, ax4)
    
def transit_depth(lamda,avg_depth):
    radius_ratio = 0.0821 - 0.0075*np.log(lamda/3500.0)                        # When lamda is in angstrom
    depth = np.square(radius_ratio)    
    return depth
   
def calculate_alpha_H(star, planet):
    index_zero = np.argwhere(np.floor(star.wavelength_array) == 3000.0)[0][0]
    index_last = np.argwhere(np.floor(star.wavelength_array) == 4000.0)[0][0]    # Before the sodium lines
    
    wavelength = star.wavelength_array[index_zero : index_last] 
    #flux = star.BB_flux_array[index_zero : index_last]
    transit_depth = planet.interpolated_depths [index_zero : index_last]   
    
    #radius_array = np.sqrt((np.square(star.radius) / flux) * transit_depth )
    plt.semilogx(np.log(wavelength), transit_depth)
    alpha_H = - star.radius * np.polyfit(np.log(wavelength), transit_depth, 1)[0]   #to get the slope * star.radius
    return alpha_H
    
def create_ExoTransmit_input(star, planet):
    cwd =  os.getcwd()
    output_dir = '/vikash/TransmissionSpectroscopy/ExoTransmit/Exo_Transmit-master'
    infile = cwd + '/' + 'userInputTemplate.txt'
    outfile = output_dir + '/' + 'userInput.txt'

    lines = []
    with open(infile, 'r') as f:
        lines = f.read().splitlines()
        
    return lines
   
    EOS = '/EOS/eos_1Xsolar_cond.dat'
    planet.set_planet_equilibrium_temperature(star.temperature, star.radius, f=0.25, A_bond=0.3)
    temperature = round(planet.temperature/100) * 100
    pressure = (const.atmosphere.value * 100)  # 1 mb
    scattering_factor = 1.0    
    
    with open(outfile,'w') as wf:
        lines[3] = output_dir
        lines[5] = '/T_P/t_p_' + str(temperature) + 'K.dat'  
        lines[7] = str(EOS)
        lines[9] = cwd + '/ExoTransmit_Spectra/' + str(planet.name) + '.dat'
        lines[11] = str(planet.surface_G)                # m/s2
        lines[13] = str(planet.radius)                 # m
        lines[15] = str(star.radius)            # m
        lines[17] = str(pressure)               # pascals , leave at 0.0 if no cloud calculations
        lines[19] = str(scattering_factor)      # default = 1; 0= if no scattering
        
        for line in lines:
            wf.write(line)
            wf.write('\n')
            
def reduce_resolution(X,Y,h):
    num_points = int(len(X)/h)
    new_X = np.linspace(X[0],X[-1],num_points)
    
    new_Y = np.zeros(num_points)
    for i in range(num_points):
        for j in range(h):
            new_Y[i] += Y[i*h + j]
            
        new_Y[i] /= h
    
#==============================================================================
#     fig = plt.figure()
#     ax1 = fig.add_subplot(2,1,1)
#     ax1.plot(X,Y)
#     ax2 = fig.add_subplot(2,1,2)
#     ax2.plot(new_X,new_Y)
#==============================================================================
    return [new_X,new_Y]
    
