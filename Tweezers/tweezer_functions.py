from scipy import constants
import numpy as np
from IonChainTools import calcPositions,lengthScale
from scipy.optimize import fsolve
import itertools
from itertools import product
import pandas as pd
import re

#Constants in SI units
eps0 = constants.epsilon_0
m = 39.9626*constants.atomic_mass
c = constants.c
e = constants.e
hbar = constants.hbar
pi = np.pi


def potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist):
    '''
    Find the dipole potential of the optical tweezers beam for
    the given set of parameters at r=0 and z=0 -- without RWA


    omega_tweezer = angular frequency of tweezer laser beam [2*Pi x Hz]
    linewidths = linewidth of the given resonant transition of ion taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+), this is a list/array of some form with all relevent transitions
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data [2*Pi x Hz]
    beam_waists = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement
    '''
    p = []
    for i in range(len(linewidths)): 
        p.append( (-3.*P_opt*(c**2.)/((omega_res[i]**3.)*(beam_waist**2.))) * (linewidths[i]/((omega_res[i] - omega_tweezer)) +
                                          linewidths[i]/(omega_res[i] + omega_tweezer)) ) 
    pot = sum(p)
    return pot


def potentialRWA(omega_tweezer,linewidths,omega_res,P_opt,beam_waist):
    '''
    Find the potential of the optical tweezers beam for
    the given set of parameters at r=0 and z=0 -- with RWA


    omega_tweezer = angular frequency of tweezer laser beam [2*Pi x Hz]
    linewidths = linewidth of the given resonant transition of ion taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+), this is a list/array of some form with all relevent transitions
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data [2*Pi x Hz]
    beam_waists = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement
    '''
    p = []
    for i in range(len(linewidths)): 
        p.append( (-3*P_opt*(c**2)/((omega_res[i]**3)*(beam_waist**2))) * (linewidths[i]/((omega_res[i] - omega_tweezer)))  ) 
    pot = sum(p)
    return pot


def scattering(omega_tweezer,linewidths,omega_res,P_opt,beam_waist):
    '''
    Find the scattering of the optical tweezers beam (at r=0 and z=0) off of a given resonance
    for the given set of parameters -- without RWA
    omega_tweezer = angular frequency of tweezer laser beam [2*Pi x Hz]
    linewidths = linewidth of the given resonant transition taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+)
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data [2*Pi x Hz]
    beam_waists = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement
    '''
    s = []
    for i in range(len(linewidths)):
        s.append(((3*(c**2)*P_opt)/(hbar *pi* (omega_res[i]**3)*(beam_waist**2))) *((omega_tweezer/omega_res[i])**3)* (((linewidths[i]/(omega_res[i] - omega_tweezer))+
                                                                            (linewidths[i]/(omega_res[i] + omega_tweezer)))**2) )
    scat = sum(s)
    return scat


def scatteringRWA(omega_tweezer,linewidths,omega_res,P_opt,beam_waist):
    '''
    Find the scattering of the optical tweezers beam (at r=0 and z=0) off of a given resonance
    for the given set of parameters -- with RWA
    omega_tweezer = angular frequency of tweezer laser beam [2*Pi x Hz]
    linewidths = linewidth of the given resonant transition taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+)
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data [2*Pi x Hz]
    beam_waists = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement
    '''
    s = []
    for i in range(len(linewidths)):
        s.append( (3*P_opt*(c**2)/(hbar *pi* (omega_res[i]**3)*(beam_waist**2))) *((omega_tweezer/omega_res[i])**3) *((linewidths[i]/(omega_res[i] - omega_tweezer))**(2)) 
        )
    scat = sum(s)
    return scat

def rayleigh_length(w0, lambda_beam):
    """
    Calculate the Rayleigh length of a beam
    inputs: 
    w0 -- minimum beam waist [m]
    lambda_beam -- wavelength of the beam in meters [m]
    output:
    Rayleigh length in meters
    """
    return (pi * w0**2) / lambda_beam

def beam_propogation(w0, z_pos, lambda_beam):
    """
    Calculate the beam propagation 
    inputs:
    w0 --  minimum beam waist [m]
    z_pos -- list of z positions [m]
    lambda_beam -- wavelength of the beam [m]
    
    output:
    beam propagation as a function of z
    """
    rayleigh = rayleigh_length(w0, lambda_beam)
    return w0 * np.sqrt(1 + (z_pos / rayleigh)**2)

def intensity(x_pos,y_pos,P0, wz):
    """
    Calculate the intensity of a Gaussian beam at a given (r,z)
    inputs:
    x_pos -- list of x positions [m]
    y_pos -- list of y positions [m]
    P0 -- power [W]
    wz -- beam waist [m] calculated from beam_propogation function
    
    returns:
    Intensity in W/m^2
    """
    return (2 * P0 / (np.pi * wz**2))  * np.exp(-2 * x_pos**2 / wz**2) * np.exp(-2 * y_pos**2 / wz**2)

def intensity_TEM10(x_pos,y_pos,w0,wz,E0,n):
    """
    Calculate the intensity of a TEM10 beam at a given (r,z)
    inputs:
    x_pos -- list of x positions [m]
    y_pos -- list of y positions [m]
    P0 -- power [W]
    w0 -- minimum beam waist at z=0 [m], same one that is used in beam_propogation function
    wz -- beam waist [m] calculated from beam_propogation function

    returns:
    Intensity in W/m^2
    """
    return ( (( c * eps0 * n / 2 ) * (E0 **2 * w0**2 ))/ (wz**2) )* (8*x_pos**2 / (wz**2)) * np.exp(-2*x_pos**2 / wz**2) * np.exp(-2*y_pos**2 / wz**2)

def half_angle_beam_divergence(M_squared,w0,lambda_beam):
    """
    Calculate the half angle beam divergence
    inputs:
    M_squared -- beam quality factor
    w0 -- beam waist [m]
    lambda_beam -- wavelength of the beam [m]
    returns:
    half angle beam divergence in radians
    """
    return M_squared * lambda_beam/ (pi* w0)

def potential_position_dependent(omega_res,linewidths,omega_tweezer,intensity):
    """
    Calculate the potential of the optical tweezer beam at a specific position (r,z)
    
    
    
    """
    #rayleigh = rayleigh_length(FWHM,lambda_beam) 
    #beam_prop = beam_propogation(FWHM,rayleigh_length,z_pos)  
    #intens = intensity(P0,FWHM,beam_propogation,r)

    p = []
    for i in range(len(linewidths)): 
        p.append( (-3.*pi*(c**2.)/(2*(omega_res[i]**3.))) \
                                       * (linewidths[i]/((omega_res[i] - omega_tweezer)) +
                                          linewidths[i]/(omega_res[i] + omega_tweezer)) * intensity )
    pot = sum(p)
    return pot

def scattering_position_dependent(omega_res, linewidths, omega_tweezer, intensity):
    s = []
    for i in range(len(linewidths)):
        s.append(((3 * c**2 * intensity) / (2 * hbar * omega_res[i]**3)) *
                 (omega_tweezer / omega_res[i])**3 * 
                 ((linewidths[i] / (omega_res[i] - omega_tweezer) +
                   linewidths[i] / (omega_res[i] + omega_tweezer))**2))
    scat = sum(s)
    return scat

def pot_derivative_with_tweeze(x, omega_rf_axial, omega_tw_radial, tweezed_ion, displacement):
    """
    derivative of the potential energy of the ion chain, use this to find positions of ions in the trap
    This one is specifically for tweezing one ion in the chain
    inputs:

    x: list of ion positions
    omega_rf_axial: the rf axial trapping frequency [Hz]
    omega_tw_radial: the tweezer radial trapping frequency [Hz]
    tweezed_ion: Ion number for the tweezed ion MAKE ME A LIST
    displacement: distance between the tweezer beam center and position of the tweezed ion MAKE ME A LIST
    """
    N = len(x)
    A = 1/2 * m * omega_rf_axial**2
    B = (e**2) /(4 * pi * eps0)
    C = 1/2 * m * omega_tw_radial**2
    
    return [A*(x[m]) 
            - sum([B / (abs(x[m] - x[n])**2) for n in range(m) if x[m] != x[n]])  # Avoid division by zero
            + sum([B / (abs(x[m] - x[n])**2) for n in range(m+1, N) if x[m] != x[n]])  # Avoid division by zero
            #MAKE ME A LOOP SO THE LIST MAKES SENSE
            + C*(x[tweezed_ion] - displacement) if m == tweezed_ion else 0  # Only apply tweezer potential to the tweezed ion
            for m in range(N)]

def pot_derivative_with_2tweeze(x, omega_rf_axial, omega_tw_radial, tweezed_ion1,tweezed_ion2, displacement1,displacement2):
    """
    derivative of the potential energy of the ion chain, use this to find positions of ions in the trap
    This one is specifically for tweezing one ion in the chain
    inputs:

    x: list of ion positions
    omega_rf_axial: the rf axial trapping frequency [2*pi*Hz]
    omega_tw_radial: the tweezer radial trapping frequency [2*pi*Hz]
    tweezed_ion: Ion number for the tweezed ion
    displacement: distance between the tweezer beam center and position of the tweezed ion
    """
    N = len(x)
    A = 1/2 * m * omega_rf_axial**2
    B = (e**2) /(4 * pi * eps0)
    C = 1/2 * m * omega_tw_radial**2
    
    return [A*(x[m]) 
            - sum([B / (abs(x[m] - x[n])**2) for n in range(m) if x[m] != x[n]])  # Avoid division by zero
            + sum([B / (abs(x[m] - x[n])**2) for n in range(m+1, N) if x[m] != x[n]])  # Avoid division by zero
            + C*(x[tweezed_ion1] + displacement1) if m == tweezed_ion1 else 0  # Only apply tweezer potential to the tweezed ion1
            + C*(x[tweezed_ion2] - displacement2) if m == tweezed_ion2 else 0  # Only apply tweezer potential to the tweezed ion2|
            for m in range(N)]

def ion_spacing(N,omega_a):
    """
    Calculating the equilibrium positions of the ions in real units, as well as the distance between each ion
    inputs:
    N = number of ions
    omega_a = rf axial trap frequency [2*Pi x Hz]

    returns:
    list where first entry is list of equilibrium positions of ions in meters and second entry is list of distances between ions in meters
    """
    A = np.zeros((N, N))
    l = lengthScale(omega_a)
    ueq = calcPositions(N)*l
    
    diff_list = []
    for x, y in zip(ueq[0::], ueq[1::]):
        diff_list.append(y-x)
    return [ueq,diff_list]

def ion_spacing_tweezers(potential_from_tweezers,ionspacing,omega_rf_axial,omega_tw_radial,tweezed_ion,displacement):
    ueq = fsolve(potential_from_tweezers,ionspacing[0],args = (omega_rf_axial,omega_tw_radial,tweezed_ion,displacement))
    diff_list = []
    for x, y in zip(ueq[0::], ueq[1::]):
        diff_list.append(y-x)
    return [ueq,diff_list]

def ion_spacing_2_tweezers(pot_derivative_with_2tweeze,ionspacing,omega_rf_axial, omega_tw_radial, tweezed_ion1,tweezed_ion2, displacement1,displacement2):
    ueq = fsolve(pot_derivative_with_2tweeze,ionspacing[0],args = (omega_rf_axial, omega_tw_radial, tweezed_ion1,tweezed_ion2, displacement1,displacement2))
    diff_list = []
    for x, y in zip(ueq[0::], ueq[1::]):
        diff_list.append(y-x)
    return [ueq,diff_list]

def omega_tweezer_r(U,beam_waist,m):
    """
    Calculating the radial tweezer trap frequency (perpendicular to laser propogation) at r=0 and z=0
    given the tweezer potential U [J],
    the beam waist of the tweezer laser beam,
    and the mass of the ion

    U = potential created from the tweezer laser beam, from potential function at r=0 and z=0
    beam_waists = beam_waists = beamwaist of the tweezer laser beam
    m = mass of ion
    """
    return ((abs(U) * 4) / (m * (beam_waist)**2))**(1/2)

def omega_tweezer_a(U,beam_waist,tweezer_wavelength,m):
    """
     Calculating the axial tweezer trap frequency (along laser propogation) at r=0 and z=0
    given the tweezer potential U [J],
    the beam waist of the tweezer laser beam,
    and the mass of the ion

       U = potential created from the tweezer laser beam, from potential function
       beam_waists = beam_waists = beamwaist of the tweezer laser beam
       tweezer_wavelength = wavelgth of the tweezer laser beam
       m = mass of ion
       """
    return ((2*abs(U)/m)**(1/2)) * 1/((pi*(beam_waist**2)/tweezer_wavelength))

def TEM10_tweezer_optical_potential_to_trap_frequency_y(linewidths, omega_res,omega_tweezer, w0, m, P0):
    p = []
    for i in range(len(linewidths)): 
        p.append(np.sqrt(
            ((8*P0/(np.exp(1)*pi*w0**4))*((-2*pi*c**2) * (2/m) ) /(2*omega_res**3) ) * (linewidths[i]/((omega_res[i] - omega_tweezer)) +
                                            linewidths[i]/(omega_res[i] + omega_tweezer))
        ))
    pot = sum(p)
    return pot

def TEM10_tweezer_optical_potential_to_trap_frequency_x(linewidths, omega_res,omega_tweezer, w0, m, P0):
    p = []
    for i in range(len(linewidths)): 
        p.append(np.sqrt(
            ((40*P0/(np.exp(1)*pi*w0**4))*((-2*pi*c**2) * (2/m) )/(m*omega_res[i]**3*w0**2) * (linewidths[i]/((omega_res[i] - omega_tweezer)) +
                                            linewidths[i]/(omega_res[i] + omega_tweezer)))
        ))
    pot = sum(p)
    return pot

def mode_calc_r(m,omega_r_combined,ueq,N):
    
    """
    
    Hessian for ions in a pseudo-potential
    
    Inputs:
    N: number of ions 
    ueq -- list of equilibrium positions of ions (m)
    m -- mass of ion (kg)
    omega_r_combined -- combined radial frequency taking into account the rf potential as well as the tweezer potentials. 
                will look like array where each entry for untweezed ion is the rf radial frequency and each entry for the
                tweezed ions is sqrt(omega_tweezer^2 + omega_r_rf^2) (2*pi*Hz)
    omega_a -- axial trapping frequency created by rf potential (2*pi*Hz)
    
    Outputs: 
    modes -- list of tuples where each tuple is a mode (frequency [Hz], eigenvector)
    """
    A = np.zeros((N, N))
    coloumb = ((e**2) / (4 * pi * eps0))
    masses = np.array([m for _ in range(N)])
    for i in range(N):
        A[i][i] = (masses[i] * omega_r_combined[i]**2 - coloumb * sum(1 / (ueq[i] - ueq[m])**3 for m in range(0, i))
           - coloumb * sum(1 / (ueq[m] - ueq[i])**3 for m in range(i + 1, N))) * masses[i]
        for j in range(0, i):
            A[i][j] = (1/(ueq[i]-ueq[j])**3) * np.sqrt(masses[i])*np.sqrt(masses[j])*(coloumb)
        for j in range(i+1, N):
            A[i][j] = (1/ (ueq[j]-ueq[i])**3) *np.sqrt(masses[i])*np.sqrt(masses[j])*(coloumb)
    eigvals, eigvecs = np.linalg.eig(A) # this gives eigenvalues and eigenvectors
    freqs =( np.sqrt(1*eigvals))/(2*pi*m) #eigenvalue = spring constant k, so freq = sqrt(e-val)/(2*pi*m)
    
    
    scaledmodes = [(f, v) for f, v in zip(freqs, eigvecs.T)]
    scaledmodes = sorted(scaledmodes, key=lambda mode: mode[0],reverse=True)
    modes = []
    for f, scaledvec in scaledmodes:
        vec = np.array([scaledvec[i]/1 for i in range(len(eigvals))])
        vec = vec / np.sqrt(vec.dot(vec))
        modes.append((f, vec))
    return modes

def mode_calc_a(m,omega_a_combined,ueq,N):
    """
    Hessian for ions in a pseudo-potential
    Inputs:
    ueq -- equilibrium positions of the ions (m)
    m -- mass of ion (kg)
    omega_a_combined -- combined radial frequency taking into account the rf potential as well as the tweezer potentials. 
                will look like array where each entry for untweezed ion is the rf radial frequency and each entry for the
                tweezed ions is sqrt(omega_tweezer^2 + omega_a_rf^2) (2*pi*Hz)
    omega_a -- axial trapping frequency created by rf potential (2*pi*Hz)
    
    Outputs: 
    modes -- list of tuples where each tuple is a mode (frequency [Hz], eigenvector)
    """

    A = np.zeros((N, N))
    coloumb = ((e**2) / (4 * pi * eps0))
    masses = np.array([m for _ in range(N)])
    for i in range(N):
        A[i][i] = (masses[i] * omega_a_combined[i]**2 + coloumb * sum(2 / (ueq[i] - ueq[m])**3 for m in range(0, i))
           + coloumb * sum(2 / (ueq[m] - ueq[i])**3 for m in range(i + 1, N))) * masses[i]
        for j in range(0, i):
            A[i][j] = (-2/(ueq[i]-ueq[j])**3) * np.sqrt(masses[i])*np.sqrt(masses[j])*(coloumb)
        for j in range(i+1, N):
            A[i][j] = (-2/ (ueq[j]-ueq[i])**3) *np.sqrt(masses[i])*np.sqrt(masses[j])*(coloumb)

    eigvals, eigvecs = np.linalg.eig(A) # this gives eigenvalues and eigenvectors
    freqs =( np.sqrt(1*eigvals))/(2*pi*m) #eigenvalue = spring constant k, so freq = sqrt(e-val)/(2*pi*m)

    scaledmodes = [(f, v) for f, v in zip(freqs, eigvecs.T)]
    scaledmodes = sorted(scaledmodes, key=lambda mode: mode[0],reverse=False)
    modes = []
    for f, scaledvec in scaledmodes:
        vec = np.array([scaledvec[i]/1 for i in range(len(eigvals))])
        vec = vec / np.sqrt(vec.dot(vec))
        modes.append((f, vec))
    return modes

def eta(mode_structure,qubit_wavelength,N):
    """input:
    mode structure as output from mode_calc_r or mode_calc_a
    N = number of ions
    qubit_wavelength = wavelength of qubit transition [m] (729e-9 for Ca)
    output:
    eta values for each mode and ion instead of just the eigenvectors
    """

    eta = []
    for mode in mode_structure:
        eta.append([mode[1][i] * (2 * pi / qubit_wavelength) * np.sqrt(hbar / (2 * m * mode[0])) for i in range(N)])
    return eta

def combined_frequencies(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a):
    '''
    takes in rf and tweezer trap frequencies and adds together frequencies in quadruture
    radial modes will be effected by either the radial and axial tweezer directions (in BladeRunner setup)
    axial modes will be effected by only tweezer radial
    
    inputs:
    N = number of ions
    tweezed_ions = list of which ions are getting tweezed
    w_tweezer_r = radial trapping frequency of tweezer [2*Pi x Hz]
    w_tweezer_a = axial trapping frequency of tweezer [2*Pi x Hz]
    w_rf_r = radial rf trapping frequency [2*Pi x Hz]
    w_rf_a = axial rf trapping frequency [2*Pi x Hz]
    
    returns: array of potential combined trapping frequencies
                [0] is radial rf and radial tweezer
                [1] is radial rf and axial tweezer
                [2] is axial rf and radial tweezer
    
    '''

    omeg_tweezer_r = np.zeros(N)
    omeg_tweezer_a = np.zeros(N)
    omeg_tweezer_r[tweezed_ions] = w_tweezer_r
    omeg_tweezer_a[tweezed_ions] = w_tweezer_a

    omeg_rf_r = w_rf_r * np.ones(N) 
    omeg_rf_a = w_rf_a * np.ones(N)

    omega_combined_rr = np.sqrt(omeg_rf_r**2 + omeg_tweezer_r**2)
    omega_combined_ra = np.sqrt(omeg_rf_r**2 + omeg_tweezer_a)
    omega_combined_ar = np.sqrt(omeg_rf_a**2 + omeg_tweezer_r**2)
    
    return np.array([omega_combined_rr,omega_combined_ra,omega_combined_ar])

def trapping_ratios(w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a):
    
    """
    takes in tweezer and rf trap frequencies and outputs various ratios of frequencies in case this turns out to be helpful
    
    inputs:
    w_tweezer_r = radial trapping frequency of tweezer [2*Pi x Hz]
    w_tweezer_a = axial trapping frequency of tweezer [2*Pi x Hz]
    w_rf_r = radial rf trapping frequency [2*Pi x Hz]
    w_rf_a = axial rf trapping frequency [2*Pi x Hz]
    
    returns: array of ratios of omegas
    """
    
    tweezer_r_to_rf_ratio = w_tweezer_r / w_rf_r
    tweezer_r_to_axial_ratio = w_tweezer_r / w_rf_a
    tweezer_a_to_rf_ratio = w_tweezer_a / w_rf_a
    
    return np.array([tweezer_r_to_rf_ratio,tweezer_r_to_axial_ratio])


def individual_freqs_to_mode_vectors(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a,ueq):
    """
    takes in output from combined_frequencies and outputs new radial modes
    
    inputs:
    N = number of ions
    tweezed_ions = list of which ions are getting tweezed
    w_tweezer_r = radial trapping frequency of tweezer [2*Pi x Hz]
    w_tweezer_a = axial trapping frequency of tweezer [2*Pi x Hz]
    w_rf_r = radial rf trapping frequency [2*Pi x Hz]
    w_rf_a = axial rf trapping frequency [2*Pi x Hz]
    
    returns: modes from mode_calc_r
    
    """
    combined_freqs = combined_frequencies(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a)
    omega_r_combined = combined_freqs[0]
    #omega_a = w_rf_a
    return mode_calc_r(m,omega_r_combined,ueq,N)

def individual_freqs_to_mode_vectors_axial(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a,ueq):
    """
    takes in output from combined_frequencies and outputs new axial modes
    
    inputs:
    N = number of ions
    tweezed_ions = list of which ions are getting tweezed
    w_tweezer_r = radial trapping frequency of tweezer [2*Pi x Hz]
    w_tweezer_a = axial trapping frequency of tweezer [2*Pi x Hz]
    w_rf_r = radial rf trapping frequency [2*Pi x Hz]
    w_rf_a = axial rf trapping frequency [2*Pi x Hz]
    
    returns: modes from mode_calc_r
    
    """
    combined_freqs = combined_frequencies(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a)
    omega_a_combined = combined_freqs[2]
    #omega_a = w_rf_a
    return mode_calc_a(m,omega_a_combined,ueq,N)

def individual_freqs_to_mode_vectors_radial_weak(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a):
    """
    I'm unsure what this function even is...
    takes in output from combined_frequencies and outputs new radial modes (but with the weak trappping from tweezers)
    
    inputs:
    N = number of ions
    tweezed_ions = list of which ions are getting tweezed
    w_tweezer_r = radial trapping frequency of tweezer [2*Pi x Hz]
    w_tweezer_a = axial trapping frequency of tweezer [2*Pi x Hz]
    w_rf_r = radial rf trapping frequency [2*Pi x Hz]
    w_rf_a = axial rf trapping frequency [2*Pi x Hz]
    
    returns: modes from mode_calc_r
    
    """
    combined_freqs = combined_frequencies(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a)
    omega_r_combined = combined_freqs[1]
    omega_a = w_rf_a
    return mode_calc_r(m,omega_r_combined,omega_a)

def tweezer_optical_potential_to_trap_frequency(tweezer_wavelength,linewidths,omega_res,P_opt,beam_waist,m,U):
    """
    takes in physical parameters of calcium ion and tweezer beam and outputs expected tweezer trap frequency
    
    Inputs:
    tweezer_wavelength = tweezer wavelength [m]
    linewidths = linewidth of the given resonant transition taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+)
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data [2*Pi x Hz]
    beam_waist = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement
    U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist) OR 
        potential_position_dependent(omega_res,linewidths,omega_tweezer,intensity)
    
    outputs:
    array of radial and axial tweezer trap frequencies [2*Pi x Hz]
    
    """
    omega_tweezer = 2*pi*c/tweezer_wavelength
    
    
    w_tweezer_r =  omega_tweezer_r(U,beam_waist,m)
    w_tweezer_a = omega_tweezer_a(U,beam_waist,tweezer_wavelength,m)
    return np.array([w_tweezer_r,w_tweezer_a])



def physical_params_to_radial_mode_vectors(N,ueq,tweezed_ions,tweezer_wavelength,linewidths,omega_res,w_rf_a,w_rf_r,P_opt,beam_waist,m,U):
    """
    takes in physical parameters of tweezer beam and calcium ion as well as rf 
    trapping parameters to output combined radial modes
    
    Inputs:
    tweezer_wavelength = tweezer wavelength [m]
    linewidths = linewidth of the given resonant transition taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+)
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data [2*Pi x Hz]
    beam_waist = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement  
    N = number of ions
    tweezed_ions = list of which ions are getting tweezed
    w_tweezer_r = radial trapping frequency of tweezer [2*Pi x Hz]
    w_tweezer_a = axial trapping frequency of tweezer [2*Pi x Hz]
    w_rf_r = radial rf trapping frequency [2*Pi x Hz]
    w_rf_a = axial rf trapping frequency [2*Pi x Hz]
    U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist) OR 
        potential_position_dependent(omega_res,linewidths,omega_tweezer,intensity)
    
    outputs:
    modes from mode_calc_r, frequencies in Hz (not angular)
    
    """
   
    omega_tweezer = 2*pi*c/tweezer_wavelength
    
    #U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist)
    w_tweezer_r =  omega_tweezer_r(U,beam_waist,m)
    w_tweezer_a = omega_tweezer_a(U,beam_waist,tweezer_wavelength,m)
    return individual_freqs_to_mode_vectors(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a,ueq)

def physical_params_to_axial_mode_vectors(N,ueq,tweezed_ions,tweezer_wavelength,linewidths,omega_res,w_rf_a,w_rf_r,P_opt,beam_waist,m,U):
    """
    takes in physical parameters of tweezer beam and calcium ion as well as rf 
    trapping parameters to output combined axial modes
    
    Inputs:
    tweezer_wavelength = tweezer wavelength [m]
    linewidths = linewidth of the given resonant transition taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+)
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data [2*Pi x Hz]
    beam_waist = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement  
    N = number of ions
    tweezed_ions = list of which ions are getting tweezed
    w_tweezer_r = radial trapping frequency of tweezer [2*Pi x Hz]
    w_tweezer_a = axial trapping frequency of tweezer [2*Pi x Hz]
    w_rf_r = radial rf trapping frequency [2*Pi x Hz]
    w_rf_a = axial rf trapping frequency [2*Pi x Hz]
    U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist) OR 
        potential_position_dependent(omega_res,linewidths,omega_tweezer,intensity)
    
    outputs:
    modes from mode_calc_r, frequencies in Hz (not angular)
    
    """
   
    omega_tweezer = 2*pi*c/tweezer_wavelength
    
    #U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist)
    w_tweezer_r =  omega_tweezer_r(U,beam_waist,m)
    w_tweezer_a = omega_tweezer_a(U,beam_waist,tweezer_wavelength,m)
    return individual_freqs_to_mode_vectors_axial(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a,ueq)

def physical_params_to_radial_mode_vectors_weak(N,tweezed_ions,tweezer_wavelength,linewidths,omega_res,w_rf_a,w_rf_r,P_opt,beam_waist,m,U):
    """
    unsure what this one is too
    takes in physical parameters of tweezer beam and calcium ion as well as rf 
    trapping parameters to output combined radial modes (from weak tweezer trap)
    
    Inputs:
    tweezer_wavelength = tweezer wavelength [m]
    linewidths = linewidth of the given resonant transition taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+)
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data [2*Pi x Hz]
    beam_waist = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement  
    N = number of ions
    tweezed_ions = list of which ions are getting tweezed
    w_tweezer_r = radial trapping frequency of tweezer [2*Pi x Hz]
    w_tweezer_a = axial trapping frequency of tweezer [2*Pi x Hz]
    w_rf_r = radial rf trapping frequency [2*Pi x Hz]
    w_rf_a = axial rf trapping frequency [2*Pi x Hz]
    U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist) OR 
        potential_position_dependent(omega_res,linewidths,omega_tweezer,intensity)
    
    outputs:
    modes from mode_calc_r, frequencies in Hz (not angular)
    
    """
   
    omega_tweezer = 2*pi*c/tweezer_wavelength
    
    #U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist)
    w_tweezer_r =  omega_tweezer_r(U,beam_waist,m)
    w_tweezer_a = omega_tweezer_a(U,beam_waist,tweezer_wavelength,m)
    return individual_freqs_to_mode_vectors_radial_weak(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a)

#this section down here is functions for the sideband cooling calculations

def tweezer_combos_full_radial(
    omega_tweezer,
    linewidths,
    omega_res,
    m,
    mode_calc_r,
    N_list,
    f_rf_r,
    P_opt,
    w0,
    max_tweezed=1,
):
    """
    Compute all tweezer combinations and corresponding radial mode frequencies,
    including power dependence and tweezed/untweezed mode separation.

    mode_calc_r is expected to return a list of tuples:
        [(freq1, eigvec1), (freq2, eigvec2), ...]
    This version ensures the dataframe has Mode0_eigvec, Mode1_eigvec, ... Mode{N-1}_eigvec
    and Mode0_freq, Mode1_freq, ... Mode{N-1}_freq for each row (missing entries filled with NaN).
    """


    # --- Normalize inputs ---
    if np.isscalar(N_list):
        N_list = [int(N_list)]
    if np.isscalar(P_opt):
        P_opt = [P_opt]

    pi = np.pi
    rows = []

    # --- Loop over number of ions ---
    for N in N_list:
        # Generate all possible tweezer combinations
        all_combos = []
        for r in range(0, max_tweezed + 1):
            all_combos.extend(itertools.combinations(range(N), r))

        # RF trap setup
        w_rf_r = f_rf_r * 2 * pi
        w_rf_r_list = np.full(N, w_rf_r)
        ueq = ion_spacing(N, f_rf_r)[0]

        # --- Loop over optical powers ---
        for P_total in P_opt:
            for tweezed_positions in all_combos:
                n_tweezed = len(tweezed_positions)
                P_per = P_total / n_tweezed if n_tweezed > 0 else 0.0

                # Compute tweezer potential for this configuration
                pot = potential(omega_tweezer, linewidths, omega_res, P_per, w0)
                w_tw_r = omega_tweezer_r(pot, w0, m)

                # Combine tweezed and untweezed radial frequencies
                combo = np.array([
                    np.sqrt(w_tw_r**2 + w_rf_r_list[i]**2) if i in tweezed_positions else w_rf_r_list[i]
                    for i in range(N)
                ])

                # --- Compute radial modes ---
                modes = mode_calc_r(m, combo, ueq, N)

                # Extract frequencies and eigenvectors
                freqs = np.array([f for f, v in modes], dtype=float) if len(modes) else np.array([], dtype=float)
                if len(modes):
                    eigvecs = np.vstack([np.ravel(v) for f, v in modes])  # shape (n_modes, N)
                else:
                    eigvecs = np.empty((0, N))

                # --- Initialize row with shared info ---
                row = {
                    "N": N,
                    "Tweezed ions": tweezed_positions,
                    "P_per_tweezer (W)": P_per,
                    "Combined radial frequencies": combo,
                }

                # --- Ensure columns for all possible modes up to N exist per row ---
                # Fill Mode{i}_freq and Mode{i}_eigvec for i in [0, N-1]
                for mode_index in range(N):
                    # frequency
                    if mode_index < len(freqs):
                        row[f"Mode{mode_index}_freq"] = float(freqs[mode_index])
                    else:
                        row[f"Mode{mode_index}_freq"] = np.nan

                    # eigenvector (length N) or NaN array
                    if mode_index < eigvecs.shape[0]:
                        row[f"Mode{mode_index}_eigvec"] = np.ravel(eigvecs[mode_index]).astype(float)
                    else:
                        # use full-length nan array to keep shape consistent
                        row[f"Mode{mode_index}_eigvec"] = np.full(N, np.nan, dtype=float)

                # --- Store completed row ---
                rows.append(row)

    return pd.DataFrame(rows)
# ...existing code...
def build_mode_series_and_combinations(df, max_modes=None):
    """
    Extract per-mode lists of (df_index, eigvec_array) from df.

    Returns a dict with:
      - mode_series: mapping mode_index -> original pandas Series (unchanged)
      - mode_lists:  mapping mode_index -> list of tuples (df_index, np.ndarray(eigvec))
    If max_modes is set, only modes with index < max_modes are returned.
    """


    # find Mode{i}_eigvec columns sorted by i
    mode_cols = sorted(
        [c for c in df.columns if re.match(r"^Mode\d+_eigvec$", c)],
        key=lambda c: int(re.match(r"Mode(\d+)_eigvec$", c).group(1)),
    )
    mode_indices = [int(re.match(r"Mode(\d+)_eigvec$", c).group(1)) for c in mode_cols]

    if max_modes is not None:
        mode_indices = [i for i in mode_indices if i < int(max_modes)]

    # keep the original Series for convenience
    mode_series = {i: df[f"Mode{i}_eigvec"] for i in mode_indices}

    # build lists of (original_index, np.array(value)) for each mode
    mode_lists = {}
    for i in mode_indices:
        col = f"Mode{i}_eigvec"
        items = []
        if col in df.columns:
            for idx in df.index:
                val = df.at[idx, col]
                # convert to a numeric numpy array (works if stored as list/ndarray/scalar)
                try:
                    arr = np.asarray(val, dtype=float)
                except Exception:
                    # fall back to object array if conversion fails
                    arr = np.atleast_1d(val)
                items.append((idx, arr))
        else:
            # column missing -> empty arrays for each row (keeps index correspondence)
            for idx in df.index:
                items.append((idx, np.array([], dtype=float)))
        mode_lists[i] = items

    return {"mode_series": mode_series, "mode_lists": mode_lists}

def condense_by_min_abs(data):
    """
    Condense each tuple (idx_group, combo, array) into
    (idx_group, combo, value) where value is the element with the smallest
    absolute magnitude, but preserve its original sign.
    """
    condensed = []
    for idx_group, combo, arr in data:
        # choose element with smallest abs() but keep real sign
        min_val = min(arr, key=lambda x: abs(x))
        condensed.append((idx_group, combo, min_val))
    return condensed

def filter_by_max_min_abs(data):
    """
    Filter condensed tuples so that only those whose stored value has the
    largest absolute magnitude remain. Original sign preserved.
    """
    if not data:
        return []
    
    max_abs = max(abs(t[2]) for t in data)
    return [t for t in data if abs(t[2]) == max_abs]




def combine_lists(*lists):
    """
    Fully general N-dimensional version that only allows
    element index combinations with unique indices.
    """

    # Number of lists (N)
    N = len(lists)

    # Length of the vectors (K)
    K = len(lists[0][0][1])

    # Collect all distinct original indices
    keys = [idx for idx, _ in lists[0]]

    # Map each list by idx for fast lookup
    idx_maps = []
    for lst in lists:
        idx_maps.append({idx: arr for idx, arr in lst})

    # Storage for output groups
    groups = {key: [] for key in keys}

    # Loop over each index group
    for key in keys:

        # Grab the vector chosen from each list for this group
        chosen = [idx_maps[m][key] for m in range(N)]

        # Sweep all element-index choices independently
        for elem_choices in product(range(K), repeat=N):

            # NEW RULE: require all unique indices
            if len(set(elem_choices)) != N:
                continue

            # Build output vector element-wise
            values = np.array([
                chosen[m][elem_choices[m]]
                for m in range(N)
            ])

            # Store tuple
            groups[key].append(
                (key, elem_choices, values)
            )

    return groups

def select_global_max_min_abs(groups, tol=1e-12):
    """
    Given a list of lists where each inner list contains tuples
    (idx_group, combo, value),
    return all tuples whose |value| equals the maximum absolute value
    across the entire dataset, allowing for floating-point tolerance.
    """

    # Flatten everything into one list of tuples
    all_tuples = [t for group in groups for t in group]

    if not all_tuples:
        return []

    # Compute global maximum |value|
    global_max = max(abs(t[2]) for t in all_tuples)

    # Collect all tuples that match this max within tolerance
    winners = [
        t for t in all_tuples
        if abs(abs(t[2]) - global_max) < tol
    ]

    return winners
