from scipy import constants
import numpy as np


#Constants in SI units
eps0 = constants.epsilon_0
m = 40.07*constants.atomic_mass
c = constants.c
e = constants.e
hbar = constants.hbar
pi = np.pi

def potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waists):
    '''
    Find the potential of the optical tweezers beam for
    the given set of parameters -- without RWA


    omega_tweezer = angular frequency of tweezer laser beam
    linewidths = linewidth of the given resonant transition
        (ex, S1/2 to D5/2)
    P_opt = total optical power of tweezer laser beam
    beam_waists = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system
    '''

    return (-3*pi*(c**2)/omega_tweezer**3) * (linewidths/((omega_res - omega_tweezer)) +
                                          linewidths/(omega_res + omega_tweezer)) * ((P_opt)/(beam_waists**2))


def potentialRWA(omega_tweezer,linewidths,omega_res,P_opt,beam_waists):
    '''
    Find the potential of the optical tweezers beam for the
    given set of parameters -- Using RWA


    omega_tweezer = angular frequency of tweezer laser beam
    linewidths = linewidth of the given resonant transition
        (ex, S1/2 to D5/2)
    P_opt = total optical power of tweezer laser beam
    beam_waists = beamwaist of the tweezer laser beam
    given its frequency and the NA of the system
    '''
    return (3*c**2/omega_tweezer**3) * (linewidths/(omega_res - omega_tweezer)) * P_opt/(beam_waists**2)


def scattering(omega_tweezer,linewidths,omega_res,P_opt,beamwaists):
    '''
    Find the scattering of the optical tweezers beam off of a given resonance
    for the given set of parameters -- without RWA
    omega_tweezer = angular frequency of tweezer laser beam
    linewidths = linewidth of the given resonant transition
        (ex, S1/2 to D5/2)
    P_opt = total optical power of tweezer laser beam
    beam_waists = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system
    '''
    return ((3*pi*(c**2))/(hbar * (omega_tweezer**3))) *((omega_tweezer/omegares)**3)* (((linewidths/(omega_res - omega_tweezer))+
                                                                            (linewidths/(omegares + omega_tweezer)))**2) * ((P_opt)/(beam_waists**2))
def epsilon(d,wx,m):
    """
    Dimensionless parameter describing characteristic energy(?) scales of our ion chain
    """

    return ((e**2)/(4*pi*eps0*(d**3)*m*(wx**2)))**(1/2)

def v(wx,wxt):
    """
    Dimensionless ratio of radial trapping frequency of tweezer /rf trap
    """

    return wxt/wx

def omega_radial(U,beam_waists,m):
    """
    Calculating the radial tweezer trap frequency
    given the tweezer potential (or trap depth) U,
    the beam waist of the tweezer laser beam,
    and the mass of the ion

    U = potential created from the tweezer laser beam
    beam_waists = beam_waists = beamwaist of the tweezer laser beam
    m = mass of ion
    """
    return ((abs(U) * 4) / (m * (beam_waists)**2))**(1/2)

def omega_axial(U,beam_waists,tweezer_wavelength,m):
    """
    Calculating the radial tweezer trap frequency
       given the tweezer potential (or trap depth) U,
       the beam waist of the tweezer laser beam,
       the wavelength of the tweezer laser beam,
       and the mass of the ion

       U = potential created from the tweezer laser beam
       beam_waists = beam_waists = beamwaist of the tweezer laser beam
       tweezer_wavelength = wavelgth of the tweezer laser beam
       m = mass of ion
       """
    return ((2*abs(U)/m)**(1/2)) * 1/((pi*(beam_waist**2)/tweezer_wavelength))