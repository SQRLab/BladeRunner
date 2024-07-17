from scipy import constants
import numpy as np
from IonChainTools import calcPositions,lengthScale

#Constants in SI units
eps0 = constants.epsilon_0
m = 39.9626*constants.atomic_mass
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


def scatteringRWA(omegatweezer,linewidths,omegares,Popt,beamwaists):
    '''Find the scattering of the optical tweezers beam off of a given resonance
    for the given set of parameters -- Using RWA
    omegatweezer = angular frequency of tweezer laser beam
    linewidths = linewidth of the given resonant transition (ex, S1/2 to D5/2)
    Detuning = omegatweezer - omegaresonance where omegaresonance is the angular frequency corresponding to linewidth
    counterrotating = omegatweezer+omegaresonance 
    Popt = total optical power of tweezer laser beam
    beamwaists = beamwaist of the tweezer laser beam given its frequency and the NA of our system'''
    return (3*c**2/(hbar * (omegatweezer**3))) * ((linewidths/(omegares - omegatweezer))**(2)) * Popt/(beamwaists**2)


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
    return ((2*abs(U)/m)**(1/2)) * 1/((pi*(beam_waists**2)/tweezer_wavelength))

def mode_calc_r(m,omega_r,omega_a):
    N = len(omega_r)
    A = np.zeros((N, N))
    l = lengthScale(omega_a)
    ueq = calcPositions(N)*l
    coloumb = ((e**2) / (4 * pi * eps0))
    masses = np.array([m for _ in range(N)])
    for i in range(N):
        A[i][i] = (masses[i] * omega_r[i]**2 - coloumb * sum(1 / (ueq[i] - ueq[m])**3 for m in range(0, i))
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

def mode_calc_a(m,omega_a):
    N = len(omega_a)
    A = np.zeros((N, N))
    l = lengthScale(omega_a)
    ueq = calcPositions(N)*l
    coloumb = ((e**2) / (4 * pi * eps0))
    masses = np.array([m for _ in range(N)])
    for i in range(N):
        A[i][i] = (masses[i] * omega_a[i]**2 + coloumb * sum(2 / (ueq[i] - ueq[m])**3 for m in range(0, i))
           + coloumb * sum(2 / (ueq[m] - ueq[i])**3 for m in range(i + 1, N))) * masses[i]
        for j in range(0, i):
            A[i][j] = (-2/(ueq[i]-ueq[j])**3) * np.sqrt(masses[i])*np.sqrt(masses[j])*(coloumb)
        for j in range(i+1, N):
            A[i][j] = (-2/ (ueq[j]-ueq[i])**3) *np.sqrt(masses[i])*np.sqrt(masses[j])*(coloumb)

    eigvals, eigvecs = np.linalg.eig(A) # this gives eigenvalues and eigenvectors
    freqs =( np.sqrt(1*eigvals))/(2*pi*m) #eigenvalue = spring constant k, so freq = sqrt(e-val)/(2*pi*m)
    print(freqs)
    
    
    scaledmodes = [(f, v) for f, v in zip(freqs, eigvecs.T)]
    scaledmodes = sorted(scaledmodes, key=lambda mode: mode[0],reverse=False)
    modes = []
    for f, scaledvec in scaledmodes:
        vec = np.array([scaledvec[i]/1 for i in range(len(eigvals))])
        vec = vec / np.sqrt(vec.dot(vec))
        modes.append((f, vec))
    return modes


