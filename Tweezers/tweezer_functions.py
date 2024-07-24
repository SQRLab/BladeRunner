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


def potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist):
    '''
    Find the potential of the optical tweezers beam for
    the given set of parameters at r=0 and z=0 -- without RWA


    omega_tweezer = angular frequency of tweezer laser beam
    linewidths = linewidth of the given resonant transition of ion taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+), this is a list/array of some form with all relevent transitions
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data
    beam_waists = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement
    '''
    p = []
    for i in range(len(linewidths)): 
        p.append( (-3*P_opt*(c**2)/((omega_res[i]**3)*(beam_waist**2))) * (linewidths[i]/((omega_res[i] - omega_tweezer)) +
                                          linewidths[i]/(omega_res[i] + omega_tweezer)) ) 
    pot = sum(p)
    return pot


def potentialRWA(omega_tweezer,linewidths,omega_res,P_opt,beam_waist):
    '''
    Find the potential of the optical tweezers beam for
    the given set of parameters at r=0 and z=0 -- with RWA


    omega_tweezer = angular frequency of tweezer laser beam
    linewidths = linewidth of the given resonant transition of ion taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+), this is a list/array of some form with all relevent transitions
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data
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
    omega_tweezer = angular frequency of tweezer laser beam
    linewidths = linewidth of the given resonant transition taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+)
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data
    beam_waists = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement
    '''
    s = []
    for i in range(len(linewidths)):
        s.append(((3*(c**2)*P_opt)/(hbar *pi* (omega_res[i]**3)*(beam_waist**2))) *((omega_tweezer/omegares[i])**3)* (((linewidths[i]/(omega_res[i] - omega_tweezer))+
                                                                            (linewidths[i]/(omegares[i] + omega_tweezer)))**2) )
    scat = sum(s)
    return scat


def scatteringRWA(omega_tweezer,linewidths,omega_res,P_opt,beam_waist):
    '''
    Find the scattering of the optical tweezers beam (at r=0 and z=0) off of a given resonance
    for the given set of parameters -- with RWA
    omega_tweezer = angular frequency of tweezer laser beam
    linewidths = linewidth of the given resonant transition taken from NIST database in angular frequency units 
        (ex, S1/2 to P1/2 and S1/2 to P3/2 for 40Ca+)
    P_opt = total optical power of tweezer laser beam
    omega_res = angular frequency of resonant transition, also based off NIST data
    beam_waists = beamwaist of the tweezer laser beam
    given its frequency and the NA of our system or from measurement
    '''
    s = []
    for i in range(len(linewidths)):
        s.append( (3*P_opt*(c**2)/(hbar *pi* (omega_res[i]**3)*(beam_waist**2))) *((omega_tweezer/omega_res[i])**3) *((linewidths[i]/(omega_res[i] - omega_tweezer))**(2)) 
        )
    scat = sum(s)
    return scat





def omega_tweezer_r(U,beam_waists,m):
    """
    Calculating the radial tweezer trap frequency (perpendicular to laser propogation) at r=0 and z=0
    given the tweezer potential U [J],
    the beam waist of the tweezer laser beam,
    and the mass of the ion

    U = potential created from the tweezer laser beam, from potential function
    beam_waists = beam_waists = beamwaist of the tweezer laser beam
    m = mass of ion
    """
    return ((abs(U) * 4) / (m * (beam_waists)**2))**(1/2)

def omega_tweezer_a(U,beam_waists,tweezer_wavelength,m):
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
    return ((2*abs(U)/m)**(1/2)) * 1/((pi*(beam_waists**2)/tweezer_wavelength))

def mode_calc_r(m,omega_r_combined,omega_a):
    
    """
    
    Hessian for ions in a pseudo-potential
    
    Inputs:
    m -- mass of ion
    omega_r_combined -- combined radial frequency taking into account the rf potential as well as the tweezer potentials. 
                will look like array where each entry for untweezed ion is the rf radial frequency and each entry for the
                tweezed ions is sqrt(omega_tweezer^2 + omega_r_rf^2)
    omega_a -- axial trapping frequency created by rf potential
    
    
    """
    N = len(omega_r_combined)
    A = np.zeros((N, N))
    l = lengthScale(omega_a)
    ueq = calcPositions(N)*l
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

def mode_calc_a(m,omega_a,omega_a_combined):
    """
    
    Hessian for ions in a pseudo-potential
    
    Inputs:
    m -- mass of ion
    omega_a_combined -- combined radial frequency taking into account the rf potential as well as the tweezer potentials. 
                will look like array where each entry for untweezed ion is the rf radial frequency and each entry for the
                tweezed ions is sqrt(omega_tweezer^2 + omega_a_rf^2)
    omega_a -- axial trapping frequency created by rf potential
    
    
    """
    N = len(omega_a)
    A = np.zeros((N, N))
    l = lengthScale(omega_a)
    ueq = calcPositions(N)*l
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
    print(freqs)
    
    
    scaledmodes = [(f, v) for f, v in zip(freqs, eigvecs.T)]
    scaledmodes = sorted(scaledmodes, key=lambda mode: mode[0],reverse=False)
    modes = []
    for f, scaledvec in scaledmodes:
        vec = np.array([scaledvec[i]/1 for i in range(len(eigvals))])
        vec = vec / np.sqrt(vec.dot(vec))
        modes.append((f, vec))
    return modes


