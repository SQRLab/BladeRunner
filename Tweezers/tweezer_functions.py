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

def rayleigh_length(FWHM,lambda_beam):
    """
    calculate the rayleigh length of a beam
    inputs: 
    
    lambda_beam -- wavelength of the beam
    
    output:
    rayleigh length
    
    """
    
    return (pi * FWHM**2)/lambda_beam

def beam_propogation(FWHM,rayleigh_length,z_pos):
    """
    calculate the beam propogation 
    inputs:
    FWHM --  the full width half max of a gaussian beam
    rayleigh_length -- calculated from rayleigh_length function
    z_pos -- list of z positions
    
    output:
    beam propogation as a function of z
    
    
    """
    
    return FWHM * np.sqrt(1 + (z_pos/rayleigh_length)**2)

def intensity(P0,FWYM,beam_propogation,r):
    """
    calculate the intensity of the beam at a given (r,z)
    inputs:
    FWHM -- the full width half max of a gaussian beam
    beam_propogation -- beam propogation as a function of z
    
    """
    
    return (2*P0/pi*FWHM**2) * (FWHM / beam_propogation)**2 * np.exp(-2*r**2 / beam_propogation**2)


def potential_position_dependent(omega_res,linewidths,omega_tweezer,intensity):
    """
    Calculate the potential of the optical tweezer beam at a specific position (r,z)
    
    
    
    """
    p = []
    for i in range(len(linewidths)): 
        p.append( (-3.*pi*(c**2.)/(2*(omega_res[i]**3.))) \
                                       * (linewidths[i]/((omega_res[i] - omega_tweezer)) +
                                          linewidths[i]/(omega_res[i] + omega_tweezer)) * intensity )
    pot = sum(p)
    return pot




def omega_tweezer_r(U,beam_waist,m):
    """
    Calculating the radial tweezer trap frequency (perpendicular to laser propogation) at r=0 and z=0
    given the tweezer potential U [J],
    the beam waist of the tweezer laser beam,
    and the mass of the ion

    U = potential created from the tweezer laser beam, from potential function
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
    N = len(omega_a_combined)
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


def individual_freqs_to_mode_vectors(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a):
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
    omega_a = w_rf_a
    return mode_calc_r(m,omega_r_combined,omega_a)

def individual_freqs_to_mode_vectors_axial(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a):
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
    omega_a = w_rf_a
    return mode_calc_a(m,omega_a,omega_a_combined)

def individual_freqs_to_mode_vectors_radial_weak(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a):
    """
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

def tweezer_optical_potential_to_trap_frequency(tweezer_wavelength,linewidths,omega_res,P_opt,beam_waist,m):
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
    
    outputs:
    array of radial and axial tweezer trap frequencies [2*Pi x Hz]
    
    """
    omega_tweezer = 2*pi*c/tweezer_wavelength
    
    U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist)
    w_tweezer_r =  omega_tweezer_r(U,beam_waist,m)
    w_tweezer_a = omega_tweezer_a(U,beam_waist,tweezer_wavelength,m)
    return np.array([w_tweezer_r,w_tweezer_a])



def physical_params_to_radial_mode_vectors(N,tweezed_ions,tweezer_wavelength,linewidths,omega_res,w_rf_a,w_rf_r,P_opt,beam_waist,m):
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
    
    outputs:
    modes from mode_calc_r, frequencies in Hz (not angular)
    
    """
   
    omega_tweezer = 2*pi*c/tweezer_wavelength
    
    U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist)
    w_tweezer_r =  omega_tweezer_r(U,beam_waist,m)
    w_tweezer_a = omega_tweezer_a(U,beam_waist,tweezer_wavelength,m)
    return individual_freqs_to_mode_vectors(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a)

def physical_params_to_axial_mode_vectors(N,tweezed_ions,tweezer_wavelength,linewidths,omega_res,w_rf_a,w_rf_r,P_opt,beam_waist,m):
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
    
    outputs:
    modes from mode_calc_r, frequencies in Hz (not angular)
    
    """
   
    omega_tweezer = 2*pi*c/tweezer_wavelength
    
    U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist)
    w_tweezer_r =  omega_tweezer_r(U,beam_waist,m)
    w_tweezer_a = omega_tweezer_a(U,beam_waist,tweezer_wavelength,m)
    return individual_freqs_to_mode_vectors_axial(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a)

def physical_params_to_radial_mode_vectors_weak(N,tweezed_ions,tweezer_wavelength,linewidths,omega_res,w_rf_a,w_rf_r,P_opt,beam_waist,m):
    """
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
    
    outputs:
    modes from mode_calc_r, frequencies in Hz (not angular)
    
    """
   
    omega_tweezer = 2*pi*c/tweezer_wavelength
    
    U = potential(omega_tweezer,linewidths,omega_res,P_opt,beam_waist)
    w_tweezer_r =  omega_tweezer_r(U,beam_waist,m)
    w_tweezer_a = omega_tweezer_a(U,beam_waist,tweezer_wavelength,m)
    return individual_freqs_to_mode_vectors_radial_weak(N,tweezed_ions,w_tweezer_r,w_tweezer_a,w_rf_r,w_rf_a)