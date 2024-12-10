""""
Copyright 2024, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government 
Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the 
California Institute of Technology. This software may be subject to U.S. export control laws. 
By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. 
User has the responsibility to obtain export licenses, or other export authority as may be required before exporting 
such information to foreign countries or providing access to foreign persons."

If you have questions regarding this, please contact the JPL Software Release Authority at x4-2458.
"""


# Sweep Modulation frequency (F_c) through a defined range (determines attenuation)
#   Plot F_s & attenuation 

import numpy as np
from scipy import special
from scipy.optimize import root_scalar
import mpmath as mp
import databases.DatabaseControl as DB
import calculation_libraries.PowerAnalysis as PA
import tether_analysis.TetherDesign as TD

import warnings

databases = DB.DatabaseControl()

### Sources used for calculations: ###
### Book 1: Transmission Lines and Lumped Circuits, Giovanni Miano & Antonio Maffucci ###
### Book 2: Analysis of Multiconductor Transmission Lines, Clayton R. Paul ###
### Zint: Practical continuous functions for the internal impedance of solid cylindrical conductors ###

### Constants for Permeability and Permittivity of free space ###
MU_0        =   4 * np.pi * 1e-7 # (H/m)
EPSILON_0   =   8.85418782e-12   # (m^-3 kg^-1 s^4 A^2)


def determine_attenuation(tether, send_path, return_path, frequency, temp):
    """Calculates the attenuation in dB for two comm wires in a tether, at a given frequency and temperature. 

    Args:
        tether (RoundTetherDesign): Tether design object
        send_path (string): Path within the tether design to find the send wire
        return_path (string): Path within the tether design to find the return wire
        frequency (int): Frequency of the comm system (Hz)
        temp (float): Temperature of the conductors (C)

    Returns:
        float: Attenuation (dB)
    """

    #TODO: Remove this once development is finished
    warnings.warn("*** THIS FUNCTIONALITY IS UNDER DEVELOPMENT AND IS STILL UNVERIFIED! Make sure to verify results! ***", )

    ### Grab the send and return wire objects ###s
    send_wire = tether.getLayerAtPath(send_path)
    return_wire = tether.getLayerAtPath(return_path)

    ### Convert to angular frequency ###
    omega = 2 * np.pi * frequency

    ### Calculate RLCG parameters from tether geometry, materials and frequency (L, C, & G are aymptotic approx.) ###
    R = resistance_parameter(send_wire, return_wire, frequency, temp)
    L = inductance_parameter(send_wire, return_wire)
    C = capacitance_parameter(send_wire, return_wire)
    G = conductance_parameter(send_wire, return_wire)

    ### Calculate the attenuation coefficient ###
    rlcg = mp.sqrt((R + 1j * omega * L )*(G + 1j * omega * C))
    atten_coeff = mp.re(rlcg)

    ### Convert attenuation from Np/m to dB/m ###
    atten_dB = atten_coeff * (20 / np.log(10))

    return atten_dB * send_wire.length

def characteristic_impedance(tether, send_path, return_path, frequency, temp):
    """Calculates the characteristic impedance for two comm wires in a tether, at a given frequency and temperature. 

    Args:
        tether (RoundTetherDesign): Tether design object
        send_path (string): Path within the tether design to find the send wire
        return_path (string): Path within the tether design to find the return wire
        frequency (int): Frequency of the comm system (Hz)
        temp (float): Temperature of the conductors (C)

    Returns:
        float: Characteristic impedance NOTE: figure out units
    """
    
    #TODO: Remove this once development is finished
    warnings.warn("*** THIS FUNCTIONALITY IS UNDER DEVELOPMENT AND IS STILL UNVERIFIED! Make sure to verify results! ***", )

    ### Grab the send and return wire objects ###s
    send_wire = tether.getLayerAtPath(send_path)
    return_wire = tether.getLayerAtPath(return_path)

    ### Convert to angular frequency ###
    omega = 2 * np.pi * frequency

    ### Calculate RLCG parameters from tether geometry, materials and frequency (L, C, & G are asymptotic approx.) ###
    R = resistance_parameter(send_wire, return_wire, frequency, temp)
    L = inductance_parameter(send_wire, return_wire)
    C = capacitance_parameter(send_wire, return_wire)
    G = conductance_parameter(send_wire, return_wire)

    ### Calculate the attenuation coefficient ###
    imped = mp.sqrt((R + 1j * omega * L )/(G + 1j * omega * C))

    return imped



def watts_to_dBm(level):
    """Converts power in Watts to dBm

    Args:
        level (float): Power (W)

    Returns:
        float: Power(dBm)
    """
    return 10 * np.log10(float(1000 * level))

def dBm_to_watts(level):
    """Converts power in dBm to Watts

    Args:
        level (float): Power in dBm

    Returns:
        float: Power (W)
    """
    return 10**(level/10) / 1000


def dB_to_ratio(level):
    """Converts ratio between two values in dB to unitless

    Args:
        level (float): Ratio between two values expressed as dB

    Returns:
        float: Unitless ratio between two values
    """
    return 10 ** (level/10)


def snr(P, f_s, P_noise):
    """Calculates Signal to Noise Ratio (SNR) per bit from the received power, switching frequency, and noise power.
    OOK specific

    Args:
        P (float): Received power level (dBm)
        f_s (float): Switching frequency (Hz)
        P_noise (float): Noise power level (dBm/Hz)

    Returns:
        float: SNR 
    """

    P_tx = dBm_to_watts(P)
    noise = dBm_to_watts(P_noise)
    E_s = P_tx/f_s
    return(E_s / noise)

    # In general Ambient temp is -170 dBm/Hz 


def ber(snr):
    """Calculates the bit error rate using the Q function

    Args:
        snr (float): Signal to Noise Ratio (SNR)

    Returns:
        float: Probability of a bit error
    """

    ### Source: https://gzipwtf.com/what-is-the-relationship-between-snr-and-ber/, video at the end, time: 13:24 ###
    # (Given as a source by Dr. Aurenice Oliveira, Michigan Technological University)
    return Q(mp.sqrt(snr))


def Q(x):
    """Q function for a normal distribution

    Args:
        x (float): Input 

    Returns:
        float: Corresponding probability of the tail of a normal distribution
    """

    ### Source: https://gzipwtf.com/what-is-the-relationship-between-snr-and-ber/, video at the end, time: 13:24 ###
    # (Given as a source by Dr. Aurenice Oliveira, Michigan Technological University)
    return 0.5*special.erfc(float(x/np.sqrt(2)))


def calc_ook_datarate(P_dB, N_0, max_ber=1e-9):
    """Calculates the datarate for OOK

    Args:
        P_dB (float): Power level at the receiver (dBm)
        N_0 (float): Noise power level (dBm)
        max_ber (float, optional): Max acceptable BER. Defaults to 1e-9.
        negative_low (bool, optional): Whether the logical low value is just negative high, or zero. Defaults to True.

    Returns:
        float: Data rate for the channel (switching frequency) (bps)
    """

    #TODO: Remove this once development is finished
    warnings.warn("*** THIS FUNCTIONALITY IS UNDER DEVELOPMENT AND IS STILL UNVERIFIED! Make sure to verify results! ***", )

    ### Convert to watts ###
    P_rx = dBm_to_watts(P_dB)
    P_noise = dBm_to_watts(N_0)

    ### Calculate switching frequency (capacity) that meets the max BER ###
    term1 = np.sqrt(2) * special.erfcinv(2*max_ber)
    return (P_rx / (P_noise * (term1 ** 2)))


def ook_datarate_analysis(tether, send_path, return_path, transmit_power, noise_power, max_ber=1e-9, temp=20):
    """Estimates the datarate of the paths in a given tether design. User enters the send/return path, transmit power,
    noise power, temperature and the max allowable BER. 

    Args:
        tether (RoundTetherDesign): Round tether design object
        send_path (string):Tether path of the send wire
        return_path (string): Tether path of the return wire
        transmit_power (float): Transmission power (dBm)
        noise_power (float): Noise power (dBm)
        max_ber (float, optional): The max allowed bit error rate. Defaults to 1e-6.
        temp (int, optional): The temperature of the wire (C). Defaults to 20.

    Raises:
        RuntimeError: A solution could not be found by the root finding

    Returns:
        float: The estimated max data rate at that BER (Mbps)
    """

    # Method makes sense, make sure to check the carrier frequency is larger than the switching frequency 
    # - Doesn't meeting the Nyquist rate 

    #TODO: Remove this once development is finished
    warnings.warn("*** THIS FUNCTIONALITY IS UNDER DEVELOPMENT AND IS STILL UNVERIFIED! Make sure to verify results! ***", )

    ### Specicy the tuple of arguments to pass to the objective function and calculated datarate upper limit ###
    argTup = (send_path, return_path, transmit_power, noise_power, tether, temp, max_ber)

    upper_limit = calc_ook_datarate(argTup[2], argTup[3], max_ber=argTup[6]) 

    ### Use bisect method to find the root ###
    solution = root_scalar(ook_objective_fucn, args=(argTup), method='bisect', x0=0, bracket=[upper_limit, 1])

    ### Return the solution if we found one, anotherwise throw an error ###
    if(solution.converged):
        return solution.root
    else:
        print(solution)
        raise RuntimeError("Was unable to find a solution to the datarate analysis!")


# only need send path
def ook_objective_fucn(x, *args):
    """Objective function for On/Off Keying Datarate analysis, used by root-finding method

    Args:
        x (float): Switching frequency (Hz)
        args: Tuple of arguments that should have the following format:
        - args[0]: send_path (string)
        - args[1]: return_path (string)
        - args[2]: transmit_power (float, dBm)
        - args[3]: noise_power (float, dBm)
        - args[4]: tether (RoundTetherDesign)
        - args[5]: temp (float, C)
        - args[6]: max_ber (float) 

    Returns:
        float: Difference between calculated and desired BER
    """
    attendB = determine_attenuation(args[4], args[0], args[1], x, args[5])

    P_rx = args[2] - attendB
    SNR = snr(P_rx, x, args[3])
    # print("Frequency: %f Hz, BER: %.21f, SNR: %f, Max BER: %.12f, difference: %.21f" % (x, ber(SNR), SNR, args[6], ber(SNR) - args[6]))

    return ber(SNR) - args[6]


def determine_characteristic_impedance():
    pass


def ook_datarate_estimate(tether, sendWire, returnWire):
    pass


def qam_datarate_estimate():
    pass


def wire_separation_distance(wire1, wire2):
    """Calculates the separation distance of two wires in mm

    Args:
        wire1 (Wire): Send path wire object
        wire2 (Wire): Receive path wire object

    Returns:
        float: Euclidean distance between the centers of each wire
    """
    return np.sqrt(np.abs((wire1.x - wire2.x)**2) + np.abs((wire1.y - wire2.y)**2))


def permittivity(k):
    """Calculates the magnetic permittivity from the relative permittivity

    Args:
        k (float): Relative permittivity

    Returns:
        float: Magnetic permittivity
    """
    return k * EPSILON_0

def permeability(mu):
    """Calculates the magnetic permeability from the relative permeability

    Args:
        mu (float): Relative permeability

    Returns:
        float: Magnetic permeability
    """
    return mu * MU_0

def conductance(rho):
    """Returns the conductance from the resistivity, for readability

    Args:
        rho (float): Resistivity

    Returns:
        float: conductance
    """
    return 1/rho

def round_wire_constant(wire_r):
    """Returns a constant used for round wires. A portion of equation 5.41 in 
    Book 1

    Args:
        wire_r (float): Radius of wire in millimneters

    Returns:
        float: Constant used for round wires (m^-1)
    """
    return 1 / (np.pi * (wire_r/1000))

def wire_dist_ratio(wire1, wire2):
    """Calculate the ratio of the distance between the wires to their size.
    Assumes both wires are the same size. 

    Args:
        wire1 (Wire): First wire 
        wire2 (Wire): Second wire

    Returns:
        float: Ratio of the distance of the wires to 2 * their radius
    """
    sep_dist = wire_separation_distance(wire1, wire2)
    return sep_dist / (2*wire1.innerLayer.outerRadius)


def wire_checks(send_wire, return_wire):
    """Enforces constraints on wires for use in comms analysis. 
    
    Wires must:
     - Have the same insulation/conductor material'
     - Have the same dimensions

    Args:
        send_wire (Wire): Send wire
        return_wire (Wire): Return wire

    Raises:
        ValueError: Wires have different insulation material
        ValueError: Wires have different insulation thicknesses
        ValueError: Wires have different conductor material
        ValueError: Wires have different conductor sizes
    """

    send_conductor = send_wire.innerLayer
    return_conductor = return_wire.innerLayer
    if(send_wire is return_wire):
        raise ValueError("Both wires are the same! Does not refer to two conductor paths!")
    if(send_wire.layerMaterial != return_wire.layerMaterial):
        raise ValueError("Both wires must have the same insulation material!")
    if(send_wire.layerThickness != return_wire.layerThickness):
        raise ValueError("Both wires must have the same thickness of insulation!")
    if(send_conductor.layerMaterial != return_conductor.layerMaterial):
        raise ValueError("Both wires must have the same conductor material!")
    if(send_conductor.layerThickness != return_conductor.layerThickness):
        raise ValueError("Both wires must have the same size conductors!")


def resistance_parameter(send_wire, return_wire, frequency, temperature=20):
    """Calculates the resistance parameter for the RLCG model

    Args:
        tether (RoundTetherDesign): Round tether design object for the design, used for the tether length
        send_wire (Wire): Send wire object for the comm line
        return_wire (Wire): Return Wire object for the comm line
        frequency (float): Frequency in Hz
        temperature (int, optional): Temperature of the conductor. Defaults to 20.

    Raises:
        ValueError: Frequency must be greater than 0

    Returns:
        float: Reistance parameter for the RLCG model in Ohm/m
    """

    ### Raise errors if wires differ ###
    wire_checks(send_wire, return_wire)

    ### Check frequency, in this case it must be nonzero for inductance ###
    if(frequency <= 0):
        raise ValueError("Frequency must be greater than zero for impedance calcs!")
    
    ac_res = PA.ac_resistance(send_wire, frequency, temp=temperature) + PA.ac_resistance(return_wire, frequency, temp=temperature)

    return ac_res / send_wire.length


    ### Grab materials for conductor ###
    conductor_mat = send_wire.innerLayer.layerMaterial

    ### Calculate the resistivity of the wire ###
    unit = databases.materialDefaultUnits["resistivity"]
    res_mult = DB.build_multiplier("ohm*m", unit)
    rho = DB.get_material_property(databases.get_material_entry(conductor_mat), "resistivity") * res_mult

    ### Calculate the magnetic permeability of the wire ###
    relative_mu = DB.get_material_property(databases.get_material_entry(conductor_mat), "relative_permeability")
    mu_c = permeability(relative_mu)    

    ### Calculate some parameters we're going to need ###
    omega = 2 * np.pi * frequency
    W = round_wire_constant(send_wire.innerLayer.outerRadius)
    depth = PA.skin_depth(rho, frequency, mu_c)

    ### Calculate the resistive portion of the impedance from Zint (Ohm/m) ###
    Rac = PA.ac_resistance(send_wire, 0, temperature) / send_wire.length

    ### Calculate the surface inductance of the wire due to skin effect, Equation 5.44 in Book 1, (H/m), Page 197 ###
    L_s = W * (MU_0 * depth) / 4

    ### Calculate the surface impedance of the wire due to skin effect, Equation 5.42 in Book 1, (Ohm/m), Page 196  ###
    Z_s =  Rac + 1j * omega * L_s
    
    ### The 2 comes from the first part of equation 5.38 in Book 1, (Ohm/m) on page 196 ###
    return 2*Z_s 
    


def inductance_parameter(send_wire, return_wire):
    """Calculates an asymptotic approximation of the inter-wire inductance in the tether. 
    
    Args:
        send_wire (Wire): Send path for the comm line
        return_wire (Wire): Return path for the comm line

    Returns:
        float: Per-unit-length inductance (H/m)
    """

    ### Raise errors if wires differ ###
    wire_checks(send_wire, return_wire)

    ### Grab materials for conductor ###
    conductor_mat = send_wire.innerLayer.layerMaterial

    ### Calculate the magnetic permeability of the wire ###
    relative_mu = DB.get_material_property(databases.get_material_entry(conductor_mat), "relative_permeability")
    mu_c = permeability(relative_mu)

    ### Calculate asymptotic external inductance of the wire, equation 4.41 in Book 2, Page 126 ###
    L = (mu_c / np.pi) * mp.acosh(wire_dist_ratio(send_wire, return_wire))

    return L

def capacitance_parameter(send_wire, return_wire):
    """Calculates an asymptotic approximation of the inter-wire capacitance in the tether.

    Assumes all material separating wire is their same insulation material

    Args:
        send_wire (Wire): Send path for the comm line
        return_wire (Wire): Return path for the comm line

    Returns:
        float: Per-unit-length capacitance, (F/m)
    """
    
    ### Raise errors if wires differ ###
    wire_checks(send_wire, return_wire)

    ### Grab materials for insulation ###
    insulation_mat = send_wire.layerMaterial

    ### Calculate permittivity of insulation ###
    dielectric_const = DB.get_material_property(databases.get_material_entry(insulation_mat), "dielectric_constant")
    epsilon = permittivity(dielectric_const)

    ### Calculate capacitance between two wires separated by insulation. Equation 4.38 in book 2, F/m, Page 124 ###
    C = (np.pi * epsilon) / np.arccosh(wire_dist_ratio(send_wire, return_wire))

    return C


def conductance_parameter(send_wire, return_wire):
    """Calculates an asymptotic approximation of the inter-wire conductance in the the tether. 

    Assumes all material separating wire is their same insulation material

    Args:
        send_wire (Wire): Send path for the comm line
        return_wire (Wire): Return path for the comm line

    Returns:
        float: Per-unit-length Conductance, (S/m)
    """
    
    ### Raise erros if wires differ ###
    wire_checks(send_wire, return_wire)

    ### Grab materials for insulation ###
    insulation_mat = send_wire.layerMaterial

    ### Calculate permittivity of insulation ###
    unit = databases.materialDefaultUnits["resistivity"]
    res_mult = DB.build_multiplier("ohm*mm", unit)
    sigma = conductance(DB.get_material_property(databases.get_material_entry(insulation_mat), "resistivity")*res_mult)

    ### Calculate conductance between two wires separated by insulation. Equation 4.42 in book 2, S/m, Page 126 ###
    G = (np.pi * sigma) / np.arccosh(wire_dist_ratio(send_wire, return_wire))

    return G


def errorMQAM01():
    pass

def etaMQAM00():
    pass
