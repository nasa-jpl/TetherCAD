""""
Copyright 2024, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government 
Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the 
California Institute of Technology. This software may be subject to U.S. export control laws. 
By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. 
User has the responsibility to obtain export licenses, or other export authority as may be required before exporting 
such information to foreign countries or providing access to foreign persons."

If you have questions regarding this, please contact the JPL Software Release Authority at x4-2458.
"""

import numpy as np
import tether_analysis.TetherDesign as TD
import databases.DatabaseControl as DB
import calculation_libraries.CommAnalysis as CA

import warnings

databases = DB.DatabaseControl()


def calculate_wire_resistance(wire, temp):
    """Calculates the DC resistance of a passed wire/conductive layer

    Args:
        wire (Wire or Layer): The Wire or Layer object whose resistance to estimate
        temp (float): The temperature of the wire in Celsius

    Raises:
        ValueError: Wire object passed without an inner layer
        ValueError: Non-electrical Wire object was passed
        ValueError: A type other than Wire or Layer was passed
        ValueError: The conductor must be of a magnetic or ferromagnetic material type
        ValueError: The conductor (Wire.innerLayer or Layer) must have an initialized length
        ValueError: No material temperature coefficient for this material
        ValueError: Temperature less than abs zero

    Returns:
        float: The estimated DC resistance of the Wire/Layer, adjusted for temperature
    """

    # If given a wire, do some checks and set the conductor accordingly
    if issubclass(type(wire), TD.Wire):
        conductor = wire.innerLayer
        if conductor is None:
            raise ValueError("Wire must have an inner layer! This layer should be a conductive material!")
        if wire.wireType != "electrical":
            raise ValueError("A non electrical wire of type: %s was passed" % wire.wireType)
    elif issubclass(type(wire), TD.Layer):  # If given a layer, set it to the conductor
        conductor = wire
    else:
        raise ValueError("Invalid type passed! Must be of type Wire or Layer")

    # Make sure the conductor is of a suitable material
    wire_material_entry = databases.get_material_entry(conductor.layerMaterial)
    material_type = DB.get_material_property(wire_material_entry, "material_type")
    if material_type != 'metal' and material_type != 'ferromagnetic':
        raise ValueError("Only a metal or ferromagnetic metal can be passed as a wire conductor! Given type: %s" % 
                         material_type)
    
    # Make sure the conductor has an initialized length
    if conductor.length == 0:
        raise ValueError("A wire/layer with an uninitialized length was passed!")

    # Grab the resistivity for the conductor and ensure it is the unit we want #
    resistivity = DB.get_material_property(databases.get_material_entry(conductor.layerMaterial), "resistivity")
    unit = databases.materialDefaultUnits["resistivity"]
    res_mult = DB.build_multiplier("ohm*mm", unit)
    resistivity *= res_mult

    # Grab the length of the conductor and convert it to mm #
    length_mult = DB.build_multiplier("m", "mm")
    length_mm = conductor.length * length_mult

    # Calculate the cross sectional area of the conductor #
    cross_section = (conductor.outerRadius**2 - conductor.innerRadius**2) * np.pi

    resistance = resistivity * (length_mm / cross_section)

    # Adjust resistance by temperature #
    if temp != 20:
        try:
            tempCoeff = databases.tempCoeffDict[conductor.layerMaterial]
        except KeyError:
            raise ValueError("Material for temperature coeff not found in the dictionary! Run at default (20C) " + 
                            "or add the material/coefficient to the temperature coefficient dictionary!")

        tempDiff = temp - 20
        if temp < -273.15:
            raise ValueError("Temperature cannot be less than abs zero!")
        elif np.abs(tempDiff) > 100:
            warnings.warn("Accuracy of calculation lowers for temperature changes of > 100C (diff: %f)" 
                          % tempDiff, UserWarning)

    

        resistance = resistance * (1 + tempCoeff * tempDiff)

    return resistance


def calculate_equivalent_resistance(wireList, temp, ac=False, freq=0):
    """Calculates the equivalent parallel resistance of wires in the passed list.

    Args:
        wireList (list): List of Wire/Layer objects that are in parallel
        temp (float): Temperature of the wire conductor in celcius.

    Raises:
        ValueError: Passed list was empty

    Returns:
        float: Effective parallel resistance of wires in the passed list
    """

    if len(wireList) < 1:
        raise ValueError("No wires are in the passed list!")

    inverse_sum = 0
    for wire in wireList:

        if(ac == False):
            resistance = calculate_wire_resistance(wire, temp)
        else:
            resistance = ac_resistance(wire, freq, temp)

        inverse_sum += (1/resistance)

    return 1/inverse_sum


def ac_resistance(wire, frequency, temp):
    """Calculates AC resistance using the RAC - TED - ML approach specified in Zint

    Args:
        send_wire (Wire): Wire object for send path
        frequency (float): Frequency (Hz)
        temp (float): Temperature (C)

    Raises:
        ValueError: Frequency must be greater than or equal to zero
        ValueError: Temperature must be greater than or equal to abs zero

    Returns:
        float: Resistance of the wire (Ohms)
    """
        
    # If given a wire, do some checks and set the conductor accordingly
    if issubclass(type(wire), TD.Wire):
        conductor = wire.innerLayer
        if conductor is None:
            raise ValueError("Wire must have an inner layer! This layer should be a conductive material!")
        if wire.wireType != "electrical":
            raise ValueError("A non electrical wire of type: %s was passed" % wire.wireType)
    elif issubclass(type(wire), TD.Layer):  # If given a layer, set it to the conductor
        conductor = wire
    else:
        raise ValueError("Invalid type passed! Must be of type Wire or Layer")

    # Make sure the conductor is of a suitable material
    wire_material_entry = databases.get_material_entry(conductor.layerMaterial)
    material_type = DB.get_material_property(wire_material_entry, "material_type")
    if material_type != 'metal' and material_type != 'ferromagnetic':
        raise ValueError("Only a metal or ferromagnetic metal can be passed as a wire conductor! Given type: %s" % 
                         material_type)
    
    # Make sure the conductor has an initialized length
    if conductor.length == 0:
        raise ValueError("A wire/layer with an uninitialized length was passed!")

    ### Make sure frequency and temperature are within usual values ###
    if(frequency < 0):
        raise ValueError("Frequency must greater than or equal to zero!")
    if (temp < -273.15):
        raise ValueError("Temperature must be greater than or equal to absolute zero!")

    ### Calculate the DC resistance of the wire in Ohm/m ###
    R_dc = calculate_wire_resistance(wire, temp)
    
    if(frequency == 0):
        return R_dc

    ### Calculate the resistivity of the wire ###
    unit = databases.materialDefaultUnits["resistivity"]
    res_mult = DB.build_multiplier("ohm*m", unit)
    rho = DB.get_material_property(wire_material_entry, "resistivity") * res_mult

    ### Calculate the magnetic permeability of the wire ###
    relative_mu = DB.get_material_property(wire_material_entry, "relative_permeability")
    mu_c = CA.permeability(relative_mu)    

    ### Calculate the AC resistance factor for the wire from Zint, RAC - TED - ML approximation, Page 31 of Zint ###
    radius = conductor.outerRadius
    xi = ac_resistance_factor(radius, rho, mu_c, frequency)

    ### Calculate AC resistance from the factor and DC resistance, Xi, from page 10 of Zint ###
    R_ac = R_dc * xi

    return R_ac


def ac_resistance_factor(radius, rho, mu_c, frequency):
    """Esimates the AC resistance factor for given parameters. This factor scales the DC resistance
    by the appropriate amount for an AC signal.

    Args:
        radius (float): Radius of the conductor
        rho (float): Resistivity of the conductor material
        mu_c (float): Magnetic Permeability of the conductor
        frequency (float): Frequency of the signal

    Returns:
        float: AC resistance factor
    """

    if(frequency == 0):
        return 1
    
    ### Calculate the skin depth ###
    delta_i = skin_depth(rho, frequency, mu_c)    

    ### Modified skin depth for better approximation at low frequencies, Equation 9.2 from Zint, page 22 ###
    delta_i_prime = delta_i * (1 - np.exp(-(radius / delta_i)))

    ### RAC - TED - ML Approximation of Xi, from Zint, page 31 ###
    z = 0.62006 * (radius / delta_i)
    z1 = np.power(z, 1.82938)
    z2 = np.power(z, -0.99457)
    A0 = 0.272481 * np.power(z1 - z2, 2)
    A1 = np.power(1 + A0, 1.0941)
    y = 0.189774 / A1
    xi = np.power(radius, 2) / (((2 * radius * delta_i_prime) - np.power(delta_i_prime, 2)) * (1 + y))

    return xi


def skin_depth(rho, f, mu_c):
    """Calculates the skin depth from a given frequency

    Args:
        rho (float): Resistivity of conductor
        f (float): Frequency (Hz)
        mu_c (float): Magnetic Permeability of conductor
    Returns:
        double: Skin depth of the conductor at this frequency, in mm 
    """

    ### Return skin depth (delta_i) as defined on page 8 of Zint ### 
    return np.sqrt(rho / (np.pi * f * mu_c)) 


def dc_power_transmission_analysis(tether, tether_voltage, desired_power, send_paths, return_paths, print_results=True, 
                                   temp=20, verbose=False, analysis_name = ""):
    """Determines the DC power transmission characteristics of a given tether by stepping up the input power until the 
       desired output power is determined.

    Args:
        tether (RoundTetherDesign): RoundTetherDesign object used to interact with tether
        tether_voltage (float): The input voltage to the tether (V)
        desired_power (float): The desired output power from the tether (W)
        send_paths (list): A list of layer paths corresponding to wires in the tether design to be used as the send 
        paths of the transmission line
        return_paths (list): A list of layer paths corresponding to wires in the tether design to be used as the return 
        paths of the transmission line
        print_results (bool): Whether to print the results of the analysis in addition to returning the 
        efficiency of the tether
        temp (int, optional): Temperature of the conductor materials in celcius. Defaults to 20.
        verbose (bool, optional): Whether to print full send/return path information when printing is enabled. 

    Raises:
        ValueError: tether argument was not an object of type RoundTetherDesign
        ValueError: send_paths list is empty
        ValueError: return_paths list is empty
        ValueError: Voltage is <= 0
        ValueError: Desired power output is <= 0
        RuntimeError: A path is in both the send/return path lists
        ValueError: None of the send path wires were found in the tether
        ValueError: None of the return path wires were found in the tether
        ValueError: Imaginary solution returned 

    Returns:
        list: List of the efficiencies of each solution
    """

    # Throw some generic errors #
    if not issubclass(type(tether), TD.RoundTetherDesign):
        raise ValueError("Tether argument must be of type RoundTetherDesign!")
    if len(send_paths) < 1:
        raise ValueError("No send paths specified!")
    if len(return_paths) < 1:
        raise ValueError("No return paths specified!")
    if tether_voltage <= 0:
        raise ValueError("Voltage must be greater than 0!")
    if desired_power <= 0:
        raise ValueError("Desired power output must be greater than or equal to 0!")

    # Check for duplicates, note that this may need to change if coaxial designs are introduced #
    for path in send_paths:
        if path in return_paths:
            raise RuntimeError("The path %s is present in both lists! A single path cannot be both a send and return!")

    # Build lists of wires for send/return paths from names #
    send_path_wires = tether.buildLayerList(send_paths)
    return_path_wires = tether.buildLayerList(return_paths)
    if(len(send_path_wires) < 1):
        raise ValueError("Send path wires were not found in passed tether design!")
    if(len(return_path_wires) < 1):
        raise ValueError("Return path wires were not found in passed tether design!")

    # Calculate the equivalent resistance for the send and return paths #
    send_resistance_eq = calculate_equivalent_resistance(send_path_wires, temp)
    return_resistance_eq = calculate_equivalent_resistance(return_path_wires, temp)
    total_eq_resistance = send_resistance_eq + return_resistance_eq

    # Formulate the problem as a quadratic equation and solve #
    coeffs = [total_eq_resistance, -tether_voltage, desired_power]
    roots = np.roots(coeffs)
    solutions = [x for x in roots if np.isreal(x)]

    # If we don't have any real solutions throw an error # 
    if len(solutions) < 1:
        current = (min(roots))
        raise ValueError("Solution cannot be found with given parameters! Tether Voltage drop is likely too high!") 

    # Sort solutions in increasing order, the first will be the more stable #
    solutions = sorted(solutions)

    efficiencyList = []

    for idx, current in enumerate(solutions):

        voltage_drop_send = current * send_resistance_eq
        voltage_drop_return = current * return_resistance_eq
        power_loss_send = current * voltage_drop_send
        power_loss_return = current  * voltage_drop_return

        input_power = current * tether_voltage
        power_transmitted = input_power - power_loss_send - power_loss_return
        

        # Calculate our tether power efficiency #
        tether_power_efficiency = power_transmitted / input_power * 100

        # If we want per-wire data, calculate it here #
        if verbose:
            send_power_dict = {}
            return_power_dict = {}
            send_current_dict = {}
            return_current_dict = {}
            send_power_disp_dict = {}
            return_power_disp_dict = {}
            for wire in send_path_wires:
                resistance = calculate_wire_resistance(wire, temp)
                wire_current = voltage_drop_send / resistance
                wire_power_loss = wire_current * voltage_drop_send
                send_power_dict[wire.layerPath] = wire_power_loss
                send_current_dict[wire.layerPath] = wire_current
                send_power_disp_dict[wire.layerPath] = wire_power_loss / tether.length
            for wire in return_path_wires:
                resistance = calculate_wire_resistance(wire, temp)
                wire_current = voltage_drop_return / resistance
                wire_power_loss = wire_current * voltage_drop_return
                return_power_dict[wire.layerPath] = wire_power_loss
                return_current_dict[wire.layerPath] = wire_current
                return_power_disp_dict[wire.layerPath] = wire_power_loss / tether.length
                
        # Print results if desired #
        if(print_results):
            
            # Print stability of solution (as long as inputs are rational, we'll always have two or 0 solutions) # 
            stabStr = "More Stable Solution"
            if idx > 0:
                stabStr = "Less Stable Solution"

            print("--- %s Power Transmission Results for: %s ---" % (analysis_name, tether.name))
            print("  Solution Stability:            %s" % stabStr)
            print("  Send Path:")
            print("     - Eq. Resistance:           %f ohms" % send_resistance_eq)
            print("     - Voltage Drop:             %f V" % voltage_drop_send)
            print("     - Power Loss:               %f W" % power_loss_send)
            if(verbose):
                print("     - Send Paths:               %s" % str(send_paths))
                print("     - Current Per Wire (A):     %s" % str(send_current_dict))
                print("     - Total Loss Per Wire (W):  %s" % str(send_power_dict))
                print("     - Power Dissipation: (W/m): %s" % str(send_power_disp_dict))
            print("  Return Path:")
            print("     - Eq. Resistance:           %f ohms" % return_resistance_eq)
            print("     - Voltage Drop:             %f V" % voltage_drop_return)
            print("     - Power Loss:               %f W" % power_loss_return)
            if(verbose):
                print("     - Return Paths:             %s" % str(return_paths))
                print("     - Current Per Wire (A):     %s" % str(return_current_dict))
                print("     - Total Loss Per Wire (W):  %s" % str(return_power_dict))
                print("     - Power Dissipation: (W/m): %s" % str(return_power_disp_dict))
            print("  General:")
            print("     - Required Input Power:     %f W" % input_power)
            print("     - Actual Output Power:      %f W" % power_transmitted)
            print("     - Round Trip Resistance:    %f ohms" % (send_resistance_eq + return_resistance_eq))
            print("     - Tether Power Loss:        %f W" % (power_loss_send + power_loss_return))
            print("     - Tether Input Voltage:     %f V" % tether_voltage)
            print("     - Tether Output Voltage:    %f V" % (tether_voltage - voltage_drop_send - voltage_drop_return))
            print("     - Tether Current:           %f A" % current)
            print("     - Tether Power Efficiency:  %.2f %%" % tether_power_efficiency)
            print("     - Power Dissipation:        %f W/m" % ((power_loss_send + power_loss_return) / tether.length))
            print("\n")

        efficiencyList.append(tether_power_efficiency)

    #TODO: add the ability to return all of these calculated parameters
    return efficiencyList

def single_phase_ac_transmission_analysis(tether, line_to_line_voltage, desired_power, send_paths, return_paths, 
                                          print_results=True, temp=20, verbose=False, analysis_name=""):
    pass