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
from tether_analysis.TetherDesign import RoundTetherDesign
from databases.DatabaseControl import build_multiplier

def average_margined_spool_capacity(innerDiameter, outerDiameter, spoolWidth, tether, packMargin=0):
    """Calculates the average capacity of a spool with a desired margin

    Args:
        innerDiameter (float): Inner diameter of the spool (mm)
        outerDiameter (float): Outer diameter of the spool (mm)
        spoolWidth (float): Width of the spool
        tether (RoundTetherDesign): Tether design object
        packMargin (int, optional): Pack margin. 0.3 adds a 30% knockdown on the capacity Defaults to 0.

    Raises:
        ValueError: A round tether design object must be passed
        ValueError: The pack margin cannot be negative

    Returns:
        float: The averaged margined capacity of the described spool in meters. 
    """

    if not issubclass(type(tether), RoundTetherDesign):
        raise ValueError("An object of type %r was passed as the tether design! Must be a subclass of RoundTetherDesign"
                          % type(tether))
    if packMargin < 0:
        raise ValueError("Packing margin cannot be less than 0%!")
    
    # Calculate the average of our helix pack and square pack (we expect to be somewhere in between)
    square_capacity = squarepack_spool_capacity(innerDiameter, outerDiameter, spoolWidth, tether)
    hex_capacity = hexpack_spool_capacity(innerDiameter, outerDiameter, spoolWidth, tether)
    average_capacity = (square_capacity + hex_capacity) / 2

    # Add a margin specified by the user to make sure we have enough length even though other factors can't be modeled #
    margined_capacity = average_capacity * (1 - packMargin)

    return margined_capacity


def squarepack_spool_capacity(innerDiameter, outerDiameter, spoolWidth, tether):
    """Approximates tether capacity on a spool by counting the length of cocentric tori along the spools width

    Args:
        innerDiameter (float): Inner diameter of the spool (mm)
        outerDiameter (float): Outer diameter of the spool (mm)
        spoolWidth (float): Width of the spool (mm)
        tether (RoundTetherDesign): Tether design object 

    Raises:
        ValueError: A round tether design object must be passed
        ValueError: Outer diameter must be >= inner diameter

    Returns:
        float: The capacity of the described spool in meters.
    """
    
    if(innerDiameter <= 0 or outerDiameter <= 0 or spoolWidth <= 0):
        raise ValueError("One or more spool dimensions are less than or equal to zero!")
    if(outerDiameter < innerDiameter):
        raise ValueError("Outer diameter is smaller than inner diameter!")

    # Grab the starting radius of our winding, and our tether radius #
    start_r = innerDiameter / 2

    # Length #
    length = 0

    # Number of cocentric Tori that can fit in radius #
    diam_cnt = int(np.floor((outerDiameter/2 - innerDiameter/2) / tether.diameter))

    # Number of tori that can fit in width #
    width_cnt = spoolWidth / tether.diameter

    # Stack cocentric toruses and multiply by the number that can fit in the width
    for i in range(0, diam_cnt):
        circumference = (tether.radius + start_r + (i * tether.diameter)) * 2 * np.pi
        length += circumference * width_cnt

    mult = build_multiplier("mm", "m")

    # Return length converted from mm to m #
    return length * mult


def hexpack_spool_capacity(innerDiameter, outerDiameter, spoolWidth, tether):
    """Approximates the tether capacity on a spool by counting the length of cocentric tori, 
    fitting subsequent layers into the gaps of the previous

    Args:
        innerDiameter (float): Inner diameter of the spool (mm)
        outerDiameter (float): Outer diameter of the spool (mm)
        spoolWidth (float): Width of the spool (mm)
        tether (RoundTetherDesign): Tether design object 

    Raises:
        ValueError: A round tether design object must be passed
        ValueError: Outer diameter must be >= inner diameter

    Returns:
        float: The capacity of the described spool in meters.
    """
    
    if(innerDiameter <= 0 or outerDiameter <= 0 or spoolWidth <= 0):
        raise ValueError("One or more spool dimensions are less than or equal to zero!")
    if(outerDiameter < innerDiameter):
        raise ValueError("Outer diameter is smaller than inner diameter!")

    # Grab the starting radius of our winding, and our tether radius #
    start_r = innerDiameter / 2

    # Length #
    length = 0

    # Calculate our height increment for each layer
    height_inc = np.sqrt(np.power(tether.diameter, 2) - np.power(tether.radius, 2))

    # Calculate the # of wraps for the OD and width #
    diam_cnt = int(np.floor((outerDiameter/2 - innerDiameter/2) / height_inc))
    width_cnt = spoolWidth / tether.diameter

    # Add the length of all of our tori, alternating between n and n-1 wraps for each layer #
    for i in range(0, diam_cnt):
        circumference = (tether.radius + start_r + (i * height_inc)) * 2 * np.pi
        w_mult = width_cnt
        if i % 2 == 1: w_mult -= 1
        length += circumference * w_mult

    mult = build_multiplier("mm", "m")

    # Return length converted from mm to m #
    return length * mult
        


def determine_spool_width(tether, innerDiameter, outerDiameter, packMargin=0, verbose=True):
    """Iteratively determines the width of a spool to fit the desired tether. Outer diameter and
    inner diameter are set by the user. 

    Args:
        tether (RoundTetherDesign): Round tether design object
        innerDiameter (float): Inner diameter of the spool (mm)
        outerDiameter (float): Outer diameter of the spool (mm)
        packMargin (int, optional): Additional margin to increase for, 0.3 is 30%. Defaults to 0.
        verbose (bool, optional): Whether to print spool dimension information. Defaults to True.

    Returns:
        float: Width of the spool in mm
    """
    
    # Start with a spool 1 tether wide
    spoolWidth = tether.diameter

    # Calculate margined spool width #
    margined_spool_cap = 0
    while margined_spool_cap < tether.length:
        spoolWidth += 0.1
        margined_spool_cap = average_margined_spool_capacity(innerDiameter, outerDiameter, spoolWidth, tether, 
                                                             packMargin=packMargin)

    if(verbose):
        print("Necessary spool width to accommodate tether with diameter %f mm and length %f m: " 
              % (tether.diameter, tether.length))
        print("  - Inner Diameter: %f mm" % innerDiameter)
        print("  - Outer Diameter: %f mm" % outerDiameter)
        print("  - Spool Width:    %f mm" % spoolWidth)
        print("  - Capacity:       %f m"  % margined_spool_cap)

    return spoolWidth


def determine_spool_od(tether, innerDiameter, width, packMargin=0, verbose=True):
    """Iteratively determines the outer diameter of a spool to fit the desired tether. Inner diameter
    and width are set by the user. 

    Args:
        tether (RoundTetherDesign): Round tether design object
        innerDiameter (float): Inner diameter of the spool (mm)
        width (float): Width of the spool (mm)
        packMargin (int, optional): Additional margin to increase for, 0.3 is 30%. Defaults to 0.
        verbose (bool, optional): Whether to print spool dimension information. Defaults to True.

    Returns:
        float: Outer diameter of the spool in mm
    """

    outerDiameter = innerDiameter

    # Calculate the spool OD margined capacity #
    margined_spool_cap = 0
    while margined_spool_cap < (tether.length):
        margined_spool_cap = average_margined_spool_capacity(innerDiameter, outerDiameter, width, tether, 
                                                             packMargin=packMargin)
        outerDiameter += 0.1

    if(verbose):
        print("Necessary spool outer diameter to accommodate tether with diameter %f mm and length %f m: " 
              % (tether.diameter, tether.length))
        print("  - Inner Diameter: %f mm" % innerDiameter)
        print("  - Outer Diameter: %f mm" % outerDiameter)
        print("  - Spool Width:    %f mm" % width)
        print("  - Capacity:       %f m"  % margined_spool_cap)

    return outerDiameter


def determine_spool_vol(outerDiameter, innerDiameter, spoolWidth):
    """Helper function to calculate the volume of a spool available for a tether. 
    Spool is described by the passed dimensions.

    Args:
        outerDiameter (float): Outer diameter of the spool in mm
        innerDiameter (float): Inner diameter of the spool in mm
        spoolWidth (float): Width of the spool in mm

    Returns:
        float: Volume of the spool in mm^3
    """
    return ((np.pi * np.power(outerDiameter/2, 2)) - (np.pi * np.power(innerDiameter/2, 2)))* spoolWidth 


