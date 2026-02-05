""""
Copyright 2024, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government 
Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the 
California Institute of Technology. This software may be subject to U.S. export control laws. 
By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. 
User has the responsibility to obtain export licenses, or other export authority as may be required before exporting 
such information to foreign countries or providing access to foreign persons."

If you have questions regarding this, please contact the JPL Software Release Authority at x4-2458.
"""

import databases.DatabaseControl as DB
import numpy as np

databases = DB.DatabaseControl()


def calculateMinBendRadius(tether):
    """Calculates the minimum bend radius of the tether in mm.
    Based on rules of thumb from:
     https://www.thefoa.org/tech/ref/install/bend_radius.html#:~:text=The%20normal%20recommendation%20for%20fiber,10%20times%20the%20cable%20diameter
     https://www.anixter.com/en_us/resources/literature/wire-wisdom/minimum-bend-radius.html

    """

    tetherRule = 8 * tether.diameter

    elecPaths = tether.findWires("electrical")
    fiberPaths = tether.findWires("optical")

    maxWireOD = 0
    maxFiberOD = 0

    for path in elecPaths:
        wire = tether.getLayerAtPath(path)
        if wire.outerDiameter > maxWireOD:
            maxWireOD = wire.outerDiameter

    for path in fiberPaths:
        fiber = tether.getLayerAtPath(path)
        if fiber.outerDiameter > maxFiberOD:
            maxFiberOD = fiber.outerDiameter

    wireRule = 12 * maxWireOD
    fiberRule = 20 * maxFiberOD # TODO add Curtis's math (ask for source)

    return max(tetherRule, wireRule, fiberRule)

def straightRunStrength(tether, output=True, verbose=False):
    """Calculates the strength of the tether in newtons, assuming
    that all materials see load immediately. 

    Returns:
        float: The strength in newtons. 
    """
    
    tetherLayers = []
    tether.layer.countLayers(list=tetherLayers)

    min_strain = np.inf
    eff_spring_const = 0

    layerDict = {}

    for layer in tetherLayers:
        
        # Grab yield stress and young's modulus #
        # material_entry = databases.get_material_entry(layer.layerMaterial)
        material_entry = layer.material_entry
        yield_stress = DB.get_material_property(material_entry, "stress")
        young_mod = DB.get_material_property(material_entry, "elastic_modulus")

        if young_mod == 0:
            continue
        
        # Calculate strain, and update min strain #
        strain = yield_stress / young_mod
        if strain < min_strain:
            min_strain = strain

        # Calculate the cross sectional area of the layer #
        mm_to_m = DB.build_multiplier("mm^2", "m^2")
        csa = (layer.outerRadius**2 - layer.innerRadius**2) * np.pi * mm_to_m * layer.fillRatio * np.sin(np.radians(layer.helixAngle))

        # Calculate the spring constant for the layer (young's * cross section) #
        gpa_to_pa = 1 / DB.build_multiplier("Pa", "GPa")
        spring_const = young_mod * csa * gpa_to_pa
        eff_spring_const += spring_const

        layerDict[layer.layerPath] = [layer.name, layer.layerPath, 0, 0, spring_const]

    break_force = min_strain * eff_spring_const

    for key, val in layerDict.items():
        layerDict[key][2] = layerDict[key][4] * min_strain
        layerDict[key][3] = min_strain

    if(output):
        print("--- Tether Straight-Run Strength Analysis ---")
        print("  Calculated Strength: {} N".format(break_force))
        print("  Tether Elongation: {} m".format(min_strain * tether.length))

        if(verbose):
            print("  Per-Layer Results:")
            for k in layerDict.keys():
                properties = layerDict[k]
                if properties[2] == 0:
                    continue
                print("  Name: {}, Path: {}".format(properties[0], properties[1]))
                print("    - Strength Contribution: {} N".format(properties[2]))
                print("    - Stretch: {} m".format(properties[3]))
                print("    - Spring Constant: {}".format(properties[4]))

    return break_force


def unhelixStrength(tether, output=True, verbose=False):
    
    """Calculates the strength of the tether, allowing helixed members to
    straighten. (Not all members see load immeidately)

    Args:
        output (bool, optional): Whether to print analysis results. Defaults to True.
        verbose (bool, optional): Whether to include per-layer results. Defaults to False.

    Returns:
        _type_: _description_
    """

    tetherLayers = []
    tether.layer.countLayers(list=tetherLayers)

    # Sort layers by ascending order of length #    
    tetherLayers.sort(key=lambda x: x.length)

    # the length the entire tether needs to stretch to in order to hit the yield strain of this member #
    yield_list = [] 
    layerDict = {}
    
    tether_break_len = 1e9
    if(output):
        print("--- Tether Unhelix-Allowed Strength Analysis ---")
    
    for layer in tetherLayers:

        layerDict[layer.layerPath] = [layer.name, layer.layerPath]

        # Grab yield stress and young's modulus #
        material_entry = databases.get_material_entry(layer.layerMaterial)
        yield_stress = DB.get_material_property(material_entry, "stress")
        young_mod = DB.get_material_property(material_entry, "elastic_modulus")

        if young_mod == 0:
            continue
        
        # Calculate effective tether length, and update min strain #
        # layer_length = layer._calcFillLength()
        layer_length = layer.length
        strain = yield_stress / young_mod
        yield_length = layer_length * (1 + strain)
        yield_list.append({"layer": layer.name, "yield_length": yield_length})
        
        if yield_length < tether_break_len:
          tether_break_len = yield_length
        
        if verbose:
          print("  Name: {}, Path: {}".format(layer.name, layer.layerPath))
          print("    - Yield Strain: {}".format(strain))
          print("    - Rest Length: {}".format(layer_length))
          print("    - Yield Length: {}".format(yield_length))
          
        
    
    # Find the breaking length of the tether (first layer to hit yield strain w.r.t. tether length) # 
    for break_length in yield_list:
      # print(break_length)
      if break_length['yield_length'] == tether_break_len:
        print("Weak Link: \"%s\" at %.3d m" % (break_length['layer'], break_length['yield_length']))
        break
    
    # tether_break_len = min(yield_length_tether)

    break_force = 0

    forceCarryingLayers = []


    # Build a list of layers that will see stretch & contribute to strength #
    for layer in tetherLayers:
        layer_length = layer._calcFillLength()
        layer_length = layer.length
        if layer_length < tether_break_len:

            # Calculate force contribution #
            material_entry = databases.get_material_entry(layer.layerMaterial)
            yield_stress = DB.get_material_property(material_entry, "stress")
            young_mod = DB.get_material_property(material_entry, "elastic_modulus")

            if young_mod == 0:
                layerDict[layer.layerPath].extend([0, 0, 0])
                continue

            forceCarryingLayers.append(layer)

            strain = yield_stress / young_mod
            mm_to_m = DB.build_multiplier("mm^2", "m^2")
            csa = (layer.outerRadius**2 - layer.innerRadius**2) * np.pi * mm_to_m * layer.fillRatio
            gpa_to_pa = 1 / DB.build_multiplier("Pa", "GPa")
            spring_const = (young_mod * csa * gpa_to_pa) / layer_length * np.sin(np.radians(layer.helixAngle))

            layerDict[layer.layerPath].extend([spring_const * (tether_break_len - layer_length), (tether_break_len - layer_length), spring_const])

            break_force += spring_const * (tether_break_len - layer_length)

        else:
            layerDict[layer.layerPath].extend([0, 0, 0])

    if(output):
        print("--- Tether Unhelix-Allowed Strength Analysis ---")
        print("  Calculated Strength: {} N".format(break_force))
        print("  Tether Elongation: {} m".format(tether_break_len - tether.length))
        
        if verbose:
            print("  Per-Layer Results:")
            for k in layerDict.keys():
                properties = layerDict[k]
                if properties[2] == 0:
                    continue
                print("  Name: {}, Path: {}".format(properties[0], properties[1]))
                print("    - Strength Contribution: {} N".format(properties[2]))
                print("    - Stretch: {} m".format(properties[3]))
                print("    - Spring Constant: {}".format(properties[4]))

    return break_force
  
def calcHelixCurvature(lead, r):
    import math
    R = r * (1 + ( lead / (2*math.pi*r))**2 )
    return R

def calcMemberLead(layer):
    cir = layer.polarCoords[0] * 2 * np.pi
    memberLead = cir * np.tan(np.deg2rad(layer.parent.helixAngle))
    return memberLead