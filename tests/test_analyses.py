""""
Copyright 2024, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government 
Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the 
California Institute of Technology. This software may be subject to U.S. export control laws. 
By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. 
User has the responsibility to obtain export licenses, or other export authority as may be required before exporting 
such information to foreign countries or providing access to foreign persons."

If you have questions regarding this, please contact the JPL Software Release Authority at x4-2458.
"""

import pytest
import pathlib
from tether_analysis.TetherDesign import *
from calculation_libraries.SpoolAnalysis import *
from calculation_libraries.PowerAnalysis import *
from calculation_libraries.CommAnalysis import *
from math import isclose
from copy import deepcopy

databases.load_from_folder(pathlib.Path.cwd() / 'test_databases')

@pytest.fixture
def simple_tether():

    layer = Layer("fep", "fep_layer", 4)

    tether = RoundTetherDesign("simple tether", layer, 100)

    return tether

@pytest.fixture
def wire_list():
    wireData = databases.wire_db_dict["Teledyne-HV-FEP-2023"]
    wireQueryDict = {"Part #" : "178-9907"}
    wireEntry = DB.searchEntry(wireData, wireQueryDict)
    wire = Wire.from_entry("electrical", wireEntry)
    wireList = []

    for i in range(0, 6):
        wireList.append(wire)

    return wireList

@pytest.fixture
def electrical_tether(wire_list):

    core = Layer("fep", "core", 0.25)
    memLayer = Layer("fep", "wires", 0.5, innerLayer=core, memberList=wire_list)
    abrasion = Layer("kapton", "abrasion", 0.15, innerLayer=memLayer)

    tether = RoundTetherDesign("simple tether", abrasion, 100)

    return tether

@pytest.fixture
def coaxial_tether():
    core = Layer("fep", "core", 0.5, color="black")
    copper1 = Layer("copper", "inner conductor", 0.25, innerLayer = core)
    ins1 = Layer("fep", "middle insulator", 0.2, innerLayer=copper1)
    copper2 = Layer("copper", "outer conductor", 0.15, innerLayer=ins1)
    outer = Layer("fep", "outer insulator", 0.1, innerLayer=copper2)

    coaxTether = RoundTetherDesign("Coaxial Tether", outer, 100)

    return coaxTether


class TestSpoolAnalysis:

    # --- Tests for squarepack spool capacity --- #
    @pytest.mark.usefixtures("simple_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_basic_squarepack(self, simple_tether):

        cap1 = squarepack_spool_capacity(40, 48, 8, simple_tether)
        assert cap1 == 0

        cap2 = squarepack_spool_capacity(311, 327, 8, simple_tether)
        cap3 = squarepack_spool_capacity(311, 327, 32, simple_tether)

        assert isclose(cap2 * 4, cap3, rel_tol=1e-9)
    
    # --- Tests for helixpack spool capacity --- #
    @pytest.mark.usefixtures("simple_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_basic_helixpack(self, simple_tether):

        assert squarepack_spool_capacity(311, 327, 8, simple_tether) == hexpack_spool_capacity(311, 327, 8, 
                                                                                                 simple_tether)
        assert squarepack_spool_capacity(311, 327, 32, simple_tether) == hexpack_spool_capacity(311, 327, 32, 
                                                                                                  simple_tether)

        assert squarepack_spool_capacity(40, 48, 8, simple_tether)== 0

    # --- Tests for average margined spool capacity --- #
    @pytest.mark.usefixtures("simple_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_margin_capacity(self, simple_tether):

        margin = 0.33
        multiplier = 1-margin

        no_margin_cap = average_margined_spool_capacity(300, 330, 100, simple_tether)
        margin_cap = average_margined_spool_capacity(300, 330, 100, simple_tether, packMargin=margin)

        assert (no_margin_cap * multiplier) == margin_cap


    @pytest.mark.usefixtures("simple_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_determine_spool_width(self, simple_tether):
        
        width1 = determine_spool_width(simple_tether, 300, 330, verbose=False)
        width2 = determine_spool_width(simple_tether, 300, 330, packMargin=0.33, verbose=False)

        assert width1 < width2

        # Test some edge cases
        with pytest.raises(ValueError):
            determine_spool_width(simple_tether, 0, 330, verbose=False)
        with pytest.raises(ValueError):
            determine_spool_width(simple_tether, -10, 330, verbose=False)
        with pytest.raises(ValueError):
            determine_spool_width(simple_tether, 360, 330, verbose=False)

    # --- Tests for determine_spool_width --- #
    @pytest.mark.usefixtures("simple_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_determine_spool_od(self, simple_tether):   
        width1 = determine_spool_od(simple_tether, 300, 330, verbose=False)
        width2 = determine_spool_od(simple_tether, 300, 330, packMargin=0.33, verbose=False)

        assert width1 < width2

        # Test some edge cases
        with pytest.raises(ValueError):
            determine_spool_od(simple_tether, 0, 330, verbose=False)
        with pytest.raises(ValueError):
            determine_spool_od(simple_tether, -10, 330, verbose=False)


class TestElectricalAndPowerAnalyses:

    # --- Tests for calculate_wire_resistance --- #
    def test_calc_wire_resistance(self):

        copper_conductor = Layer("copper", "copper conductor", 0.5)
        copper_wire = Wire("electrical", "fep", "copper wire", 0.15, innerLayer=copper_conductor)
        copper_wire.assignLayerLengths(1)

        with pytest.warns(UserWarning):
            one_meter_r = np.round(calculate_wire_resistance(copper_wire, 20), 3)

        assert isclose(one_meter_r, 0.021, rel_tol=1e-3)

        with pytest.warns(UserWarning):
            r_cold = calculate_wire_resistance(copper_wire, -20)
            r_hot = calculate_wire_resistance(copper_wire, 40)

        assert r_cold < one_meter_r < r_hot

        aluminum_conductor = Layer("aluminum", "aluminum conductor", 0.5)
        aluminum_wire = Wire("electrical", "fep", "aluminum wire", 0.15, innerLayer=aluminum_conductor)
        aluminum_wire.assignLayerLengths(1)

        with pytest.warns(UserWarning):
            one_meter_r_al = np.round(calculate_wire_resistance(aluminum_wire, 20), 3) 

        assert one_meter_r_al > one_meter_r

    # --- Test errors for calculate_wire_resistance --- #
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_calc_wire_resistance_errors(self):
        layer = Layer("fep", "layer", 1)

        conductor = Layer("copper", "conductor", 0.5)


        wire1 = Wire("electrical", "fep", "wire", 0.1)
        wire1.assignLayerLengths(1)

        wire2 = Wire("electrical", "fep", "wire", 0.1, innerLayer=conductor)
        wire2.assignLayerLengths(0)

        wire3 = Wire("electrical", "fep", "wire", 0.1, innerLayer=layer)

        silica = Layer("silica", "silica", 0.1)
        optWire = Wire("optical", "fep", "wire", 0.1, innerLayer=silica)

        # Test that a wire without an inner layer causes an error
        with pytest.raises(ValueError, match="Wire must have an inner layer! This layer should be a conductive material!"):
            calculate_wire_resistance(wire1, 20)

        # Test that a non-electrical wire causes an error
        with pytest.raises(ValueError, match="A non electrical wire of type: %s was passed" % optWire.wireType):
            calculate_wire_resistance(optWire, 20)

        # Test that a type other than wire or layer causes an error
        with pytest.raises(ValueError, match="Invalid type passed! Must be of type Wire or Layer"):
            calculate_wire_resistance(2.0, 20)

        # Test that the conductor is of the right material types
        with pytest.raises(ValueError, match="Only a metal or ferromagnetic metal can be passed as a wire conductor! Given type: %s" % "polymer"):
            calculate_wire_resistance(wire3, 20)

        # Test that the conductor layer has an initialized length
        with pytest.raises(ValueError, match="A wire/layer with an uninitialized length was passed!"):
            calculate_wire_resistance(conductor, 20)
            
        # Test the above works for Wire too
        with pytest.raises(ValueError,  match="A wire/layer with an uninitialized length was passed!"):    
            calculate_wire_resistance(wire2, 20)

        # Temporarily clear all temp coefficients
        tempDict = dict(databases.tempCoeffDict)
        databases.tempCoeffDict.clear()
        conductor.assignLayerLengths(1)

        # Test that an error is thrown if we adjust temp for a material without a coefficient
        with pytest.raises(ValueError, match=re.escape("Material for temperature coeff not found in the dictionary! "
                                                       + "Run at default (20C) or add the material/coefficient to the temperature coefficient dictionary!")):
            calculate_wire_resistance(conductor, 30)

        # Restore temp coefficients
        databases.tempCoeffDict = tempDict

        # Test that an error is thrown if the temperature is less than abs zero
        with pytest.raises(ValueError, match="Temperature cannot be less than abs zero!"):
            calculate_wire_resistance(conductor, -400)


    # --- Test calculate eq resistance --- #
    @pytest.mark.usefixtures("wire_list")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_eq_resistance(self, wire_list):
        
        for wire in wire_list:
            wire.assignLayerLengths(100)

        res_list = []

        with pytest.raises(ValueError):
            calculate_equivalent_resistance(res_list, 20)

        eqResistance = 0
        last = np.inf
        for wire in wire_list:
            res_list.append(wire)

            eqResistance = calculate_equivalent_resistance(res_list, 20)
            assert eqResistance < last
            last = eqResistance
        assert np.round(eqResistance, 2) == 2.91


        res_list = []
        eqResistance = 0
        last = np.inf
        for wire in wire_list:
            res_list.append(wire)

            eqResistance = calculate_equivalent_resistance(res_list, 20, ac=True, freq=1e9)
            assert eqResistance < last
            last = eqResistance
        assert np.round(eqResistance, 2) == 124.17


    # --- Test AC Resistance --- #
    @pytest.mark.usefixtures("electrical_tether", "coaxial_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_ac_resistance(self, electrical_tether, coaxial_tether):

        # Test against a known value (wire and coax)
        coaxRes = ac_resistance(coaxial_tether.getLayerAtPath("L2"), 1e3, 20)
        assert(isclose(coaxRes, 17.41, rel_tol=1e-2))
        wireRes = ac_resistance(electrical_tether.getLayeratPath("L2M1L1"), 1e3, 20)
        assert(isclose(wireRes, 19.27, rel_tol=1e-2))

        # Test that resistance scales with frequency
        freq = 1;
        lastRes = 0;
        for i in range(1, 11, 1):
            wireRes = ac_resistance(electrical_tether.getLayeratPath("L2M1L1"), freq, 20)
            assert(wireRes > lastRes)
            lastRes = wireRes
            freq *= 1000

        # Test that freq of 1 is greater than dc 
        coaxRes = ac_resistance(coaxial_tether.getLayerAtPath("L2"), 1, 20)
        dcRes = calculate_wire_resistance(coaxial_tether.getLayerAtPath("L2"), 20)
        assert(coaxRes > dcRes)

        # Test that resistance at 0 frequency matches dc
        coaxRes = ac_resistance(coaxial_tether.getLayerAtPath("L2"), 0, 20)
        assert(coaxRes == dcRes)


    # --- Test AC Resistance Errors --- #
    @pytest.mark.usefixtures("electrical_tether", "coaxial_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_ac_resistance(self, electrical_tether):
    
        freq = 1e3

        layer = Layer("fep", "layer", 1)

        conductor = Layer("copper", "conductor", 0.5)

        wire1 = Wire("electrical", "fep", "wire", 0.1)
        wire1.assignLayerLengths(1)

        wire2 = Wire("electrical", "fep", "wire", 0.1, innerLayer=conductor)
        wire2.assignLayerLengths(0)

        wire3 = Wire("electrical", "fep", "wire", 0.1, innerLayer=layer)

        silica = Layer("silica", "silica", 0.1)
        optWire = Wire("optical", "fep", "wire", 0.1, innerLayer=silica)

        # Test that a wire without an inner layer causes an error
        with pytest.raises(ValueError, match="Wire must have an inner layer! This layer should be a conductive material!"):
            ac_resistance(wire1, freq, 20)

        # Test that a non-electrical wire causes an error
        with pytest.raises(ValueError, match="A non electrical wire of type: %s was passed" % optWire.wireType):
            ac_resistance(optWire, freq, 20)

        # Test that a type other than wire or layer causes an error
        with pytest.raises(ValueError, match="Invalid type passed! Must be of type Wire or Layer"):
            ac_resistance(2.0, freq, 20)

        # Test that the conductor is of the right material types
        with pytest.raises(ValueError, match="Only a metal or ferromagnetic metal can be passed as a wire conductor! Given type: %s" % "polymer"):
            ac_resistance(wire3, freq, 20)

        # Test that the conductor layer has an initialized length
        with pytest.raises(ValueError, match="A wire/layer with an uninitialized length was passed!"):
            ac_resistance(conductor, freq, 20)
            
        # Test the above works for Wire too
        with pytest.raises(ValueError,  match="A wire/layer with an uninitialized length was passed!"):    
            ac_resistance(wire2, freq, 20)

        # Temporarily clear all temp coefficients
        tempDict = dict(databases.tempCoeffDict)
        databases.tempCoeffDict.clear()
        conductor.assignLayerLengths(1)

        # Test that an error is thrown if we adjust temp for a material without a coefficient
        with pytest.raises(ValueError, match=re.escape("Material for temperature coeff not found in the dictionary! "
                                                       + "Run at default (20C) or add the material/coefficient to the temperature coefficient dictionary!")):
            ac_resistance(conductor, freq, 30)

        # Restore temp coefficients
        databases.tempCoeffDict = tempDict

        # Test that an error is thrown if the temperature is less than abs zero
        with pytest.raises(ValueError, match="Temperature must be greater than or equal to absolute zero!"):
            ac_resistance(conductor, freq, -400)

        with pytest.raises(ValueError, match="Frequency must greater than or equal to zero!"):
           ac_resistance(electrical_tether.getLayerAtPath("L2M1L1"), -1, 20)

        

    # --- Test DC Power Analysis ###
    @pytest.mark.usefixtures("electrical_tether", "coaxial_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_dc_power_analysis(self, electrical_tether, coaxial_tether):
        
        wireList = electrical_tether.findWires("electrical")
        send_paths = [wireList[0], wireList[2], wireList[4]]
        return_paths = [wireList[1], wireList[3], wireList[5]]

        # Test against a known figure from a tether design
        eff = dc_power_transmission_analysis(electrical_tether, 500, 500, send_paths, return_paths, False)
        assert np.round(eff[0]) == 97
        assert eff[0] > eff[1]

        # Test using layers in a coaxial cable 
        eff = dc_power_transmission_analysis(coaxial_tether, 500, 100, ["L2"], ["L4"], False)
        assert(isclose(eff[0], 99.86, rel_tol=1e-2))
        assert eff[0] > eff[1]

        # Test very large voltages
        eff = dc_power_transmission_analysis(coaxial_tether, 5e20, 100, ["L2"], ["L4"], False)
        assert eff[0] == 100.0
        assert eff[0] > eff[1]        
    

    # --- Test Errors for DC Power Analysis --- #
    @pytest.mark.usefixtures("electrical_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_dc_power_analysis_failures(self, electrical_tether):

        wireList = electrical_tether.findWires("electrical")
        send_paths = [wireList[0], wireList[2], wireList[4]]
        return_paths = [wireList[1], wireList[3], wireList[5]]

        # Test standard cases #
        with pytest.raises(ValueError, match="Tether argument must be of type RoundTetherDesign!"):
            dc_power_transmission_analysis(1, 100, 100, send_paths, return_paths, False)
        with pytest.raises(ValueError, match="Voltage must be greater than 0!"):
            dc_power_transmission_analysis(electrical_tether, 0, 100, send_paths, return_paths, False)
        with pytest.raises(ValueError, match="Desired power output must be greater than or equal to 0!"):
            dc_power_transmission_analysis(electrical_tether, 100, 0, send_paths, return_paths, False)
        with pytest.raises(ValueError, match="No send paths specified!"):
            dc_power_transmission_analysis(electrical_tether, 500, 100, [], return_paths, False)
        with pytest.raises(ValueError, match="No return paths specified!"):
            dc_power_transmission_analysis(electrical_tether, 500, 100, send_paths, [], False)
        with pytest.raises(AttributeError):
            alt_sends = ["L23"]
            dc_power_transmission_analysis(electrical_tether, 500, 100, alt_sends, return_paths, False)
        with pytest.raises(AttributeError):
            alt_recvs = ["L24"]
            dc_power_transmission_analysis(electrical_tether, 500, 100, send_paths, alt_recvs, False)
        with pytest.raises(RuntimeError):
            send_paths_dupe = [wireList[0], wireList[1], wireList[2], wireList[4]]
            dc_power_transmission_analysis(electrical_tether, 500, 100, send_paths_dupe, return_paths, False)
        with pytest.raises(RuntimeError):
            return_paths_dupe = [wireList[0], wireList[1], wireList[3], wireList[5]]
            dc_power_transmission_analysis(electrical_tether, 500, 100, send_paths, return_paths_dupe, False)
        with pytest.raises(ValueError, 
                           match=r"Solution cannot be found with given parameters! Tether Voltage drop is likely too" +
                            " high!"):
            dc_power_transmission_analysis(electrical_tether, 1, 100, send_paths, return_paths, False)

    # --- Test DC Power Analysis ###
    @pytest.mark.usefixtures("electrical_tether", "coaxial_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_single_phase_ac_transmission_analysis(self, electrical_tether):
        
        wireList = electrical_tether.findWires("electrical")
        send_paths = [wireList[0], wireList[2], wireList[4]]
        return_paths = [wireList[1], wireList[3], wireList[5]]

        # Test against a known figure from a tether design
        eff = single_phase_ac_transmission_analysis(electrical_tether, 500, 500, send_paths, return_paths, False)
        assert np.round(eff[0]) == 97
        assert eff[0] > eff[1]

        # Test using layers in a coaxial cable 
        eff = single_phase_ac_transmission_analysis(coaxial_tether, 500, 100, ["L2"], ["L4"], False)
        assert(isclose(eff[0], 99.86, rel_tol=1e-2))
        assert eff[0] > eff[1]

        # Test very large voltages
        eff = single_phase_ac_transmission_analysis(coaxial_tether, 5e20, 100, ["L2"], ["L4"], False)
        assert eff[0] == 100.0
        assert eff[0] > eff[1]        
    

    # --- Test Errors for DC Power Analysis --- #
    @pytest.mark.usefixtures("electrical_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_single_phase_ac_transmission_analysis(self, electrical_tether):

        wireList = electrical_tether.findWires("electrical")
        send_paths = [wireList[0], wireList[2], wireList[4]]
        return_paths = [wireList[1], wireList[3], wireList[5]]

        # Test standard cases #
        with pytest.raises(ValueError, match="Tether argument must be of type RoundTetherDesign!"):
            dc_power_transmission_analysis(1, 100, 100, send_paths, return_paths, False)
        with pytest.raises(ValueError, match="Voltage must be greater than 0!"):
            dc_power_transmission_analysis(electrical_tether, 0, 100, send_paths, return_paths, False)
        with pytest.raises(ValueError, match="Desired power output must be greater than or equal to 0!"):
            dc_power_transmission_analysis(electrical_tether, 100, 0, send_paths, return_paths, False)
        with pytest.raises(ValueError, match="No send paths specified!"):
            dc_power_transmission_analysis(electrical_tether, 500, 100, [], return_paths, False)
        with pytest.raises(ValueError, match="No return paths specified!"):
            dc_power_transmission_analysis(electrical_tether, 500, 100, send_paths, [], False)
        with pytest.raises(AttributeError):
            alt_sends = ["L23"]
            dc_power_transmission_analysis(electrical_tether, 500, 100, alt_sends, return_paths, False)
        with pytest.raises(AttributeError):
            alt_recvs = ["L24"]
            dc_power_transmission_analysis(electrical_tether, 500, 100, send_paths, alt_recvs, False)
        with pytest.raises(RuntimeError):
            send_paths_dupe = [wireList[0], wireList[1], wireList[2], wireList[4]]
            dc_power_transmission_analysis(electrical_tether, 500, 100, send_paths_dupe, return_paths, False)
        with pytest.raises(RuntimeError):
            return_paths_dupe = [wireList[0], wireList[1], wireList[3], wireList[5]]
            dc_power_transmission_analysis(electrical_tether, 500, 100, send_paths, return_paths_dupe, False)
        with pytest.raises(ValueError, 
                           match=r"Solution cannot be found with given parameters! Tether Voltage drop is likely too" +
                            " high!"):
            dc_power_transmission_analysis(electrical_tether, 1, 100, send_paths, return_paths, False)
        

class TestGeometryAnalyses:
    
    # Start by testing straight run strength:
    # - Test against a known value, to identify regressions in future changes
    # - Test strength does not scale with length
    # - Test layers of length zero give a strength of zero
    @pytest.mark.usefixtures("electrical_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_straight_run_strength(self, electrical_tether):

        # Test a known value from a well working version
        strength = electrical_tether.straightRunStrength(output=False)
        assert np.round(strength, 3) == 215.832

        # Make sure this does not change with length
        electrical_tether.layer.assignLayerLengths(10e3)
        strength1 = electrical_tether.straightRunStrength(output=False)
        assert np.round(strength, 3) == np.round(strength1, 3)

    # Test unhelix strength:
    @pytest.mark.usefixtures("wire_list")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_unhelix_strength(self, wire_list):
        
        # Test a known value from a well working version
        core = Layer("fep", "core", 0.25)
        memLayer = Layer("fep", "wires", 0.5, innerLayer=core, memberList=wire_list)
        abrasion = Layer("kapton", "abrasion", 0.15, innerLayer=memLayer)
        tether = RoundTetherDesign("simple tether", abrasion, 100)
        strength = tether.unhelixStrength(output=False)
        assert np.round(strength, 3) == 121.024

        # Change helix angle to 90, ensure that it matches the straight run strength 
        memLayer = Layer("fep", "wires", 0.5, innerLayer=core, memberList=wire_list, helixAngle=90)
        abrasion = Layer("kapton", "abrasion", 0.15, innerLayer=memLayer)
        tether = RoundTetherDesign("simple tether", abrasion, 100)
        str1 = tether.unhelixStrength(output=False)
        str2 = tether.straightRunStrength(output=False)
        assert np.round(str1, 6) == np.round(str2, 6)

    # Test the mass analysis (we test the individual mass functions elsewhere so just perform basic testing here)
    @pytest.mark.usefixtures("electrical_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_mass(self, electrical_tether):
        
        # Test against a known value and verify that the breakdown is correct 
        mass, massList = electrical_tether.calculateMass(breakdown=True)
        assert np.round(mass, 3) == 1164.068
        massSum = 0
        for item in massList:
            massSum += item[2]
        assert massSum == mass

        # Test that all layer lengths of zero results in zero mass
        electrical_tether.layer.assignLayerLengths(0)
        electrical_tether.length = 0
        mass, massList = electrical_tether.calculateMass(breakdown=True)
        assert mass == 0
        for item in massList:
            assert item[2] == 0

    # Test that the fill ratio lowers strength
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_fillRatio_strength(self):
        # Create some layers with different fill ratios
        layer1 = Layer("fep", "layer1", 1, fillRatio=1)
        layer2 = Layer("fep", "layer1", 1, fillRatio=0.5)
        layer3 = Layer("fep", "layer1", 1, fillRatio=0.25)

        # Create some tethers using each layer
        tet1 = RoundTetherDesign("tet1", layer1, 10)
        tet2 = RoundTetherDesign("tet2", layer2, 10)
        tet3 = RoundTetherDesign("tet3", layer3, 10)

        # Calculate the unhelix-allowed strength for each tether
        unhStr1 = tet1.unhelixStrength(output=False)
        unhStr2 = tet2.unhelixStrength(output=False)
        unhStr3 = tet3.unhelixStrength(output=False)

        # Calculate the straight-run strength for each tether
        srStr1 = tet1.straightRunStrength(output=False)
        srStr2 = tet2.straightRunStrength(output=False)
        srStr3 = tet3.straightRunStrength(output=False)

        # Make sure they all have the expected order of max->min strength
        assert (unhStr1 > unhStr2 > unhStr3)
        assert (srStr1 > srStr2 > srStr3)

        # Make sure the relationship is what we expect for each
        assert ((unhStr1 * 0.5) == unhStr2)
        assert ((unhStr1 * 0.25) == unhStr3)
        assert ((srStr1 * 0.5) == srStr2)
        assert ((srStr1 * 0.25) == srStr3)



class TestCommAnalyses:
    
    @pytest.mark.usefixtures("electrical_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_resistance_parameter(self, electrical_tether):

        # Make sure we test against a known well-working version
        wireList = electrical_tether.findWires("electrical")
        send_wire = electrical_tether.getLayerAtPath(wireList[0])
        return_wire = electrical_tether.getLayerAtPath(wireList[3])
        print("Length: %f" % send_wire.length)
        R = resistance_parameter(send_wire, return_wire, 1e3, temperature=20)
        assert isclose(R, 0.349378631, rel_tol=1e-6)

        # Ensure that this scales with frequency
        R1 = resistance_parameter(send_wire, return_wire, 1e9, temperature=20)
        assert R1 > R

    
    @pytest.mark.usefixtures("electrical_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_inductance_parameter(self, electrical_tether):

        # Make sure we test against a known well-working version
        wireList = electrical_tether.findWires("electrical")
        send_wire = electrical_tether.getLayerAtPath(wireList[0])
        return_wire = electrical_tether.getLayerAtPath(wireList[3])
        L = inductance_parameter(send_wire, return_wire)
        assert isclose(L * 1e7, 7.882, rel_tol=1e-3)
    
    @pytest.mark.usefixtures("electrical_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_capacitance_parameter(self, electrical_tether):

        # Make sure we test against a known well-working veQrsion
        wireList = electrical_tether.findWires("electrical")
        send_wire = electrical_tether.getLayerAtPath(wireList[0])
        return_wire = electrical_tether.getLayerAtPath(wireList[3])
        C = capacitance_parameter(send_wire, return_wire)
        assert np.round(C * 1e11, 3) == 2.936
    
    @pytest.mark.usefixtures("electrical_tether")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_conductance_parameter(self, electrical_tether):

        # Make sure we test against a known well-working version
        wireList = electrical_tether.findWires("electrical")
        send_wire = electrical_tether.getLayerAtPath(wireList[0])
        return_wire = electrical_tether.getLayerAtPath(wireList[3])
        G = conductance_parameter(send_wire, return_wire)
        assert np,round(G * 19, 3) == 1.625

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_comm_wirechecks(self):

        conductor1 = Layer("copper", "c1", 0.5,)
        conductor2 = Layer("aluminum", "c2", 0.5)
        conductor3 = Layer("copper", "c3", 0.75)

        wire0 = Wire("electrical", "kapton", "w0", 0.25, innerLayer=conductor1)    
        wire1 = Wire("electrical", "fep", "w1", 0.25, innerLayer=conductor1)
        wire2 = Wire("electrical", "fep", "w2", 0.15, innerLayer=conductor1)
        wire3 = Wire("electrical", "fep", "w3", 0.25, innerLayer=conductor2)
        wire4 = Wire("electrical", "fep", "w4", 0.25, innerLayer=conductor3)

        # Test that error is thrown if the same wire is passed twice
        with pytest.raises(ValueError, match="Both wires are the same! Does not refer to two conductor paths!"):
            resistance_parameter(wire0, wire0, 1)
        with pytest.raises(ValueError, match="Both wires are the same! Does not refer to two conductor paths!"):
            inductance_parameter(wire0, wire0)       
        with pytest.raises(ValueError, match="Both wires are the same! Does not refer to two conductor paths!"):
            capacitance_parameter(wire0, wire0)
        with pytest.raises(ValueError, match="Both wires are the same! Does not refer to two conductor paths!"):
            conductance_parameter(wire0, wire0)

        # Test that error is thrown on wires with different insulation material
        with pytest.raises(ValueError, match="Both wires must have the same insulation material!"):
            resistance_parameter(wire0, wire1, 1)
        with pytest.raises(ValueError, match="Both wires must have the same insulation material!"):
            inductance_parameter(wire0, wire1)       
        with pytest.raises(ValueError, match="Both wires must have the same insulation material!"):
            capacitance_parameter(wire0, wire1)
        with pytest.raises(ValueError, match="Both wires must have the same insulation material!"):
            conductance_parameter(wire0, wire1)

        # - Value error if insulation thickneses differ
        with pytest.raises(ValueError, match="Both wires must have the same thickness of insulation!"):
            resistance_parameter(wire2, wire1, 1)
        with pytest.raises(ValueError, match="Both wires must have the same thickness of insulation!"):
            inductance_parameter(wire2, wire1)       
        with pytest.raises(ValueError, match="Both wires must have the same thickness of insulation!"):
            capacitance_parameter(wire2, wire1)
        with pytest.raises(ValueError, match="Both wires must have the same thickness of insulation!"):
            conductance_parameter(wire2, wire1)

        # - Value error if conductor materials differ
        with pytest.raises(ValueError, match="Both wires must have the same conductor material!"):
            resistance_parameter(wire3, wire1, 1)
        with pytest.raises(ValueError, match="Both wires must have the same conductor material!"):
            inductance_parameter(wire3, wire1)       
        with pytest.raises(ValueError, match="Both wires must have the same conductor material!"):
            capacitance_parameter(wire3, wire1)
        with pytest.raises(ValueError, match="Both wires must have the same conductor material!"):
            conductance_parameter(wire3, wire1)

        # - Value error if conductor sizes differ
        with pytest.raises(ValueError, match="Both wires must have the same size conductors!"):
            resistance_parameter(wire4, wire1, 1)
        with pytest.raises(ValueError, match="Both wires must have the same size conductors!"):
            inductance_parameter(wire4, wire1)       
        with pytest.raises(ValueError, match="Both wires must have the same size conductors!"):
            capacitance_parameter(wire4, wire1)
        with pytest.raises(ValueError, match="Both wires must have the same size conductors!"):
            conductance_parameter(wire4, wire1)

        # - Value error Negative frequency is passed
        with pytest.raises(ValueError, match="Frequency must be greater than zero for impedance calcs!"):
            resistance_parameter(wire1, deepcopy(wire1), 0)


