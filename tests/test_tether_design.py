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
from tether_analysis.TetherDesign import *
from databases.DatabaseControl import *
from math import isclose

# Grab database object and make sure we are using the test databases for these tests #
databaseCtrl = DatabaseControl()
databaseCtrl.load_from_folder(pathlib.Path.cwd() / 'test_databases')


### Define test fixtures ###
@pytest.fixture
def basic_layer():
    return Layer("kapton", "kapton_layer", 1.1)

@pytest.fixture
def member_list():
    
    mem1 = Wire("optical", "silica", "basic_fo", 0.75)

    conductor = Layer("copper", "copper_layer", 1)
    insulator = Wire("electrical", "fep", "basic_wire", 0.4, innerLayer=conductor)

    memList = []
    for i in range(0, 4):
        memList.append(mem1)
        memList.append(insulator)

    return memList

@pytest.fixture
def member_layer(member_list):
        return Layer("kapton", "kapton_layer", 1.5, memberList=member_list)

@pytest.fixture
def wire_entry():

    return DB.searchEntry(databaseCtrl.wire_db_dict["Teledyne-HV-FEP-2023"], {"Part #" : "178-8073"})

@pytest.fixture
def sm_optical_entry():

    return DB.searchEntry(databaseCtrl.optical_db_dict["LindenPhotonics-RadHard-2023"], {"Part #" : "LINDEN-SPE-7209"})

@pytest.fixture
def mm_optical_entry():

    return DB.searchEntry(databaseCtrl.optical_db_dict["LindenPhotonics-RadHard-2023"], {"Part #" : "LINDEN-SPE-7088"})

@pytest.fixture
def multiple_member_list():
    wire_outer = Layer("etfe", "etfe_layer", 1)
    fiber_outer = Layer("hytrel_dupont_5526", "hytrel_layer", 0.25)

    wireList = []

    for i in range(0, 6):
        wireList.append(wire_outer)
        wireList.append(fiber_outer)

    return wireList

@pytest.fixture
def small_large_member_list():
    wire_outer = Layer("etfe", "etfe_layer", 1)
    fiber_outer = Layer("hytrel_dupont_5526", "hytrel_layer", 0.05)

    wireList = []

    for i in range(0, 6):
        wireList.append(wire_outer)
        wireList.append(fiber_outer)

    return wireList

@pytest.fixture
def large_qty_member_list():
    wire_outer = Layer("etfe", "etfe_layer", 1)
    fiber_outer = Layer("hytrel_dupont_5526", "hytrel_layer", 0.05)

    wireList = []

    for i in range(0, 18):
        wireList.append(wire_outer)
        wireList.append(fiber_outer)

    return wireList

@pytest.fixture
def very_large_member_list():
    wire_outer = Layer("etfe", "etfe_layer", 1)
    fiber_outer = Layer("hytrel_dupont_5526", "hytrel_layer", 25)

    wireList = []

    for i in range(0, 18):
        wireList.append(wire_outer)
        wireList.append(fiber_outer)

    return wireList

@pytest.fixture
def complex_tree():
    angle=60

    wire_inner = Layer("copper", "copper_layer", 2)
    wire_outer = Wire("electrical", "etfe", "etfe_ins_wire", 1, innerLayer=wire_inner)

    fiber_inner = Layer("silica", "silica_layer", 1.5)
    fiber_outer = Wire("optical", "hytrel_dupont_5526", "hytrel_layer", 0.25, innerLayer=fiber_inner)

    wireList = []

    for i in range(0, 3):
        wireList.append(wire_outer)
        wireList.append(fiber_outer)

    core = Layer("nickel", "core_component", 1.5, color="black")
    initStr = Layer("aluminum", "aluminum_layer", 1, innerLayer=core, color="white")
    memLayer = Layer("fep", "mem_layer", 1, innerLayer=initStr, memberList=wireList, helixAngle=angle)

    cableList = []
    for i in range(0, 3):
        cableList.append(memLayer)

    act_core = Layer("aluminum", "aluminum_core", 3, color="white")
    strength_layer = Layer("vectran_kuraray_um", "inner_strength_layer", 1, innerLayer=act_core)
    cableLayer = Layer("polyurethane", "cable_layer", 1, innerLayer=strength_layer, memberList=cableList, 
                       helixAngle=angle)
    outer_strength_layer = Layer("vectran_kuraray_um", "outer_strength_layer", 1, innerLayer=cableLayer)
    abrasionLayer = Layer("etfe", "abrasion_layer", 2.5, innerLayer=outer_strength_layer, color="olivedrab")

    outer_mem_layer = Layer("fep", "outer_mem_layer", 10, innerLayer=strength_layer, 
                            memberList=[abrasionLayer, abrasionLayer, abrasionLayer, abrasionLayer])

    return outer_mem_layer


### Define test suite for Layer ###
class TestLayer:

    ### --- Tests for the constructor --- ###
    def test_constructor_standalone(self):
        layer = Layer("etfe", "etfe_layer", 1.1)

        assert layer.layerMaterial == "etfe"
        assert layer.innerDiameter == 0
        assert layer.innerRadius == 0
        assert layer.outerDiameter == 2.2
        assert layer.outerRadius == 1.1

    @pytest.mark.usefixtures("basic_layer")
    def test_constructor_nested(self, basic_layer):

        
        outLayer = Layer("fep", "fep_layer", 0.5, innerLayer=basic_layer)

        assert outLayer.innerDiameter == basic_layer.outerDiameter
        assert outLayer.innerRadius == basic_layer.outerRadius
        assert outLayer.outerDiameter == basic_layer.outerDiameter + (2*outLayer.layerThickness)
        assert outLayer.outerRadius == basic_layer.outerRadius + outLayer.layerThickness
        assert outLayer.innerLayer == basic_layer

    @pytest.mark.usefixtures("basic_layer", "member_list")
    def test_constructor_members(self, basic_layer, member_list):
        
        layer = Layer("etfe", "etfe_layer", 1, innerLayer=basic_layer, memberList=member_list)

        maxDiam = 0
        for mem in member_list:
            if mem.outerRadius > maxDiam:
                maxDiam = mem.outerDiameter

        assert layer.outerRadius >= mem.outerDiameter + basic_layer.outerRadius
        assert len(layer.memberList) == len(member_list)
        assert not layer.check_member_collisions()

    def test_invalid_thickness(self):
           
        with pytest.raises(ValueError):
            Layer("etfe", "etfe_layer", 0)
        
        with pytest.raises(ValueError):
            Layer("kapton", "kapton_layer", -1)

    def test_fill_ratio(self):

        layer1 =  Layer("kapton", "kapton_layer", 1.1, fillRatio=0.9)
        layer2 =  Layer("kapton", "kapton_layer", 1.1)

        mass1 = layer1.calculateMass(10)
        mass2 = layer2.calculateMass(10) * 0.9

        assert isclose(mass1, mass2, rel_tol=1e-6)

        with pytest.raises(ValueError, match=r"Fill ratio must be between \(0, 1]!"):
            Layer("kapton", "kapton_layer", 1.1, fillRatio=-1)
        with pytest.raises(ValueError, match=r"Fill ratio must be between \(0, 1]!"):
            Layer("kapton", "kapton_layer", 1.1, fillRatio=0)
        with pytest.raises(ValueError, match=r"Fill ratio must be between \(0, 1]!"):
            Layer("kapton", "kapton_layer", 1.1, fillRatio=1.1)
        

    @pytest.mark.usefixtures("basic_layer")
    def test_warn_copy_off(self, basic_layer):

        with pytest.warns(UserWarning):
            Layer("fep", "fep_layer", 0.5, innerLayer=basic_layer, copy=False)


    # --- Test Define Layer Geometry --- #
    @pytest.mark.usefixtures("basic_layer", "member_list")
    def test_geometry_updated(self, basic_layer, member_list):

        layerDiam = basic_layer.outerDiameter
        startRad = basic_layer.memberStartRadius
        innDiam = basic_layer.innerDiameter
        thick = basic_layer.layerThickness


        for mem in member_list:
            basic_layer.memberList.append(deepcopy(mem))

        basic_layer.define_layer_geometry()

        assert basic_layer.outerDiameter > layerDiam
        assert basic_layer.memberStartRadius > startRad
        assert innDiam == basic_layer.innerDiameter
        assert thick < basic_layer.layerThickness
        assert not basic_layer.check_member_collisions()


    @pytest.mark.usefixtures("basic_layer", "member_list")
    def test_iteration_limit(self, basic_layer, member_list):
        
        basic_layer.iter_limit = 0

        for mem in member_list:
            basic_layer.memberList.append(deepcopy(mem))

        with pytest.raises(RuntimeError):
            basic_layer.define_layer_geometry()

    @pytest.mark.usefixtures("multiple_member_list", "small_large_member_list")
    def test_helix_very_small_members(self, multiple_member_list, small_large_member_list):

        # These two designs should have the same OD #
        core1 = Layer("steel", "core component", 0.5, color="black")
        initStr1 = Layer("aluminum", "aluminum_layer", 1, innerLayer=core1, color="white")
        memLayer1 = Layer("fep", "fep_layer", 1, innerLayer=initStr1, memberList=multiple_member_list)

        core2 = Layer("steel", "core component", 0.5, color="black")
        initStr2 = Layer("aluminum", "aluminum_layer", 1, innerLayer=core2, color="white")
        memLayer2 = Layer("fep", "fep_layer", 1, innerLayer=initStr2, memberList=small_large_member_list)

        assert memLayer1.outerDiameter == memLayer2.outerDiameter

        maxMemDiam = 0
        for mem in memLayer1.memberList:
            if maxMemDiam < mem.outerDiameter:
                maxMemDiam = mem.outerDiameter
        assert memLayer1.layerThickness == maxMemDiam

        maxMemDiam = 0
        for mem in memLayer2.memberList:
            if maxMemDiam < mem.outerDiameter:
                maxMemDiam = mem.outerDiameter
        assert memLayer2.layerThickness == maxMemDiam

        assert not memLayer1.check_member_collisions()
        assert not memLayer2.check_member_collisions()

    @pytest.mark.usefixtures("very_large_member_list")
    def test_very_large_members(self, very_large_member_list):
        core1 = Layer("steel", "core_component", 0.5, color="black")
        initStr1 = Layer("aluminum", "aluminum_layer", 1, innerLayer=core1, color="white")
        memLayer1 = Layer("fep", "fep_layer", 1, innerLayer=initStr1, memberList=very_large_member_list, step=0.01)

        maxMemDiam = 0
        for mem in memLayer1.memberList:
            if maxMemDiam < mem.outerDiameter:
                maxMemDiam = mem.outerDiameter
        assert memLayer1.layerThickness >= maxMemDiam

        assert not memLayer1.check_member_collisions()

    @pytest.mark.usefixtures("large_qty_member_list")
    def test_large_qty_members(self, large_qty_member_list):
        core1 = Layer("steel", "core_component", 0.5, color="black")
        initStr1 = Layer("aluminum", "aluminum_layer", 1, innerLayer=core1, color="white")
        memLayer1 = Layer("fep", "fep_layer", 1, innerLayer=initStr1, memberList=large_qty_member_list)

        maxMemDiam = 0
        for mem in memLayer1.memberList:
            if maxMemDiam < mem.outerDiameter:
                maxMemDiam = mem.outerDiameter
        assert memLayer1.layerThickness >= maxMemDiam

        assert not memLayer1.check_member_collisions()



    # --- Test update member positions --- #
    @pytest.mark.usefixtures("member_layer")
    def test_member_positions_updated(self, member_layer):
        
        member_layer.x += 1
        member_layer.y -= 1
        member_layer.polarCoords = cart_to_polar([member_layer.x, member_layer.y])

        memPosList = []

        for member in member_layer.memberList:
            memPosList.append([member.x, member.y])

        member_layer.update_member_positions()

        for idx, member in enumerate(member_layer.memberList):

            assert memPosList[idx][0] + 1 == member.x
            assert memPosList[idx][1] - 1 == member.y

    # --- Test check_member_collisions --- #
    @pytest.mark.usefixtures("member_layer")
    def test_collision_check(self, member_layer):
        
        assert not member_layer.check_member_collisions()

        member_layer.memberList[0].x = member_layer.memberList[1].x
        member_layer.memberList[0].y = member_layer.memberList[1].y

        assert member_layer.check_member_collisions()


    # --- Test calculateMass --- #
    def test_mass_correctness(self):
        
        layer1 = Layer("copper", "copper_layer", 1)
        mass = layer1.calculateMass(1)

        ### Make sure the mass is equal to 1m of 2mm diam copper wire, to the thousandth place ###
        assert round(mass, 3) == 28.149

        layer2 = Layer("copper", "copper_layer", 1)
        layer3 = Layer("copper", "copper_layer", 1, innerLayer=layer2)

        mass = layer3.calculateMass(1)

        ### Make sure the mass of nested layers is equal to 1mm of 4mm diam copper wire, to the thousandth place ###
        assert round(mass, 3) == 112.595


    @pytest.mark.usefixtures("multiple_member_list")
    def test_helix_angles(self, multiple_member_list):
        core = Layer("eglass", "eglass_layer", 0.5)

        ### Check errors thrown for invalid helix angles ###
        with pytest.raises(ValueError):
            layer = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=0)
            layer.calculateMass(1)
        with pytest.raises(ValueError):
            layer = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=91)
            layer.calculateMass(1)
        with pytest.raises(ValueError):
            layer = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=-4)
            layer.calculateMass(1)

        ### Check for error if we have an angle with a lead that's too small ###
        with pytest.raises(ValueError):
            layer = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=1)
            layer.calculateMass(1)

        ### Make sure increasing the angle results in decreasing mass ###
        memLayer1 = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=89)
        memLayer2 = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=90)
        memLayer3 = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=30)

        ### Test both a small and large difference in helix angle ###
        assert memLayer1.calculateMass(1) > memLayer2.calculateMass(1)
        assert memLayer3.calculateMass(1) > memLayer1.calculateMass(1)

    @pytest.mark.usefixtures("multiple_member_list")
    def test_breakdown_list(self, multiple_member_list):

        ### Setup layers and count layers ###
        core = Layer("eglass", "eglass_layer", 0.5)
        layer = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=60)
        numLayers = layer.countLayers()

        ### Calculate masses to compare ###
        breakDown = []
        mass = layer.calculateMass(1, breakDownList=breakDown)
        altMass = 0
        for item in breakDown:
            altMass += item[2]

        ### Make sure the # of layers and the total masses are consistent ###
        assert len(breakDown) == numLayers
        assert mass == altMass

    def test_length_edge_cases(self):

        layer1 = Layer("copper", "copper_layer", 1)

        assert round(layer1.calculateMass(1e6), 3) == 28148670.176
        assert layer1.calculateMass(0) == 0
        with pytest.raises(ValueError):
            layer1.calculateMass(-1)

    # --- Test  assign layer lengths --- #
    @pytest.mark.usefixtures("multiple_member_list")
    def test_assign_lengths(self, multiple_member_list):
        core = Layer("eglass", "eglass_layer", 0.5)
        layer = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=90)

        ### Check that the layers all have length of one when helix angle is 90 ###
        layer.assignLayerLengths(1)
        assert layer.length == 1
        assert layer.innerLayer.length == 1
        for mem in layer.memberList:
            assert mem.length == 1

        ### Check that member's length increases when the length is greater than one
        layer2 = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=60)
        layer2.assignLayerLengths(1)
        for mem in layer2.memberList:
            assert mem.length > 1

        ### Make sure an error is thrown if we are passed a negative value ###
        with pytest.raises(ValueError):
            layer2.assignLayerLengths(-1)


    # --- Test check_layer_collisions --- #
    @pytest.mark.usefixtures("multiple_member_list")
    def test_check_layer_collisions(self, multiple_member_list):

        core = Layer("eglass", "eglass_layer", 0.5)
        layer1 = Layer("kapton", "kapton_layer", 1, innerLayer=core, memberList=multiple_member_list, helixAngle=60)

        layer1.memberList[2].polarCoords[0] += 1

        with pytest.raises(ValueError):
            layer1.check_layer_collisions()


### Define test suite for Wire ###
class TestWire:

    # --- Tests for wire constructor --- #
    def test_invalid_wireType(self):
        with pytest.raises(ValueError):
            Wire("impossible", "material", "basic_wire", 3)

    def test_basic_wire(self):

        conductor = Layer("copper", "copper_layer", 0.25)
        wire = Wire("electrical", "fep", "basic_wire", 0.5, innerLayer=conductor)

        assert wire.outerRadius == 0.75
        assert wire.innerRadius == 0.25
        assert wire.outerDiameter == 1.5
        assert wire.innerDiameter == 0.5

    # --- Test from_entry class method --- #

    def test_invalid_wiretype_entry(self):
        with pytest.raises(ValueError):
            Wire.from_entry("impossible", None)

    @pytest.mark.usefixtures("wire_entry")
    def test_wire_from_entry(self, wire_entry):

        with pytest.warns(UserWarning):
            wire = Wire.from_entry("electrical", wire_entry)

        assert wire.outerDiameter == 1.4
        assert wire.innerLayer.outerDiameter == 0.8

        assert wire.layerMaterial == "fep"
        assert wire.innerLayer.layerMaterial == "copper"

    def test_wire_entry_from_thickness(self):
        
        temp1 = databaseCtrl.wire_db_dict["Teledyne-HV-FEP-2023"].iloc[10, 5]
        temp2 = databaseCtrl.wire_db_dict["Teledyne-HV-FEP-2023"].iloc[10, 6]
        databaseCtrl.wire_db_dict["Teledyne-HV-FEP-2023"].iloc[10, 5] = 0.3
        databaseCtrl.wire_db_dict["Teledyne-HV-FEP-2023"].iloc[10, 6] = None

        entry = searchEntry(databaseCtrl.wire_db_dict["Teledyne-HV-FEP-2023"], {"Part #" : "178-8073"})

        with pytest.warns(UserWarning):
            wire = Wire.from_entry("electrical", entry)

        assert wire.outerDiameter == 1.4
        assert wire.innerLayer.outerDiameter == 0.8

        assert wire.layerMaterial == "fep"
        assert wire.innerLayer.layerMaterial == "copper"

        databaseCtrl.wire_db_dict["Teledyne-HV-FEP-2023"].iloc[10, 5] = temp1
        databaseCtrl.wire_db_dict["Teledyne-HV-FEP-2023"].iloc[10, 6] = temp2

    @pytest.mark.usefixtures("sm_optical_entry")
    def test_from_sm_optical_entry(self, sm_optical_entry):

        fo = Wire.from_entry("optical", sm_optical_entry)

        assert fo.outerDiameter == 2.2
        assert isclose(fo.innerLayer.outerDiameter, 125 * 1e-3, rel_tol=1e-15)
        assert isclose(fo.innerLayer.innerLayer.outerDiameter, 9 * 1e-3, rel_tol=1e-15)

        assert fo.layerMaterial == "LCP"
        assert fo.innerLayer.layerMaterial == "silica"
        assert fo.innerLayer.innerLayer.layerMaterial == "silica"

    @pytest.mark.usefixtures("mm_optical_entry")
    def test_from_mm_optical_entry(self, mm_optical_entry):
        
        fo = Wire.from_entry("optical", mm_optical_entry)

        assert fo.outerDiameter == 0.9
        assert isclose(fo.innerLayer.outerDiameter, 125 * 1e-3, rel_tol=1e-15)
        assert isclose(fo.innerLayer.innerLayer.outerDiameter, 50 * 1e-3, rel_tol=1e-15)

        assert fo.layerMaterial == "LCP"
        assert fo.innerLayer.layerMaterial == "silica"
        assert fo.innerLayer.innerLayer.layerMaterial == "silica"

    @pytest.mark.usefixtures("wire_entry")
    def test_missing_wire_entry_values(self, wire_entry):

        entry1 = deepcopy(wire_entry)
        entry1["insulation_od"] = None
        entry1["insulaion_thickness"] = None
        with pytest.raises(ValueError):
            Wire.from_entry("electrical", entry1)

        entry2 = deepcopy(wire_entry)
        entry2["conductor_od"] = None
        with pytest.raises(ValueError):
            Wire.from_entry("electrical", entry2)

        entry3 = deepcopy(wire_entry)
        entry3["insulation_material"] = None
        with pytest.raises(ValueError):
            Wire.from_entry("electrical", entry3)

        entry4 = deepcopy(wire_entry)
        entry4["conductor_material"] = None
        with pytest.raises(ValueError):
            Wire.from_entry("electrical", entry4)

    @pytest.mark.usefixtures("sm_optical_entry", "mm_optical_entry")
    def test_missing_optical_entry_values(self, sm_optical_entry, mm_optical_entry):

        entry1 = deepcopy(sm_optical_entry)
        entry1["type"] = None
        with pytest.raises(ValueError):
            Wire.from_entry("optical", entry1)

        entry2 = deepcopy(mm_optical_entry)
        entry2["core_od"] = None
        with pytest.raises(ValueError):
            Wire.from_entry("optical", entry2)

        entry3 = deepcopy(mm_optical_entry)
        entry3["cladding_od"] = None
        with pytest.raises(ValueError):
            Wire.from_entry("optical", entry3)

        entry4 = deepcopy(sm_optical_entry)
        entry4["jacket_od"] = None
        with pytest.raises(ValueError):
            Wire.from_entry("optical", entry4)

        entry5 = deepcopy(mm_optical_entry)
        entry5["jacket_material"] = None
        with pytest.raises(ValueError):
            Wire.from_entry("optical", entry5)

    def test_invalid_optical_type(self):

        temp1 = databaseCtrl.optical_db_dict["LindenPhotonics-RadHard-2023"].iloc[0, 1]
        databaseCtrl.optical_db_dict["LindenPhotonics-RadHard-2023"].iloc[0, 1] = "impossible"
        
        entry = searchEntry(databaseCtrl.optical_db_dict["LindenPhotonics-RadHard-2023"], 
                            {"Part #" : "LINDEN-SPE-7079"})

        with pytest.raises(ValueError):
            Wire.from_entry("optical", entry)

        databaseCtrl.optical_db_dict["LindenPhotonics-RadHard-2023"].iloc[0, 1] = temp1


### Test suite for RoundTetherDesign ###
class TestRoundTetherDesign:

    # --- Test Constructor --- #
    @pytest.mark.usefixtures("complex_tree")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_round_tether_constructor(self, complex_tree):

        tether = RoundTetherDesign("tether", complex_tree, 100)

        pathSet = set()
        layerCount = tether.layer.countLayers()
        layerList = []
        build_mem_list(tether.layer, layerList)

        assert len(layerList) == layerCount

        for layer in layerList:
            pathSet.add(layer.layerPath)
            assert layer.length >= tether.length
            
        assert len(pathSet) == layerCount

        with pytest.raises(ValueError):
            RoundTetherDesign("tether", complex_tree, -1)

    # --- Test get layer at path --- #
    @pytest.mark.usefixtures("complex_tree")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_get_layer_at_path(self, complex_tree):

        tether = RoundTetherDesign("tether", complex_tree, 1)
        layer = tether.getLayerAtPath("L1M2L3M2L1M3L2")

        expLayer = tether.layer.memberList[1].innerLayer.innerLayer.memberList[1].memberList[2].innerLayer

        assert layer is expLayer

        expLayer.layerPath = "impossible"

        with pytest.raises(ValueError):
            tether.getLayerAtPath("L1M2L3M2L1M3L2")
        expLayer.layerPath = "L1M2L3M2L1M3L2"

        with pytest.raises(Exception):
            tether.getLayerAtPath("L13")
        with pytest.raises(Exception):
            tether.getLayerAtPath("L1M5")


    # --- Test findWires --- #
    @pytest.mark.usefixtures("complex_tree", "basic_layer")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_find_wires(self, complex_tree, basic_layer):

        tether = RoundTetherDesign("tether", complex_tree, 1)
        elecList = tether.findWires("electrical")
        optList = tether.findWires("optical")

        assert len(elecList) == 36
        assert len(optList) == 36

        with pytest.raises(ValueError):
            tether.findWires("impossible")

        tether = RoundTetherDesign("tether", basic_layer, 1)
        wireList1 = tether.findWires("electrical")
        wireList2 = tether.findWires("optical")

        assert len(wireList1) == len(wireList2) == 0

    # --- Test buildLayerLists --- #
    @pytest.mark.usefixtures("complex_tree")
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_build_layer_list(self, complex_tree):

        tether = RoundTetherDesign("tether", complex_tree, 100)
        memList = []
        build_mem_list(tether.layer, memList)

        pathList = []
        for layer in memList:
            pathList.append(layer.layerPath)

        layers = tether.buildLayerList(pathList)
    
        for layer in layers:
            assert layer.layerPath in pathList

        assert len(layers) == len(pathList)


# --- Helpers --- #
def build_mem_list(layer, list):

    list.append(layer)

    for mem in layer.memberList:
        build_mem_list(mem, list)

    if layer.innerLayer is not None:
        build_mem_list(layer.innerLayer, list)


