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
from databases.DatabaseControl import *
import pandas as pd
import pathlib
from math import isclose

# Grab database object and make sure we are using the test databases for these tests #
databaseCtrl = DatabaseControl()
databaseCtrl.load_from_folder(pathlib.Path.cwd() / 'test_databases')

@pytest.fixture
def wire_database():
    path = pathlib.Path.cwd() / 'test_databases'
    db = pd.read_csv(path / 'Teledyne-HV-FEP-2023_WireDB.csv').replace(np.nan, None, regex=True)
    return ['Teledyne-HV-FEP-2023', db]

@pytest.fixture
def optical_database():
    path = pathlib.Path.cwd() / 'test_databases'
    db = pd.read_csv(path / 'LindenPhotonics-RadHard-2023_OpticalDB.csv').replace(np.nan, None, regex=True)
    return ['LindenPhotonics-RadHard-2023', db]

@pytest.fixture
def material_database():
    path = pathlib.Path.cwd() / 'test_databases'
    db = pd.read_csv(path / 'material_DB.csv').replace(np.nan, None, regex=True)
    return db

@pytest.fixture
def material_entry():
    return databaseCtrl.get_material_entry("metal")

### Test DatabaseControl class ###

class TestDatabaseControl:

    # --- Functions to test the constructor --- #
    def test_material_db_found(self):
        assert databaseCtrl.material_db is not None

    def test_wire_db_found(self):
        assert len(databaseCtrl.wire_db_list) > 0
        assert len(databaseCtrl.wire_db_dict) > 0

    def test_optical_db_found(self):
        assert len(databaseCtrl.optical_db_list) > 0
        assert len(databaseCtrl.optical_db_dict) > 0

    # --- Function to test that databases were properly processed --- #
    @pytest.mark.usefixtures("wire_database")
    def test_wire_convert(self, wire_database):

        for col in wire_database[1].columns:
            
            if col not in databaseCtrl.wireDefaultUnits.keys():
                continue

            default_unit = databaseCtrl.wireDefaultUnits[col]
            unit = wire_database[1][col][0]
            
            multiplier = build_multiplier(unit, default_unit)

            for idx, val in enumerate(wire_database[1][col][1:]):
                if val is None:
                    continue
                else:

                    setVal = databaseCtrl.wire_db_dict[wire_database[0]][col].iloc[idx]
                    val = float(val)
                    assert (setVal / val) == multiplier

    @pytest.mark.usefixtures("optical_database")
    def test_optical_convert(self, optical_database):

        for col in optical_database[1].columns:
            
            if col not in databaseCtrl.opticalDefaultUnits.keys():
                continue

            default_unit = databaseCtrl.opticalDefaultUnits[col]
            unit = optical_database[1][col][0]
            
            multiplier = build_multiplier(unit, default_unit)

            for idx, val in enumerate(optical_database[1][col][1:]):
                if val is None:
                    continue

                setVal = databaseCtrl.optical_db_dict[optical_database[0]][col].iloc[idx]
                val = float(val)
                assert (setVal / val) == multiplier

    @pytest.mark.usefixtures("material_database")
    def test_material_convert(self, material_database):
        
        for col in material_database.columns:

            if col not in databaseCtrl.materialDefaultUnits.keys():
                continue

            default_unit = databaseCtrl.materialDefaultUnits[col]

            for idx, val in enumerate(material_database[col]):

                if val is None:
                    continue

                unit_col = col + "_unit"

                setVal = databaseCtrl.material_db[col].iloc[idx]
                setUnit = databaseCtrl.material_db[unit_col].iloc[idx]
                unit = material_database[unit_col].iloc[idx]

                multiplier = build_multiplier(unit, default_unit)
                val = float(val)
                
                assert setUnit == default_unit
                assert (val * multiplier) == setVal
    
    @pytest.mark.usefixtures("wire_database")
    def test_invalid_wire_type(self, wire_database):
        
        with pytest.raises(ValueError):
            databaseCtrl.process_database(wire_database[1], wire_type="impossible")

    @pytest.mark.usefixtures("wire_database")
    def test_missing_unit_process(self, wire_database):
        
        temp = wire_database[1].iloc[0, 6]
        wire_database[1].iloc[0, 6] = None

        with pytest.raises(ValueError):
            databaseCtrl.process_database(wire_database[1])

        wire_database[1].iloc[0, 6] = temp

    ## Build override material tests ##
    #TODO: not yet implemented

    ## Add override material tests ##
    #TODO: not yet implemented

    # --- Functions to test get material entry --- #
    def test_no_material_entry(self):
        mat = "impossible_material"

        with pytest.raises(ValueError, match="Passed material: %s not found in material database!" % mat):
            databaseCtrl.get_material_entry(mat)

    def test_empty_material_type(self):

        ### set material type to None ###
        temp = databaseCtrl.material_db.iloc[7, 1]
        databaseCtrl.material_db.iloc[7, 1] = None

        with pytest.raises(ValueError):
            databaseCtrl.get_material_entry("acrylic_low")

        databaseCtrl.material_db.iloc[7, 1] = temp

    def test_missing_filled(self):
        metal_entry = databaseCtrl.get_material_entry("metal")

        with pytest.warns(UserWarning):
            copper_entry = databaseCtrl.get_material_entry("copper")

            assert copper_entry["min_temp"].item() == metal_entry["min_temp"].item()
            assert copper_entry["max_temp"].item() == metal_entry["max_temp"].item()
            assert copper_entry["dielectric_strength"].item() == metal_entry["dielectric_strength"].item()
            assert copper_entry["dielectric_constant"].item() == metal_entry["dielectric_constant"].item()

    def test_invalid_representative_entry_val(self):

        # Set polymer representative material and acrlyic_high density to None
        temp1 = databaseCtrl.material_db.iloc[0, 2]
        temp2 = databaseCtrl.material_db.iloc[7, 2]
        databaseCtrl.material_db.iloc[0, 2] = None
        databaseCtrl.material_db.iloc[7, 2] = None

        with pytest.raises(ValueError):
            databaseCtrl.get_material_entry("acrylic_low")

        databaseCtrl.material_db.iloc[0, 2] = temp1
        databaseCtrl.material_db.iloc[7, 2] = temp2

    def test_invalid_representative_entry_unit(self):

        # Set polymer representative material density unit and acrlyic_high density to None
        temp1 = databaseCtrl.material_db.iloc[0, 3]
        temp2 = databaseCtrl.material_db.iloc[7, 2]
        databaseCtrl.material_db.iloc[0, 3] = None
        databaseCtrl.material_db.iloc[7, 2] = None

        with pytest.raises(ValueError):
            databaseCtrl.get_material_entry("acrylic_low")

        databaseCtrl.material_db.iloc[0, 3] = temp1
        databaseCtrl.material_db.iloc[7, 2] = temp2

    def test_warn_fill_false(self):

        with pytest.warns(UserWarning):
            databaseCtrl.get_material_entry("acrylic_high", fill_missing=False)


### Test helper functions ###

class TestDatabaseHelpers:

    # --- Functions to test searchEntry --- #
    @pytest.mark.usefixtures("optical_database")
    def test_valid_entry_returned(self, optical_database):

        db = databaseCtrl.process_database(optical_database[1], wire_type="optical")

        entry1 = db[(db["Part #"] == "LINDEN-SPE-7210") & (db["type"] == "multimode")]
        entry2 = searchEntry(databaseCtrl.optical_db_dict["LindenPhotonics-RadHard-2023"],
                              {"Part #" : "LINDEN-SPE-7210", "type":"multimode"})

        assert entry1.equals(entry2)

    def test_empty_returned(self):

        entry = searchEntry(databaseCtrl.optical_db_dict["LindenPhotonics-RadHard-2023"], 
                            {"Part #" : "LINDEN-SPE-7210", "type":"othermode"})

        assert entry.empty

    def test_error_empty_query_dict(self):
        with pytest.raises(ValueError):
            searchEntry(databaseCtrl.optical_db_dict["LindenPhotonics-RadHard-2023"], {})

    # --- Functions to test get_material_property --- #
    @pytest.mark.usefixtures("material_entry")
    def test_get_material_property(self, material_entry):
        expected = 0.0089588
        received = get_material_property(material_entry, "density")

        assert isclose(expected, received, rel_tol=1e-9)

    @pytest.mark.usefixtures("material_entry")
    def test_invalid_property(self, material_entry):
        with pytest.raises(KeyError):
            get_material_property(material_entry, "impossible")
                
    # --- Functions to test replace_with_mult --- #
    def test_correct_mults(self):
        
        retStr1 = replace_with_mult("kg", "g", "kg")
        retStr2 = replace_with_mult("mm", "m", "mm")
        retStr3 = replace_with_mult("MPa", "GPa", "MPa / F")

        assert retStr1 == "1000.0"
        assert retStr2 == "0.001"
        assert retStr3 == "0.001 / F"

    # --- Functions to test get_conv_func --- #
    def test_conv_func(self):
        assert mm_convert == get_conv_func("mm")
        assert m_convert == get_conv_func("m")
        assert km_convert == get_conv_func("km")
        assert v_convert == get_conv_func("V")
        assert kg_convert == get_conv_func("kg")
        assert g_convert == get_conv_func("g")
        assert gpa_convert == get_conv_func("GPa")
        assert ohm_convert == get_conv_func("ohm")
        assert C_convert == get_conv_func("C")
        assert F_convert == get_conv_func("F")
        assert bar_convert == get_conv_func("bar")

    def test_invalid_unit(self):
        with pytest.raises(ValueError):
            get_conv_func("impossible")

    # --- Functions to test build_multiplier --- #
    def test_correct_multiplier(self):
        unit_str = "kg * cm^2 / 12 + F"
        default_str = "g * m^2 / 12 + C"

        multiplier = build_multiplier(unit_str, default_str)
        assert isclose((1000 * 0.01**2 / 12 + 0.55555555555), multiplier, rel_tol=1e-9)
    
    def test_diff_unit_num(self):
        with pytest.raises(ValueError):
            build_multiplier("g/m", "m")

    def test_diff_format(self):
        with pytest.raises(ValueError):
            build_multiplier("g*m", "g/m")

    # --- Functions to test standard temp conversion functions --- #
    def test_celcius_conversions(self):
        assert celcius_convert(12, "C") == 12
        assert isclose(celcius_convert(43.2, "F"), 6.222, rel_tol=1e-3)
        assert isclose(celcius_convert(43.2, "K"), -229.95, rel_tol=1e-2)
        with pytest.raises(ValueError):
            celcius_convert(12, "impossible")

    def test_fahrenheit_conversions(self):
        assert fahrenheit_convert(12, "F") == 12
        assert isclose(fahrenheit_convert(43.2, "C"), 109.76, rel_tol=1e-2)
        assert isclose(fahrenheit_convert(43.2, "K"), -381.91, rel_tol=1e-2)
        with pytest.raises(ValueError):
            fahrenheit_convert(12, "impossible")

    def test_kelvin_conversions(self):
        assert kelvin_convert(12, "K") == 12
        assert isclose(kelvin_convert(43.2, "C"), 316.35, rel_tol=1e-2)
        assert isclose(kelvin_convert(43.2, "F"), 279.3722, rel_tol=1e-4)
        with pytest.raises(ValueError):
            kelvin_convert(12, "impossible")

    def test_invalid_temp(self):
        with pytest.raises(ValueError):
            tempConversion(12, "C", "impossible")

    
    # --- Functions to test unit multiplier functions --- #
    def test_C_convert(self):
        assert C_convert("F") == (5/9)
        assert C_convert("C") == 1
        with pytest.raises(ValueError):
            C_convert("impossible")

    def test_F_convert(self):
        assert F_convert("C") == (1 / (5/9))
        assert F_convert("F") == 1
        with pytest.raises(ValueError):
            F_convert("impossible")

    def test_mm_convert(self):
        assert mm_convert("mm") == 1
        assert mm_convert("m") == 1000
        assert mm_convert("km") == 1e6
        assert mm_convert("micron") == 0.001
        assert mm_convert("cm") == 10
        assert mm_convert("in") == 25.4
        assert mm_convert("ft") == 304.8
        assert mm_convert("yd") == 914.4
        assert mm_convert("mi") == 1.609e+6

        with pytest.raises(ValueError):
            mm_convert("impossible")
            
    def test_m_convert(self):
        assert m_convert("mm") == 0.001
        assert m_convert("m") == 1
        assert m_convert("km") == 1000
        assert m_convert("micron") == 1e-6
        assert m_convert("cm") == 0.01
        assert m_convert("in") == 0.0254
        assert m_convert("ft") == 0.3048
        assert m_convert("yd") == 0.9144
        assert m_convert("mi") == 1609.34

        with pytest.raises(ValueError):
            m_convert("impossible")

    def test_km_convert(self):
        assert km_convert("mm") == 1e-6
        assert km_convert("m") == 1e-3
        assert km_convert("km") == 1
        assert km_convert("micron") == 1e-9
        assert km_convert("cm") == 1e-5
        assert km_convert("in") == 2.54e-5
        assert km_convert("ft") == 0.0003048
        assert km_convert("yd") == 0.0009144
        assert km_convert("mi") == 1.60934

        with pytest.raises(ValueError):
            km_convert("impossible")
    
    def test_v_convert(self):

        assert v_convert("V") == 1
        assert v_convert("kV") == 1000
        assert v_convert("MV") == 1e6
        assert v_convert("mV") == 0.001

        with pytest.raises(ValueError):
            v_convert("impossible")

    def test_kg_convert(self):

        assert kg_convert("g") == 0.001
        assert kg_convert("kg") == 1
        assert kg_convert("mg") == 1e-6
        assert kg_convert("lb") == 0.453592
        assert kg_convert("oz") == 0.0283495

        with pytest.raises(ValueError):
            kg_convert("impossible")

    def test_g_convert(self):

        assert g_convert("g") == 1
        assert g_convert("kg") == 1000
        assert g_convert("mg") == 0.001
        assert g_convert("lb") == 453.592
        assert g_convert("oz") == 28.3495

        with pytest.raises(ValueError):
            g_convert("impossible")

    def test_bar_convert(self):

        assert bar_convert("bar") == 1
        assert bar_convert("MPa") == 10
        assert bar_convert("GPa") == 10000
        assert bar_convert("Pa") == 1e-5
        assert bar_convert("torr") == 0.00133322
        assert bar_convert("psi") == 0.0689476

        with pytest.raises(ValueError):
            bar_convert("impossible")

    def test_gpa_convert(self):

        assert gpa_convert("GPa") == 1
        assert gpa_convert("MPa") == 0.001
        assert gpa_convert("Pa") == 1e-9
        assert gpa_convert("bar") == 0.0001
        assert gpa_convert("psi") == 6.89476e-6
        assert gpa_convert("torr") == 1.3332e-7

        with pytest.raises(ValueError):
            gpa_convert("impossible")

    def test_ohm_convert(self):

        assert ohm_convert("ohm") == 1
        assert ohm_convert("Tohm") == 1e12
        assert ohm_convert("Gohm") == 1e9
        assert ohm_convert("Mohm") == 1e6
        assert ohm_convert("kohm") == 1e3
        assert ohm_convert("mohm") == 1e-3
        assert ohm_convert("microohm") == 1e-6
        assert ohm_convert("nohm") == 1e-9
        assert ohm_convert("pohm") == 1e-12

        with pytest.raises(ValueError):
            ohm_convert("impossible")

    def test_is_float(self):

        assert is_float("0.9")
        assert not is_float("impossible")

    