""""
Copyright 2024, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government 
Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the 
California Institute of Technology. This software may be subject to U.S. export control laws. 
By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. 
User has the responsibility to obtain export licenses, or other export authority as may be required before exporting 
such information to foreign countries or providing access to foreign persons."

If you have questions regarding this, please contact the JPL Software Release Authority at x4-2458.
"""

import os
import pathlib
import re
import pandas as pd
import warnings
import numpy as np


class DatabaseControl:
    """Database control class (singleton) that is used to manage databases of wires and materials./ 
    """

    _self = None
    _initialized = False

    # Only one instance is allowed across the entire process #
    def __new__(cls, **kwargs):
        if cls._self is None:
            cls._self = super().__new__(cls)
        
        return cls._self

    def __init__(self, path=None, reinit=False):
        """Database constructor, loads all of the databases if can find within the specified or default path

        Args:
            path (str, optional): Path in which to look for the material and wire databases. Defaults to None.
        """

        # Singleton initializtion logic        
        if not self._initialized:
            self._initialized = True
        else:
            if reinit == False:
              return
        
        self.wire_db_list = []
        self.wire_db_dict = {}

        self.conductor_db_list = []
        self.conductor_db_dict = {}

        self.optical_db_list = []
        self.optical_db_dict = {}

        # Temperature units #
        self.temp_units = ["C", "F"]

        # Default units for standard columns for wire, optical, and material databases #
        self.wireDefaultUnits = {"insulation_thickness":"mm", "insulation_od":"mm", "conductor_od":"mm",
                                 "weight/length":"kg/km", "Voltage(DC)":"V"}
        self.opticalDefaultUnits = {"core_od":"mm", "cladding_od":"mm", "jacket_od":"mm", "weight/length":"kg/km"}
        self.materialDefaultUnits = {"density":"g/mm^3", "strength":"GPa", "stress":"GPa", "elastic_modulus":"GPa", "thermal_cond":"W/(m*K)", "min_temp":"C", 
                                     "max_temp":"C", "dielectric_strength":"V/mm", "resistivity":"ohm*mm"}

        # Set up a dictionary of the expected property columns in the material database #
        # self.material_prop_list = ["density", "strength", "modulus", "brittle_temp", "max_temp", "dielectric_strength", 
        #                            "dielectric_constant", "resistivity", "permeability"]

        # Temp coefficients for common conductor materials #
        self.tempCoeffDict = {"silver":3.8e-3, "copper":3.9e-3, "gold":3.4e-3, "aluminum":3.9e-3, "tungsten":4.5e-3, 
                              "iron":5e-3, "nickel":6e-3}

        # Either search for databases in the databases folder package directory, or a user given path #
        
        if path == None:
            self.dbPath = pathlib.Path(__file__).resolve().parent
            print(self.dbPath)
        else:
            self.dbPath = pathlib.Path(path)

        # Load material database #
        self.material_db_list = []
        self.material_db_dict = {}

        # Look for valid databases in the database path #
        for files in os.listdir(self.dbPath):
            
            # Use the path to read the .csv into a pandas dataframe #
            if(files.endswith(".csv") and "_" in files):
                filename = re.findall(r'^(.+?)_', files)[0]
                location = self.dbPath / files
                print("Reading database %s" % filename)
                curr_database = pd.read_csv(location).replace(np.nan, None, regex=True)
            else:
                continue

            # Process as a wire database #
            if files.endswith('_ConductorDB.csv'):

                # Process Database
                updated_db = self.process_wire_database(curr_database)

                # Add database to internal wire record #
                self.conductor_db_list.append(filename)
                self.conductor_db_dict[filename] = updated_db

            # Process as a wire database #
            if files.endswith('_WireDB.csv') and 'Template' not in files:

                # Process Database as a wire database#
                updated_db = self.process_wire_database(curr_database)

                # Add database to internal wire record #
                self.wire_db_list.append(filename)
                self.wire_db_dict[filename] = updated_db

            # Process as an optical database #
            elif files.endswith('_OpticalDB.csv') and 'Template' not in files:

                # Process Database as an optical database #
                updated_db = self.process_wire_database(curr_database, wire_type="optical")

                # Add database to internal optical record #
                self.optical_db_list.append(filename)
                self.optical_db_dict[filename] = updated_db

            # Process as a material database #
            elif files.endswith('material_DB.csv'): 
                db = self.process_material_database(curr_database)
                
                # Replace any NaNs with None (Stuff can be weird when loading/accessing dataframes after initialization)
                db.replace(np.nan, None, inplace=True)
                
                self.material_db_list.append(filename)
                self.material_db_dict[filename] = db

            else:
                continue

    def load_from_folder(self, path):
        self.__init__(path=path)

    def process_material_database(self, db):
        """Processes the data from the created pandas dataframe, converting
        units and values to the appropriate defaults for the material database.

        Args:
            db (Pandas.DataFrame): The dataframe created when the material database .csv file is read into the system

        Returns:
            Pandas.DataFrame: A sliced dataframe containing properly converted values. 
        """
        
        cols = db.columns

        for col in cols:

            # Ignore the unitand name columns #
            if "_unit" in col or "_name" in col:
                continue

            # Grabs default unit for the column #
            if col not in self.materialDefaultUnits.keys():
                continue
            defaultUnit = self.materialDefaultUnits[col]
            
            if (unit_col := col + "_unit") in cols:
                
                # Temperature requires a different conversion method #
                if defaultUnit in self.temp_units: 
                    
                    for idx, val in enumerate(db[col]):
                        unit = db[unit_col][idx]

                        # Skip values where the unit is empty #
                        if unit is None:
                            continue

                        if(is_float(val)):
                            db.loc[idx, col] = tempConversion(val, unit, defaultUnit)
                            db.loc[idx, unit_col] = defaultUnit

                else:

                    for idx, val in enumerate(db[col]):
                        unit = db[unit_col][idx]

                        # Skip values where unit is empty, values should be empty as well #
                        if unit is None:
                            continue

                        multiplier = build_multiplier(unit, defaultUnit)

                        if(is_float(val)):
                            db.loc[idx, col] = float(val) * multiplier
                            db.loc[idx, unit_col] = defaultUnit

        return db.iloc[:, :-2]
            


        
    def process_wire_database(self, db, wire_type="wire"):
        """Processes wire databases to have the correct units and fill any missing
        information. 

        Args:
            db (pandas.DataFrame): Dataframe of the corresponding database read from its .csv
            wire_type (str, optional): The type of wire this represents. Defaults to "wire".

        Raises:
            ValueError: Invalid type of wire database specified
            ValueError: Missing units

        Returns:
            pandas.DataFrame: A sliced dataframe with all of its units properly converted
        """

        cols = db.columns

        defaultUnitList = None

        if wire_type == "wire":
            defaultUnitList = self.wireDefaultUnits
        elif wire_type == "optical":
            defaultUnitList = self.opticalDefaultUnits
        else:
            raise ValueError("Invalid type of database specified!")

        # Loop through columns
        for col in cols:

            if col not in defaultUnitList.keys():
                continue
            defaultUnit = defaultUnitList[col]

            # Grab the appropriate unit from the column, throwing an error if we expect one but none is present #
            unit = db[col][0]
            if unit is None:
                raise ValueError("Got None for unit when a default is present!")

            multiplier = build_multiplier(unit, defaultUnit)

            # Modify units using the multiplier #
            for idx, val in enumerate(db[col]):

                if(is_float(val)):
                    db.loc[idx, col] = float(val) * multiplier
                
                
        # Return a database that's been trimmed #
        return db.loc[1:]
    

    #TODO in some future work
    def build_override_material(self, material, entry):
        pass
    def add_override_material(self, material_entry):
        pass

    def find_material(self, material):
      for filename,db in self.material_db_dict.items():
          # Grab the desired material entry #
          foundEntry = searchEntry(db, {"material_name":material})  
          if not foundEntry.empty:
            break
      
      if foundEntry.empty:
        return None, None
      else:
        return foundEntry, filename
      
    def add_temp_material(self, material_entry):
      temp_name = 'temp_db'
      temp_db = self.material_db_dict.get(temp_name, None)
      
      if temp_db is None:
        temp_db = pd.DataFrame(material_entry)
        self.material_db_dict[temp_name] = temp_db
        self.material_db_list.append(temp_name)
      else:
        searchEntry(temp_db, {'material_name': material_entry['material_name']})
        pd.concat([temp_db, material_entry], ignore_index=True)
        

    def get_material_entry(self, material, fill_missing=True):
        """Returns a corresponding material entry from the material database. 

        Args:
            material (str): Name of the material within the material database
            fill_missing (bool, optional): Whether to fill missing information with a representative material. 
            Defaults to True.

        Raises:
            ValueError: Material not found in the material database
            ValueError: Invalid type specified for representative material
            ValueError: Material is missing representative material type
            ValueError: Material is missing a unit for one of its values

        Returns:
            pandas.DataFrame: DataFrame containing the values of the specified material, and any needed properties 
            filled from a representative material type. 
        """
        
        
        foundEntry, _ = self.find_material(material)
        
        if foundEntry is None:          
          raise ValueError("Passed material: %s not found in material databases!" % material)
        
        # Create a copy to avoid messing with the actual dataframe #
        entry = foundEntry.copy(deep=True)

        # Replace any values if necessary and issue a warning #
        if fill_missing:
          
            # Grab the fill entry based on the material type #
            material_type = entry["material_type"].item()
            # type_entry = searchEntry(self.material_db, {"material_name":material_type})
            type_entry, _ = self.find_material(material_type)
  
            if type_entry.empty:
                raise ValueError("Invalid material type of %s specified for %s!" % (material_type, material))
            
            replaced_list = []
            for col in entry:
  
                # Ignore unit columns #
                if "_unit" in col:
                    continue
  
                # If we don't have a value for this portion #
                if entry[col].iloc[0] is None:
                    
                    # Add it to the list of missing types #
                    replaced_list.append(col)
                    unit_col = col + "_unit"
  
                    reprEntryUnit = None                
  
                    # Double check the material entry isn't missing this value #
                    if (reprEntry := type_entry[col].item()) is None:
                        raise ValueError("Material type: %s is missing a value for: %s" % (material_type, col))
  
                    # Only add a unit if we have a unit column for it #
                    if unit_col in self.material_db_dict[self.material_db_list[0]].columns:
                        if (reprEntryUnit := type_entry[unit_col].item()) is None:
                            raise ValueError("Material type: %s is missing a unit for: %s" % (material_type, unit_col))
  
                    # Set the entry to the representative value #
                    entry[col] = reprEntry
  
                    # As long as we have a unit to update, add it #
                    if reprEntryUnit is not None:
                        entry[unit_col] = reprEntryUnit

            if len(replaced_list) > 0:
                repl_str = ", ".join(replaced_list)

                # NOTE: This process is a bit of a bandaid fix, continuation wasn't playing nice with warnings #
                warningStr1 = "The following material properties for %s" % material
                warningStr2 = " were pulled from a representative material:" 
                warningStr3 = "(%s): %s" % (material_type, repl_str)
                warnings.warn(warningStr1 + warningStr2 + warningStr3)
        else:
            warnings.warn("fill_missing is False, entry may be missing values!")
                
        return entry



def searchEntry(database, queryDict):
    """Searches for entries in the passed database that match parameters specified in the key/value pairs of query dict

    Args:
        database (pandas.Dataframe): Database object which to search
        queryDict (dictionary): A set of key value pairs representing column/value used to search the database
    
    Returns:
        pandas.dataframe: Dataframe of entry or entries that match the specified criteria in the queryDict 
    """

    database.replace(np.nan, None, inplace=True)

    if not queryDict:
        raise ValueError("Empty query dictionary passed!")

    query = ""
    for idx, (key, val) in enumerate(queryDict.items()):
        query += "(database[\"%s\"] == %r)" % (key, val)
        if query != "" and idx != len(queryDict.items()) - 1:
            query += " & "
    return database.loc[eval(query)]

def get_material_property(material_entry, property):
    """Get a specific material property from a material entry. 
    Used for readability. 

    Args:
        material_entry (pandas.DataFrame): DataFrame containing material property information. 
        property (str): Property to grab from the DataFrame

    Returns:
        object: Material property (can have a number of different possible types)
    """
    return material_entry[property].item()
    
def replace_with_mult(unit, defaultUnit, mult_str):
    """Replace units in a given string to their corresponding multiplier values to the 
    specified default units. 

    Args:
        unit (str): Individual unit portion from the string
        defaultUnit (str): Default unit to convert the individual to
        mult_str (str): Entire unit string being worked on

    Returns:
        str: The string with instances of unit replaced with the multiplier to convert it to default unit
    """
    conv = get_conv_func(defaultUnit)
    unit_multiplier = conv(unit)
    return mult_str.replace(unit, str(unit_multiplier))

def build_multiplier(unit_str, default_str):
    """Builds a multiplier to convert the unit string to the default. 
    For example: g/ohm * cm can be converted to kg/ohm * in.

    Args:
        unit_str (str): Unit string to convert from
        default_str (str): Unit string to convert to

    Raises:
        ValueError: Unit string has a different number of units
        ValueError: Unit string has a different format

    Returns:
        float: Multiplier to use to multiply the passed unit string to the default unit string
    """

    # Pick out each individual unit from the unit str #
    mult_str = unit_str.replace("^", "**")
    reg_str = "[A-Za-z]+"
    unitList = re.findall(reg_str, unit_str)
    defaultList = re.findall(reg_str, default_str)

    # Check the # of units are the same #
    if len(unitList) != len(defaultList):
        raise ValueError("The # of individual units in the passed string and default string do not match!")
    
    # Check the format is the same #
    unitTestStr = re.sub(reg_str, "-", unit_str)
    defTestStr = re.sub(reg_str, "-", default_str)
    if unitTestStr != defTestStr:
        raise ValueError("The format of %s does not match %s!" % (unit_str, default_str))

    for idx, unit in enumerate(unitList):
        if unit == "m":
            continue # Do meters last to avoid messing things up
        elif unit == "V":
            continue
        elif unit == "Pa":
            continue
        else:
            defaultUnit = defaultList[idx]
            mult_str = replace_with_mult(unit, defaultUnit, mult_str)

    if "m" in unitList:
        mult_str = replace_with_mult("m", defaultList[unitList.index("m")], mult_str)
    if "V" in unitList:
        mult_str = replace_with_mult("V", defaultList[unitList.index("V")], mult_str)
    if "Pa" in unitList:
        mult_str = replace_with_mult("Pa", defaultList[unitList.index("Pa")], mult_str)   

    multiplier = eval(mult_str)

    return multiplier

def tempConversion(temp, unit, defaultUnit):
    """Individual temp conversion function

    Args:
        temp (float): Temperature value
        unit (str): Unit for the temperature value
        defaultUnit (str): Unit to convert the temperature into

    Raises:
        ValueError: Invalid temperature unit passed

    Returns:
        float: Converted temperature
    """
    if defaultUnit == "C":
        return celcius_convert(temp, unit)
    elif defaultUnit == "F":
        return fahrenheit_convert(temp, unit)
    elif defaultUnit == "K":
        return kelvin_convert(temp, unit)
    else:
        raise ValueError("Invalid default temperature unit passed: %s! Valid temp units are C, F, K" % defaultUnit)

def celcius_convert(temp, unit):
    """Conversion function from celcius to other units

    Args:
        temp (float): Temperature in celcius
        unit (str): Unit to convert to

    Raises:
        ValueError: Unrecognized unit passed

    Returns:
        float: Converted Temperature
    """

    newTemp = 0

    if unit == "C":
        newTemp = temp
    elif unit == "F":
        newTemp = (temp - 32) * (5/9)
    elif unit == "K":
        newTemp = (temp - 273.15)
    else:
        raise ValueError("Unrecognized unit passed: %s! Valid units are: C, F, K" % unit)
    
    return newTemp

def fahrenheit_convert(temp, unit):
    """Conversion function from fahrenheit to other units

    Args:
        temp (float): Temperature in fahrenheit
        unit (str): Unit to convert to

    Raises:
        ValueError: Unrecognized unit passed

    Returns:
        float: Converted Temperature
    """
    newTemp = 0

    if unit == "F":
        newTemp = temp
    elif unit == "C":
        newTemp = (temp * (9/5)) + 32
    elif unit == "K":
        newTemp = ((temp - 273.15) * (9/5)) + 32
    else:
        raise ValueError("Unrecognized unit passed: %s! Valid units are: C, F, K" % unit)
    
    return newTemp

def kelvin_convert(temp, unit):
    """Conversion function from kelvin to other units

    Args:
        temp (float): Temperature in kelvin
        unit (str): Unit to convert to

    Raises:
        ValueError: Unrecognized unit passed

    Returns:
        float: Converted Temperature
    """
    newTemp = 0

    if unit == "K":
        newTemp = temp
    elif unit == "C":
        newTemp = temp + 273.15
    elif unit == "F":
        newTemp = (((temp - 32) * (5/9)) + 273.15)
    else:
        raise ValueError("Unrecognized unit passed: %s! Valid units are: C, F, K" % unit)

    return newTemp

def get_conv_func(defaultUnit):
    """Returns the function used to convert a unit to the passed default unit

    Args:
        defaultUnit (str): Unit we want to convert to

    Raises:
        ValueError: Default unit is not found (no implementation of a conversion function for it)

    Returns:
        function: Function to use to convert 
    """
    
    func = None
    if defaultUnit == "mm":
        func = mm_convert
    elif defaultUnit == "m":
        func = m_convert
    elif defaultUnit == "km":
        func = km_convert
    elif defaultUnit == "V":
        func = v_convert
    elif defaultUnit == "kg":
        func = kg_convert
    elif defaultUnit == 'g':
        func = g_convert
    elif defaultUnit == "GPa":
        func = gpa_convert
    elif defaultUnit == "ohm":
        func = ohm_convert
    elif defaultUnit == "C":
        func = C_convert
    elif defaultUnit == "F":
        func = F_convert
    elif defaultUnit == "K":
        func = K_convert
    elif defaultUnit == "bar":
        func = bar_convert
    elif defaultUnit == "W":
        func = w_convert
    else:
        raise ValueError("Invalid default unit specified: %s !" % defaultUnit)
    
    return func

def C_convert(unit):
    """Converts to celcius from a unit when part of a larger unit representation (no 32 offset)

    Args:
        unit (str): Unit to conver to

    Raises:
        ValueError: Unrecognized unit was passed

    Returns:
        float: multiplier to convert unit to celcius
    """
    multiplier = -1

    # If we are using this function, then C is part of a composite and can use a multiplier #

    if (unit == "F"):
        multiplier = 5/9
    elif(unit == "C"):
        multiplier = 1
    else:
        raise ValueError("Unrecognized unit passed: %s! Valid units are: C, F" % unit)
    
    return multiplier

def F_convert(unit):
    """Converts to fahrenheit from a unit when part of a larger unit representation (no 32 offset)

    Args:
        unit (str): Unit to conver to

    Raises:
        ValueError: Unrecognized unit was passed

    Returns:
        float: multiplier to convert unit to fahrenheit
    """

    multiplier = -1

    # If we are using this function, then C is part of a composite and can use a multiplier #

    if (unit == "C"):
        multiplier = 1 / (5/9)
    elif(unit == "F"):
        multiplier = 1
    else:
        raise ValueError("Unrecognized unit passed: %s! Valid units are: C, F" % unit)
    
    return multiplier

def K_convert(unit):
    if(unit == "K"):
        multiplier = 1
    else:
        ### NOTE: This is a bandaid fix for now
        raise ValueError("Unrecognized unit passed: %s! Valid units are: K" % unit)
    
    return multiplier

def mm_convert(unit):
    """Function to return multiplier from a passed unit to mm

    Args:
        unit (str): Unit to convert from

    Raises:
        ValueError: Unit to convert from has no entry

    Returns:
        float: Multiplier to convert from unit to mm
    """
    multiplier = -1

    # Check unit and return multiplier to convert it into millimeters #
    if(unit == 'mm'):
        multiplier = 1
    elif(unit == 'micron'):
        multiplier = 0.001
    elif(unit == 'cm'):
        multiplier = 10
    elif(unit == 'm'):
        multiplier = 1000
    elif(unit == 'km'):
        multiplier = 1000000
    elif(unit == 'in'):
        multiplier = 25.4
    elif(unit == 'ft'):
        multiplier = 304.8
    elif(unit == 'yd'):
        multiplier = 914.4
    elif(unit == 'mi'):
        multiplier = 1.609e6
    else:
        raise ValueError("Unrecognized unit passed: %s! Valid units are: mm, cm, m, km, in, ft, yd" % unit)

    return multiplier


def m_convert(unit):
    """Function to return multiplier from a passed unit to m

    Args:
        unit (str): Unit to convert from

    Raises:
        ValueError: Unit to convert from has no entry

    Returns:
        float: Multiplier to convert from unit to m
    """
    multiplier = -1

    # Check unit and return multiplier to convert it into meters #
    if(unit == 'mm'):
        multiplier = 0.001
    elif(unit == 'cm'):
        multiplier = 0.01
    elif(unit == 'micron'):
        multiplier = 1e-6
    elif(unit == 'm'):
        multiplier = 1
    elif(unit == 'km'):
        multiplier = 1000
    elif(unit == 'in'):
        multiplier = 0.0254
    elif(unit == 'ft'):
        multiplier = 0.3048
    elif(unit == 'yd'):
        multiplier = 0.9144
    elif(unit == 'mi'):
        multiplier = 1609.34
    else:
       raise ValueError("Unrecognized unit passed %s! Valid units are: mm, cm, m, km, in, ft, yd, mi" % unit)

    return multiplier

def km_convert(unit):
    """Function to return multiplier from a passed unit to km

    Args:
        unit (str): Unit to convert from

    Raises:
        ValueError: Unit to convert from has no entry

    Returns:
        float: Multiplier to convert from unit to km
    """
    multiplier = -1

    if(unit == 'mm'):
        multiplier = 1e-6
    elif(unit == 'micron'):
        multiplier = 1e-9
    elif(unit == 'cm'):
        multiplier = 1e-5
    elif(unit == 'm'):
        multiplier = 1e-3
    elif(unit == 'km'):
        multiplier = 1
    elif(unit == 'ft'):
        multiplier = 0.0003048
    elif(unit == 'yd'):
        multiplier = 0.0009144
    elif(unit == 'mi'):
        multiplier = 1.60934
    elif(unit == 'in'):
        multiplier = 2.54e-5
    else:
       raise ValueError("Unrecognized unit passed %s! Valid units are: mm, cm, m, km, in, ft, yd, mi" % unit)

    return multiplier
    

def v_convert(unit):
    """Function to return multiplier from a passed unit to V

    Args:
        unit (str): Unit to convert from

    Raises:
        ValueError: Unit to convert from has no entry

    Returns:
        float: Multiplier to convert from unit to V
    """
    multiplier = -1

    if(unit == 'V'):
        multiplier = 1
    elif(unit == 'mV'):
        multiplier = 0.001
    elif(unit == 'kV'):
        multiplier = 1000
    elif(unit == 'MV'):
        multiplier = 1e6
    else:
        raise ValueError("Unrecognized unit passed: %s ! Valid units are: V, mV, kV" % unit)

    return multiplier

def w_convert(unit):
    """Function to return multiplier from a passed unit to W

    Args:
        unit (str): Unit to convert from

    Raises:
        ValueError: Unit to convert from has no entry

    Returns:
        float: Multiplier to convert from unit to W
    """

    multiplier = -1

    if(unit == 'W'):
        multiplier = 1
    elif(unit == 'mW'):
        multiplier = 0.001
    elif(unit == 'kW'):
        multiplier = 1000
    elif(unit == 'MW'):
        multiplier = 1e6
    else:
        raise ValueError("Unrecognized unit passed: %s ! Valid units are: W, mW, kW, MW" % unit)
    
    return multiplier

def kg_convert(unit):
    """Function to return multiplier from a passed unit to kg

    Args:
        unit (str): Unit to convert from

    Raises:
        ValueError: Unit to convert from has no entry

    Returns:
        float: Multiplier to convert from unit to kg
    """
    multiplier = -1

    # Check unit and return multiplier to convert it to grams #
    if(unit == 'g'):
        multiplier = 1e-3
    elif(unit == 'mg'):
        multiplier = 1e-6
    elif(unit == 'kg'):
        multiplier = 1
    elif(unit == 'lb'):
        multiplier = 0.453592
    elif(unit == 'oz'):
        multiplier = 0.0283495
    else:
        raise ValueError("Unrecognized unit passed %s! Valid units are: mg, g, kg, lb, oz" % unit)

    return multiplier

def g_convert(unit):
    """Function to return multiplier from a passed unit to g

    Args:
        unit (str): Unit to convert from

    Raises:
        ValueError: Unit to convert from has no entry

    Returns:
        float: Multiplier to convert from unit to g
    """

    multiplier = -1

    if(unit == 'g'):
        multiplier = 1
    elif(unit == 'mg'):
        multiplier = 1e-3
    elif(unit == 'kg'):
        multiplier = 1e3
    elif(unit == 'lb'):
        multiplier = 453.592
    elif(unit == 'oz'):
        multiplier = 28.3495
    else:
        raise ValueError("Unrecognized unit passed %s! Valid units are mg, g, kg, lb, oz" % unit)
    
    return multiplier

def bar_convert(unit):
    """Function to return multiplier from a passed unit to bar

    Args:
        unit (str): Unit to convert from

    Raises:
        ValueError: Unit to convert from has no entry

    Returns:
        float: Multiplier to convert from unit to bar
    """
    multiplier = -1

    if(unit == 'bar'):
        multiplier = 1
    elif(unit == 'MPa'):
        multiplier = 10
    elif(unit == 'GPa'):
        multiplier = 10000
    elif(unit == 'Pa'):
        multiplier = 1e-5
    elif(unit == 'psi'):
        multiplier = 0.0689476
    elif(unit == 'torr'):
        multiplier = 0.00133322
    else:
        raise ValueError("Unrecgonized unit passed %s! Valid units are GPa, MPa, Pa, bar, psi, torr" % unit)
    
    return multiplier

def gpa_convert(unit):
    """Function to return multiplier from a passed unit to GPa

    Args:
        unit (str): Unit to convert from

    Raises:
        ValueError: Unit to convert from has no entry

    Returns:
        float: Multiplier to convert from unit to GPa
    """
    multiplier = -1

    if(unit == 'GPa'):
        multiplier = 1
    elif(unit == 'MPa'):
        multiplier = 1e-3
    elif(unit == 'Pa'):
        multiplier = 1e-9
    elif(unit == 'bar'):
        multiplier = 1e-4
    elif(unit == 'psi'):
        multiplier = 6.89476e-6
    elif(unit == 'torr'):
        multiplier = 1.3332e-7
    else:
        raise ValueError("Unrecgonized unit passed %s! Valid units are GPa, MPa, Pa, bar, psi, torr" % unit)

    return multiplier

def ohm_convert(unit):
    """Function to return multiplier from a passed unit to Ohm

    Args:
        unit (str): Unit to convert from

    Raises:
        ValueError: Unit to convert from has no entry

    Returns:
        float: Multiplier to convert from unit to Ohm
    """
    multiplier = -1

    if (unit == 'ohm'):
        multiplier = 1
    elif (unit == 'Tohm'):
        multiplier = 1e12
    elif (unit == 'Gohm'):
        multiplier = 1e9
    elif (unit == 'Mohm'):
        multiplier = 1e6
    elif (unit == 'kohm'):
        multiplier = 1e3
    elif (unit == 'mohm'):
        multiplier = 1e-3
    elif (unit == 'microohm'):
        multiplier = 1e-6
    elif (unit == 'nohm'):
        multiplier = 1e-9
    elif (unit == 'pohm'):
        multiplier = 1e-12
    else:
        raise ValueError("Unrecognized unit passed %s! Valid units are ohm, Tohm, Gohm, Mohm, kohm, mohm, \
                          microohm, nohm, pohm" % unit)
    
    return multiplier

def is_float(string):
    """Returns whether the passed string can be interpreted as a float, capturing 
    the error that can occur from trying directly. 

    Args:
        string (str): String to check

    Returns:
        bool: Whether the string can be cast as a float
    """
    try:
        float(string)
        return True
    except (TypeError, ValueError) as error:
        return False