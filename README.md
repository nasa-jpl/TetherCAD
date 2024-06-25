
# TetherCAD README

This page gives a quick start on how to use TetherCAD. For more detailed usage and documentation please see the GitHub's wiki page!

# Installation

TetherCAD can be installed using the setup.py file included. Simply clone the repository, move to the repository folder and run:
`python3 -m setup.py install` to install the library directly on your system or within a virtual environment. 

# Quick Start

TetherCAD is best used within standalone Python scripts or Jupyter notebooks. The following subsections show examples of defining simple geometries and using them in analyses. The operation of the classes and functions used in these examples are explained in more detail during the Usage and Design sections of this wiki page. 

## Example Geometry
A piece of simple example code building a tether is shown below:

### Code:
```
from tethercad.tether_analysis.TetherDesign import Layer, Wire, RoundTetherDesign

# Construct a kapton insulated copper wire
conductor = Layer("copper", "Conductor", 0.25)
insulator = Wire("electrical", "kapton", "Kapton Wire", 0.15, innerLayer=conductor)

# Arrange the wires in a list
wireList = [insulator, insulator, insulator, insulator]

# Build the layers of the tether
core = Layer("fep", "Core", 0.25, color="black")
memLayer = Layer("fep", "Wires", 0.1, innerLayer=core, memberList=wireList)
strength_layer = Layer("vectran_kuraray_ht", "strength_layer", 0.15, innerLayer=memLayer)
abrasionLayer = Layer("ptfe", "abrasion", 0.2, innerLayer=strength_layer)

# Pass to a tether design object
tether = RoundTetherDesign("Example", abrasionLayer, 100)

tether.illustrate()
tether.tetherDetails()
```
### Output:

![image](https://github.jpl.nasa.gov/storage/user/10303/files/aec29926-7b1d-4345-9cf6-d9e9554fd085)

![Screenshot 2023-08-09 at 4 39 27 PM](https://github.jpl.nasa.gov/storage/user/10303/files/3e72ec52-3462-471a-98e1-f07c554145d2)


### Explanation:

This code starts by importing necessary classes and then by building a simple Kapton insulated wire. This wire has an outer layer Kapton which is 0.15mm thick, and an inner conductor layer of Copper which is 0.25mm thick. 

It then creates a list representing four of these wires for use within the tether design. 

Next the code sets up the core layer, the inner most layer of the tether. It is a layer made of FEP, with a thickness of 0.25mm, and has its color set to black to distinguish it from other FEP layers. 

A layer of FEP containing members (default colored) is then placed around the core layer, with a thickness of 0.1mm The list of wires are passed as an argument, and the layer is sized appropriately to fit them evenly spaced. The  thickness for this layer is increased from 0.1mm to a larger quantity to support the size of the wires. 

The last two layers: strength and abrasion, are placed around the existing layers following the same format to finish the cross sectional design. 

This design is then given as an input to a RoundTetherDesign object, which provides additional context such as the length, and computes the mass. Here we pass this tether object a name of "Example", the uppermost layer within our tree: `abrasion` and give it a length of 100m.

Using the new tether object, the tether in its entirety can be visualized and we can print out its outer diameter, radius and mass. 

## Example Analysis

A piece of simple example code running a power analysis is shown below. **Note:** This code uses the tether design shown in a previous example.  

### Code:

```
from tethercad.calculation_libraries.PowerAnalysis import dc_power_transmission_analysis

pathList = tether.findWires("electrical")
send_paths = [pathList[0], pathList[2]]
return_paths = [pathList[1], pathList[3]]

efficiency = dc_power_transmission_analysis(tether, 500, 1000, send_paths, return_paths, True)
```
### Output:

![Screenshot 2023-08-09 at 6 15 05 PM](https://github.jpl.nasa.gov/storage/user/10303/files/427ad8e3-6139-40a1-8de2-adf6adf7f011)

### Explanation

First, this code imports the DC power analysis function: dc_power_transmission_analysis.

It then identifies the "paths" of all electrical wires within the tethers. When a layer tree is given to a RoundTetherDesign object, it initializes a unique path string for each layer in the tree. This path is used to identify the layer by other code such as analyses. 

These paths are split into lists to indicate which wires should be used for the send and return paths of the transmission line. 

The dc_power_transmission_analysis is then called, passing tether object, a line input voltage of 500V, a desired power transmission of 1kW, the lists of send/return paths, and a boolean whether to output transmission information. 

# Running Tests

The tests for the system are built using pytest, and can be run from within the tethercad folder using the command:
`python3 -m pytest tests/<test_filename>.py

As of writing there are three test files: 
* test_analyses.py
* test_databases.py
* test_tether_design.py

# Copyright Information

Copyright 2024, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government 
Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the 
California Institute of Technology. This software may be subject to U.S. export control laws. 
By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. 
User has the responsibility to obtain export licenses, or other export authority as may be required before exporting 
such information to foreign countries or providing access to foreign persons."

If you have questions regarding this, please contact the JPL Software Release Authority at x4-2458.


