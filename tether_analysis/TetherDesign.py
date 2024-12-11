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
import matplotlib
from matplotlib import pyplot as plt
from copy import deepcopy
import warnings
import re


import databases.DatabaseControl as DB

databases = DB.DatabaseControl()
        

class Layer():

    def __init__(self, material, name, thickness, innerLayer=None, memberList=None, copy=True, color=None, 
                 helixAngle=90, fillRatio=1.0, iter_limit=100000, step=0.001) -> None:
        """Constructor for Layer class.

        Args:
            material (str): The material name for this layer
            name (str): Name of this layer. Cannot be "".
            thickness (float): The thickness of this layer, specified in millimeters
            innerLayer (Layer, optional): Inner layer wrapped by this layer. Defaults to None.
            memberList (Layer, optional): A list of layers (which may contain others) 
            that should be helixed within this layer. Defaults to None.
            copy (bool, optional): Whether to deepcopy nested or helixed layers. 
            Will throw a warning if disabled. Defaults to True.
            color (str, optional): Desired color of the layer, overriding the material color. 
            All CSS named colors are valid. Defaults to None.
            helixAngle(float, optional): Angle at which the fill and any members in the member list are helixed. Defaults to 0.
            iter_limit (int, optional): Iteration limit of the fitting algorithm, 
            to avoid an infinite loop. Defaults to 100000.
            step (float, optional): The amount to step the layer size by when defining layer geometry. Defaults to 1e-3.

        Raises:
            ValueError: Must have a thickness of greater than 0
            ValueError: A name must be specified for this layer
            ValueError: Fill ratio must be between 0 and 1
        """

        if thickness <= 0:
            raise ValueError("Thickness must be a positive nonzero value!")
        if name is None or name == "":
            raise ValueError("A name must be specified for this layer!")

        self.innerRadius = 0                # (mm)
        self.memberStartRadius = 0          # (mm)
        self.outerRadius = 0                # (mm)
        self.innerDiameter = 0              # (mm)
        self.outerDiameter = 0              # (mm)
        self.layerMass = 0                  # (kg)
        self.layerThickness = thickness     # (mm)
        self.helixAngle = helixAngle        # (degrees)
        self.length = 0                     # (m) Used by RoundTetherDesign
        self.layerID = 0                    
        self.x = 0
        self.y = 0
        self.polarCoords = [0, 0]


        if fillRatio <= 0 or fillRatio > 1:
            raise ValueError("Fill ratio must be between (0, 1]!")
        
        self.fillRatio = fillRatio

        # Material of layer #
        self.layerMaterial = material

        # A name of the layer, can be used to discern in the details printout #
        self.name = name
        self.layerPath = ""

        self.color = color

        # Warn that this can be dangerous if copy is disabled #
        if not copy:
            warnings.warn("Notice: disabling copy can have unintended consequences on the overall design")

        # Initialize member list to be empty and inner layer to be none #
        self.memberList = []
        self.memberPosition = []

        # Set the inner layer #
        if copy:
            self.innerLayer = deepcopy(innerLayer)
        else:
            self.innerLayer = innerLayer

        # Do some member specific actions if any were passed #
        if(memberList is not None):

            if copy:
                for member in memberList:
                    self.memberList.append(deepcopy(member))
            else:
                self.memberList = memberList

            # Make size layer thickness to be enough to cover the largest member #
            self.minThickness = 0
            for member in self.memberList:

                if 2 * member.outerRadius > self.minThickness:
                    self.minThickness = member.outerDiameter

            if self.minThickness > self.layerThickness:
                self.layerThickness = self.minThickness

        # Initialize control for spacing members within the layer #
        self.equalSpacing = True
        self.equalSpacingAngle = 0
        self.extraAngle = 0
        self.nonNestedCount = 0
        self.iter_limit = iter_limit

        # Define our geometry on initialization #
        self.define_layer_geometry(step=step)

        self.check_layer_collisions()


    def _setLayerDimensions(self):
        """Define the dimensions of the layer from the thickness and inner radius. 
        """
        if self.innerLayer is not None:
            self.innerRadius = self.innerLayer.outerRadius
            self.memberStartRadius = self.innerRadius

        self.outerRadius = self.innerRadius + self.layerThickness
        self.innerDiameter = 2 * self.innerRadius
        self.outerDiameter = 2 * self.outerRadius


    def define_layer_geometry(self, step=0.001):
        """Defines the geometry of the layer, sizing things to be appropriate if necessary. 

        Args:
            step (float, optional): The step size used to enlarge the layer if necessary (mm). Defaults to 0.001.
        """

        self._setLayerDimensions()

        # Make sure we have layers in the layer list #
        if (len(self.memberList) < 1):
            return

        space = self._scale_start_radius(step)
       
        # If desired, set the extra angle allocation for each member #
        self.extraAngle = (space / len(self.memberList))

        # Add any extra space, and then reassign cartesian coordinates for matplotlib #
        for member in self.memberList:

                member.polarCoords[1] += self.extraAngle * self.memberList.index(member)
                coords = polar_to_cart(member.polarCoords)
                member.x = coords[0]
                member.y = coords[1]

        maxDist = self._solve_collisions()

        self._resize_layer(maxDist)

        for member in self.memberList:
            memPos = [member.x, member.y]
            self.memberPosition.append(memPos)

        # Update position of sub members based on this member's coords (important for helixed Layers)
        self.update_member_positions()


    def _scale_start_radius(self, step):
        """Scales the logical start radius for members until they can all be wrapped radially in within this layer

        Args:
            step (float): The size to increase the radius by until members fit (mm)

        Raises:
            RuntimeError: Whether the iteration limit has been hit. (To prevent long runtimes and infinite loops)

        Returns:
            float: left over radial space after fitting
        """
        # Initialize some starter variables #
        fit = False
        iter = 0

        # Keep increasing the inner radius until the members will fit #
        while (not fit):

            # Raise error if this takes > iter_limit iterations #
            if iter  == self.iter_limit:
                raise RuntimeError("Hit iteration limit! Try increasing the iteration limit or the step size!")
            else:
                iter += 1

            # Available space is 2 rad
            space = 2 * np.pi
            currAngle = 0

            # Calculate the necessary theta for each member and update the max if appropriate #
            for member in self.memberList:

                if self.memberList.index(member) == 0:
                    # Calculate the XY position for the centerpoint of the member, and rotate by our angle #
                    member.polarCoords = [self.memberStartRadius + member.outerRadius, currAngle]
                    self._assign_cart_from_polar(member)
                    continue

                prevMember = self.memberList[self.memberList.index(member) - 1]

                # Calculate sides of triangle #
                sideA = self.memberStartRadius + member.outerRadius
                sideB = self.memberStartRadius + prevMember.outerRadius
                sideC = member.outerRadius + prevMember.outerRadius

                # Use the triangle to create an appropriate angle allocation for this member #
                theta = get_angle_C(sideA, sideB, sideC)

                # Calculate the cartesian x,y based on polar coords #
                currAngle += theta
                member.polarCoords = [sideA, currAngle]
                self._assign_cart_from_polar(member)

                # Remove the allocated space from what is available #
                space -= theta

                # If we are the end, also perform the same process a second time, with the first member 
                if self.memberList.index(member) == len(self.memberList) - 1:
                    firstMember = self.memberList[0]
                    sideA = self.memberStartRadius + firstMember.outerRadius
                    sideB = self.memberStartRadius + member.outerRadius
                    sideC = firstMember.outerRadius + member.outerRadius
                    space -= get_angle_C(sideA, sideB, sideC)

            # As long as we haven't used more than 2 rad, everything fits #
            if space >= 0:
                fit = True
            # Otherwise, increment the starting radius (thickness of the layer and check again) #
            else:
                self.memberStartRadius += step

        return space
    

    def _solve_collisions(self):
        """Solves member/member collisions within a layer by scaling their X/Y position within the layer

        Returns:
            float: The max radial distance occupied by any member (radius from center + member radius)
            Used to resize the layer thickness if needed. 
        """

        # Calculate the max scale factor between colliding members #
        maxScale = 1
        for i in range(0, len(self.memberList)):
            for j in range(0, len(self.memberList)):
                
                # Ignore self to self comparisons #
                if i == j:
                    continue
                
                # Grab the positions of each member to each other, checking if they are at least tangent #
                mem1 = self.memberList[i]
                mem2 = self.memberList[j]
                p1 = [mem1.x, mem1.y]
                p2 = [mem2.x, mem2.y]
                if check_tangent(p1, p2, mem1.outerRadius, mem2.outerRadius):
                    continue

                # Determine the scale factor needed to make them at least tangent #
                actDist = np.linalg.norm(np.array(p1) - np.array(p2))
                reqDist = mem1.outerRadius + mem2.outerRadius
                scaleFactor = reqDist / actDist

                # Update our max scale factor from each member #
                if scaleFactor > maxScale:
                    maxScale = scaleFactor

        # Apply scale factor to each new member #
        maxDist = 0
        for member in self.memberList:
            member.polarCoords[0] *= maxScale
            coords = polar_to_cart(member.polarCoords)
            member.x = coords[0]
            member.y = coords[1]

            dist = member.polarCoords[0] + member.outerRadius
            if dist > maxDist:
                maxDist = dist

        return maxDist
    
    def _resize_layer(self, maxDist):
        """Resizes the layer based on the max distance occupied by a member

        Args:
            maxDist (float): The max radial distance occupied by any member (radius from center + member radius)
        """

        # Calculate our outer radius #
        self.outerRadius = self.innerRadius + self.layerThickness

        # Readjust layer thickness for the largest member if necessary #
        if self.outerRadius < maxDist:
            self.outerRadius = maxDist
            self.layerThickness = self.outerRadius - self.innerRadius
            self.outerDiameter = 2 * self.outerRadius
            self.innerDiameter = 2 * self.innerRadius


    def update_member_positions(self):
        """Updates any members within any internal Layers with this Layer's coordinates, called by 
        the encapsulating tether design
        """

        # Update the positions of helixed members, and then any layers they encapsulate #
        for idx, member in enumerate(self.memberList):

            memPos = self.memberPosition[idx]
            member.x = self.x + memPos[0]
            member.y = self.y + memPos[1]
            member.polarCoords = cart_to_polar([member.x, member.y])
            member.update_member_positions()

        # Update the position of any layer this one encapsulates, and any layers it encapsulates #
        if self.innerLayer is not None:
            self.innerLayer.x = self.x
            self.innerLayer.y = self.y
            self.innerLayer.polarCoords = cart_to_polar([self.innerLayer.x, self.innerLayer.y]) 
            self.innerLayer.update_member_positions()

 
    def check_member_collisions(self, verbose=False):
        """Returns whether there are any collisions between members within the same layer

        Args:
            verbose (bool, optional): Whether to print member positions which collide. Defaults to False.

        Returns:
            bool: Whether any collisions were found
        """
        collisionCount = 0
        for mem1 in self.memberList:
            for mem2 in self.memberList:

                if mem1 is mem2:
                    continue
                else:
                    if not check_tangent([mem1.x, mem1.y], [mem2.x, mem2.y], mem1.outerRadius, mem2.outerRadius):
                        collisionCount += 1

                        if verbose:
                            print("Member at: (%f, %f) with radius: %f collides with member at (%f, %f) with radius: %f" 
                                  % (mem1.x, mem1.y, mem1.outerRadius, mem2.x, mem2.y, mem2.outerRadius))
            
        if collisionCount >= 1:
            if verbose:
                print("%d Collisions detected!" % collisionCount)
            return True
        
        return False
    
    def check_layer_collisions(self, output_on_error=False, verbose=False):
        """Returns whether any members within this layer occupies space outside of the layer's outer radius. 
        (Does not check inner radius)

        Args:
            output_on_error (bool, optional): Whether to Draw the tether if a collision is detected. Defaults to False. 
            verbose (bool, optional): Whether to print specific collision information. Defaults to False.

        Raises:
            ValueError: Whether a layer/member collision is detected
        """

        collisionsPresent = False

        for mem in self.memberList:

            outer = mem
            inner = mem.innerLayer

            # Check once here 
            if(checkOverlap(self, outer)):
                collisionsPresent = True
                if verbose: print("Collision between member at: %s and %s detected" % (self.name, outer.name))

            while(inner is not None):

                if(checkOverlap(self, outer)):
                    collisionsPresent = True
                    if verbose: print("Collision between member at: %s and %s detected" % (self.name, outer.name))
                if(checkOverlap(self, inner)):
                    collisionsPresent = True
                    if verbose: print("Collision between member at: %s and %s detected" % (self.name, inner.name))
                if(checkOverlap(outer, inner)):
                    collisionsPresent = True
                    if verbose: print("Collision between member at: %s and %s detected" % (outer.name, inner.name))

                outer = inner
                inner = outer.innerLayer

            
            if collisionsPresent:
                if output_on_error:
                    self.illustrate()
                raise ValueError("A collision between two layers was detected!") 
                

    def _assign_cart_from_polar(self, member):
        """Protected function that sets a member's cartesian coordinates from its polar coordinates

        Args:
            member (_type_): _description_
        """
        coords = polar_to_cart(member.polarCoords)
        member.x = coords[0]
        member.y = coords[1]


    def layerDetails(self, recursive=False, tabNum=0, output=True, breakdownList=None):
        """Prints the details of this layer, calling it for the next layer if desired. 

        Args:
            recursive (bool, optional): Whether to print the details of the layer below. Defaults to False.
            tabNum (int, optional): Number of tabs to indent the details of this layer. Defaults to 0.
        """
        if output:
            startStr = ''.join("   " for x in range(0, tabNum))
            print(startStr + "Name: %s, Path (if set): %s" % (self.name, self.layerPath))
            print(startStr + "Layer Dimensions:")
            print(startStr + " - Inner Radius: %f mm" % self.innerRadius)
            print(startStr + " - Outer Radius: %f mm" % self.outerRadius)
            print(startStr + " - Layer Thickness: %f mm" % self.layerThickness)
            print(startStr + "Layer Length: %f m (This will be zero until assigned by a tether design object)" % self.length)
            print(startStr + "Material: %s" % self.layerMaterial)
            if(len(self.memberList) > 0):
                print(startStr + " - Equal Spacing: " + str(self.equalSpacing))
                print(startStr + " - Member Number: %d" % len(self.memberList))
                print(startStr + " - Member Start Radius: %f" %self.memberStartRadius)
        
        if breakdownList is not None:
            breakdownList.append([self.name, self.layerPath, self.layerMaterial, self.innerRadius, self.outerRadius, self.layerThickness, (self.x, self.y), self.length])

        if recursive and self.innerLayer is not None:
            self.innerLayer.layerDetails(recursive=True, tabNum=tabNum+2, output=output, breakdownList=breakdownList)

        if recursive:
            for member in self.memberList:
                member.layerDetails(recursive=True, tabNum=tabNum+2, output=output, breakdownList=breakdownList)    

    def illustrate(self, ax=None, borderless=False, figdim=5):
        """Illustrates this layer using matplotlib

        Args:
            ax (matplotlib.Axes, optional): Passed axes object on which to illustrate. Defaults to None.
            borderless (bool, optional): Whether to print a border around each layer. Defaults to False.
            figdim (int, optional): Size of the illustrated image in inches. Defaults to 5.
        """

        showAfter = False

        # If this is the entry point for the illustrate call, creat a figure and axis object # 
        if ax is None:
            fig, ax = plt.subplots(num=None, figsize=(figdim, figdim), dpi=150)
            showAfter = True

        # Set our color from the material property #
        if self.color is None:
            material_entry = DB.searchEntry(databases.material_db, {'material_name':self.layerMaterial.lower()})
            if material_entry.empty:
                warnings.warn("Unable to find specified material: %s! Using white as the color!" % 
                              self.layerMaterial.lower())
                self.color = "white"
            else:
                self.color = material_entry["color"].item()

        # Draw this layer #
        width = 1
        if (borderless):
            width = 0

        circle = matplotlib.patches.Circle(
            (self.x, self.y), radius=self.outerRadius, facecolor=self.color, edgecolor="black", linewidth=width)
        ax.add_patch(circle)

        # Draw any helixed members #
        for member in self.memberList:
            member.illustrate(ax, borderless)

        # Draw an inner member #
        if self.innerLayer is not None:
            self.innerLayer.illustrate(ax, borderless)

        # Show the plot if this call originated the drawing #
        if showAfter:
            lim = self.outerRadius * 1.25
            plt.xlim([-lim, lim])
            plt.ylim([-lim, lim])

            plt.show()

    def calculateMass(self, length, breakDownList=None):
        """Calculates the mass of the layer and all layers it encompasses

        Args:
            length (float): Length of the layer in meters. 
            breakDownList (list, optional): A list in which to input mass breakdown information. Defaults to None.

        Returns:
            float: The mass of the layer and all layers it encompasses
        """

        self.assignLayerLengths(length)
        return self._calcLayerMass(breakDownList=breakDownList)

    def _calcLayerMass(self, breakDownList=None):
        """Performs the mass calculation for this layer, calling its corresponding function for all layers below it. 

        Args:
            breakDownList (list, optional): List in which to input mass breakdown information. Defaults to None.

        Raises:
            ValueError: Specified material could not be found in the material database

        Returns:
            float: the mass of this layer and all layers it encompasses
        """

        # Grab the density of the material and the multiplier for the length, throwing an error if not found #
        material_entry = DB.searchEntry(databases.material_db, {'material_name':self.layerMaterial.lower()})
        if material_entry.empty:
            raise ValueError("Unable to find specified material: %s!" % self.layerMaterial)

        # Grab the material's density and grab multiplier to go between mm and m
        density = material_entry["density"].item()
        m_to_mm_mult = DB.build_multiplier("m", "mm")

        # Set length variables in terms of mm and m for readability #
        length_mm = self.length * m_to_mm_mult 
        fill_length = self.length 

        # Length adjustment for helixed fill 
        if(self.helixAngle != 90):
            if(self.helixAngle == 0):
                raise ValueError("The helix angle must have a positive, nonzero value!")

            # Calculate "lead" of the helixed fill # 
            startCircumference = (self.outerRadius - self.innerRadius) * np.pi
            fillLead = startCircumference * np.tan(np.deg2rad(self.helixAngle))

            # Calculate the length of the helixed fill #
            numArcs = length_mm / fillLead 
            helixLen = np.sqrt(fillLead ** 2 + startCircumference ** 2)
            length_mm = helixLen * numArcs

        totalVolume = (self.outerRadius**2 - self.innerRadius**2) * np.pi * length_mm * self.fillRatio

        layerMass = 0

        for member in self.memberList:
            layerMass += member._calcLayerMass(breakDownList=breakDownList)
            memVolume = member.length * np.pi * member.outerRadius ** 2 * m_to_mm_mult
            totalVolume -= memVolume

        layerMass += totalVolume * density

        # If a list was passed, add this layer's mass to it #
        if breakDownList is not None:
            breakDownList.append([self.name, self.layerPath, totalVolume * density, self.layerMaterial, self.length])

        # Add the mass of the layer(s) below in this one #
        if self.innerLayer is not None:
            layerMass += self.innerLayer._calcLayerMass(breakDownList=breakDownList)

        return layerMass
    

    def assignLayerLengths(self, length):
        """Assigns lengths to all layers encapsulated by this one in the tree. 

        Args:
            length (float): Base length to assign

        Raises:
            ValueError: Length cannot be less than 0
            ValueError: Helix angle of members must within this layer must be positive and nonzero
            ValueError: Helix angle is too small to result in a lead that is >= each member's diameter. 
        """

        if length < 0:
            raise ValueError("A length of less than 0 was passed!")

        # Assign the length of this layer and any wrapped layers to what was passed #
        self.length = length 

        if self.innerLayer is not None:
            self.innerLayer.assignLayerLengths(length)

        # If we have members, calculate the length from this length and their helix #
        if len(self.memberList) > 0:

            m_to_mm_mult = DB.build_multiplier("m", "mm")
            mm_to_m_mult = DB.build_multiplier("mm", "m")

            length_mm = self.length * m_to_mm_mult
            length_m = length

             # Check the helix angle #
            if self.helixAngle <= 0 or self.helixAngle > 90:
                raise ValueError("The helix angle must have a positive, nonzero value!")
            
             # No helixing is performed if the angle is 90, member length is what was passed #
            if self.helixAngle == 90:
                actLength = length_m

                for member in self.memberList:
                    member.assignLayerLengths(actLength)

            else:
                
                # Find the largest lead, all helixed members have to have the same lead #
                maxLead = 0
                for member in self.memberList:
                    
                    # Calcualte helix parameters #
                    startCircumference = member.polarCoords[0] * 2 * np.pi
                    memberLead = startCircumference * np.tan(np.deg2rad(self.helixAngle))

                    # Update the max lead #
                    if memberLead > maxLead:
                        maxLead = memberLead

                # Use each member's lead to calculate the actual length #
                for member in self.memberList:
                    
                    if maxLead < member.outerDiameter:
                        raise ValueError("Helix angle of %.2f deg is too small for this member's diameter: \
                                         %f mm, it has a lead of: %f mm! " % 
                                         (self.helixAngle, member.outerDiameter, memberLead))

                    # Calculate the # of distinct arcs in the tether, and adjust the length for each arc #
                    numArcs = length_mm / maxLead 
                    helixLen = np.sqrt(maxLead ** 2 + startCircumference ** 2)
                    actLength = helixLen * numArcs * mm_to_m_mult

                    member.assignLayerLengths(actLength)


                
    def countLayers(self, list=None):
        """Returns the # of layers that this layer encapsulates. 

        Args:
            list (list, optional): List passed which is used to keep track of layers. Defaults to None.

        Returns:
            int: The number of layers below this one.
        """

        layerList = None
        if list is None:
            layerList = []
        else:
            layerList = list

        layerList.append(self)

        for mem in self.memberList:
            mem.countLayers(list=layerList)
        if self.innerLayer is not None:
            self.innerLayer.countLayers(list=layerList)

        return len(layerList)
    
    def findWires(self, wireList, wireType):
        """Finds and builds a list of the layer paths for all wires mathcing the passed wire type. 

        Args:
            wireList (list): List in which to input wire paths. 
            wireType (str): String representing the wire type. Must be "electrical" or "optical"

        Raises:
            ValueError: An invalid wire type is passed
        """

        if wireType != "electrical" and wireType != "optical":
            raise ValueError("Invalid wire type passed: %s!" % wireType)

        if issubclass(type(self), Wire) and self.wireType == wireType:
            wireList.append(self.layerPath)
        else:
            for mem in self.memberList:
                if issubclass(type(mem), Wire) and mem.wireType == wireType:
                    wireList.append(mem.layerPath)    
                else:
                    mem.findWires(wireList, wireType)

            if self.innerLayer is not None:
                self.innerLayer.findWires(wireList, wireType)


    def _assignLayerPath(self, layerNum, prefix=""):
        """Assigns a path to the layer in the context of the overall tree. Should not
        be called by the user. 

        Args:
            layerNum (int): The number of layers deep we are.
            prefix (str, optional): The path to prefix this one. Defaults to "".
        """

        self.layerPath = prefix + "L%d" % layerNum
        
        for idx, mem in enumerate(self.memberList):
            mem._assignLayerPath(1, prefix=self.layerPath + "M%d" % (idx + 1))
        
        if self.innerLayer is not None:
            self.innerLayer._assignLayerPath(layerNum + 1, prefix=prefix)

    def printLayerPaths(self):
        """Prints the names of all layers from this point downwards in the tree.
        """
        print(self.layerPath)

        for mem in self.memberList:
            mem.printLayerPaths()

        if self.innerLayer is not None:
            self.innerLayer.printLayerPaths()


    def __eq__(self, other) : 
        """Allows comparison of object to others using the == operator. Checks for equality not identicality. 

        Args:
            other (Layer): Object to compare against

        Returns:
            bool: Whether ALL object attributes equal the other
        """
        return self.__dict__ == other.__dict__


class Wire(Layer):

    def __init__(self, wireType, insulation_material, name, insulation_thickness, innerLayer=None, color=None):
        """Wire constructor. This should be the outermost layer of the wire.

        Args:
            wireType (str): Type of wire this object represents. Must be "electrical" or "optical".
            insulation_material (str): Name of the insulation/jacket material.
            name (str): Name of the layer this wire object represents (insulation usually).
            insulation_thickness (float): Thickness of the insulation/jacket in mm
            innerLayer (Layer, optional): Rest of the wire (defined by layers) that this encapsulates. Defaults to None.
            color (str, optional): Color of the outermost layer. All named CSS colors are valid. Defaults to None.

        Raises:
            ValueError: Thrown if an invalid wireType is passed.
        """

        if wireType != "electrical" and wireType != "optical":
            raise ValueError("Invalid wire type specified: %s!" % wireType)

        self.wireType = wireType

        super().__init__(insulation_material, name, insulation_thickness, innerLayer, None, True, color)

    @classmethod
    def from_entry(cls, wireType, entry):
        """Class method build the Wire object from a wire or optical fiber entry. 

        Args:
            wireType (str): Type of wire this object represents. Must be "electrical" or "optical".
            entry (pd.Dataframe): Dataframe of the wire entry, must follow the format for wires/optical fibers. 

        Raises:
            ValueError: Thrown if an invalid wireType is passed.
        """
        if wireType == "electrical":

            return(cls._set_from_electrical_entry(wireType, entry))

        elif wireType == "optical":
            
            return(cls._set_from_optical_entry(wireType, entry))

        else:
            raise ValueError("Attempting to set from entry with an invalid wire type specified: %s!" % wireType)
        
    @classmethod
    def _set_from_electrical_entry(cls, wireType, entry):
        """Return an electrical Wire object defined from an entry.

        Args:
            wireType (str): Type of wire this object represents. Must be "electrical" or "optical".
            entry (pd.Dataframe): Dataframe of the wire entry, must follow the format for wires/optical fibers.

        Raises:
            ValueError: Insulation OD and Insulation thickness not specified by wire entry.
            ValueError: Conductor OD was not specified by the wire entry.
            ValueError: Insulation material was not specified by the wire entry.
            ValueError: Conductor material was not specified by the wire entry.

        Returns:
            Wire: Wire object built using the constructor. 
        """

        # Break out values #
        conductor_od = entry["conductor_od"].item()
        insulator_od = entry["insulation_od"].item()
        insulator_thick = entry["insulation_thickness"].item()
        conductor_mat = entry["conductor_material"].item()
        insulator_mat = entry["insulation_material"].item()

        # Throw errors if any values are missing #
        if insulator_od is None and insulator_thick is None:
            raise ValueError("Neither insulation OD nor insulation thickness specified by wire entry!")
        if conductor_od is None:
            raise ValueError("No conductor OD was specified by the wire entry!")
        if insulator_mat is None:
            raise ValueError("No insulation material was specified by the wire entry!")
        if conductor_mat is None:
            raise ValueError("No conductor material was specified by the wire entry!")

        # Calculate insulation OD from thickness if not present #
        if insulator_od is None:
            insulator_od = conductor_od + (2 * insulator_thick)
        elif insulator_thick is None:
            insulator_thick = (insulator_od - conductor_od) / 2

        # Determine whether the conductor material exists in the database #
        conductor_mat_entry = None
        try:
            conductor_mat_entry = databases.get_material_entry(conductor_mat)
        except ValueError:
            pass

        # Determine whether the insulator material exists in the database
        insulator_mat_entry = None
        try:
            insulator_mat_entry = databases.get_material_entry(insulator_mat)
        except ValueError:
            pass

        # Set the conductor color, using a default if not found #
        if conductor_mat_entry is not None:
            conductor_color = DB.get_material_property(conductor_mat_entry, "color")
        else:
            conductor_color = "goldenrod"

        # Set the insulator color, using a default if not found #
        if insulator_mat_entry is not None:
            insulator_color = DB.get_material_property(insulator_mat_entry, "color")
        else:
            insulator_color = "dimgrey"

        #TODO: check for overrides, & generate override string if found
        #Search material DB for override string, if not found create a new material entry and add it

        conductor = Layer(conductor_mat, "%s Conductor" % conductor_mat, conductor_od / 2, color=conductor_color)

        return cls(wireType, insulator_mat, "%s Insulated Wire" % insulator_mat, insulator_thick, innerLayer=conductor, 
                   color=insulator_color)
    

    @classmethod
    def _set_from_optical_entry(cls, wireType, entry):
        """Return an optical Wire object defined from an entry

        Args:
            wireType (str): Type of wire this object represents. Must be "electrical" or "optical".
            entry (pd.Dataframe): Dataframe of the wire entry, must follow the format for wires/optical fibers.

        Raises:
            ValueError: Fiber type is not specified by the wire entry.
            ValueError: Core OD is not specified by a multimode entry.
            ValueError: Cladding OD is not specified by a multimode entry.
            ValueError: Incorrect wire type specified.
            ValueError: Jacket OD is not specified by optical fiber entry. 
            ValueError: Jacket material is not specified by optical fiber entry.

        Returns:
            Wire: Wire object built by using the constructor. 
        """
        
        # Set some industry standard constants (from the FOA) #
        singlemode_core = 0.009
        singlemode_clad = 0.125

        # Break out values #
        fiber_type = entry["type"].item()
        core_od = entry["core_od"].item()
        cladding_od = entry["cladding_od"].item()
        jacket_od = entry["jacket_od"].item()
        jacket_material = entry["jacket_material"].item()

        # Throw errors if any necessary values are missing, set a few constants if necessary #
        if fiber_type is None:
            raise ValueError("No fiber type is specified by the entry!")
        elif fiber_type == "multimode":
            if core_od is None:
                raise ValueError("No core OD specified for multimode type fiber!")
            if cladding_od is None:
                raise ValueError("No cladding OD specified for multimode type fiber!")
        elif fiber_type == "singlemode":
            core_od = singlemode_core
            cladding_od = singlemode_clad
        else:
            raise ValueError("Incorrect type specified: %s!" % fiber_type)
    
        if jacket_od is None:
            raise ValueError("No jacket OD specified for fiber entry!")
        if jacket_material is None:
            raise ValueError("No jacket material is specified for fiber entry!")

        # Determine whether the conductor material exists in the database #
        jacket_mat_entry = None
        try:
            jacket_mat_entry = databases.get_material_entry(jacket_material)
        except ValueError:
            pass

        # Set jacket color to a default if its entry is not found in the material database #
        if jacket_mat_entry is not None:
            jacket_color = DB.get_material_property(jacket_mat_entry, "color")
        else:
            jacket_color = "cornflowerblue"

        #TODO: check for overrides, & generate override string if found
        #Search material DB for override string, if not found create a new material entry and add it

        # Calculate thicknesses #
        core_thickness = core_od / 2
        clad_thickness = (cladding_od - core_od) / 2
        jacket_thickness = (jacket_od - cladding_od) / 2

        # Set up core/cladding layers #
        core = Layer('silica', "silica Core", core_thickness)
        cladding = Layer('silica', "silica Cladding", clad_thickness, innerLayer=core)

        return cls(wireType, jacket_material, "%s Jacket Fiber Optic" % jacket_material, jacket_thickness, 
                   innerLayer=cladding, color=jacket_color)
    
    def __eq__(self, other) : 
        """Allows comparison of object to others using the == operator. Checks for equality not identicality. 

        Args:
            other (Layer): Object to compare against

        Returns:
            bool: Whether ALL object attributes equal the other
        """
        return self.__dict__ == other.__dict__


# A class representing the overall tether design. User interacts most with this class #
class RoundTetherDesign():
    """Object meant to describe a round tether design. Contains overall information about the tether 
    and contains the tree of layers used to describe the tethers specific geoemtry. Used as input into
    analyses.
    """

    def __init__(self, name, layer, length) -> None:
        """Constructor for RoundTetherDesign class. This class encapsulates a layer tree with other
        information that is used to help inform Analyses. 

        Args:
            name (str): Name of the RoundTetherDesign object.
            layer (Layer): Outermost layer of the tree. 
            length (float): Length of the tether in meters.

        Raises:
            ValueError: Thrown if the tether length is negative
        """

        self.layer = layer

        # Object name #
        self.name = name

        # Lengthwise Properties #
        if length < 0:
            raise ValueError("Cannot have a negative length value!")
        self.length = length          # (m)

        # Set names for the wires #
        self.assignLayerPaths()

        # Geometry and physical properties #
        self.diameter = self.layer.outerDiameter       # (mm)
        self.radius = self.layer.outerRadius           # (mm)
        
        # Calculate lengths and estimate mass, MBR, and strength
        self._calculateMechanicalProperties()

    def _calculateMechanicalProperties(self):
        """ Calculates mechanical properties of tether on construction or on updates. 
        """
        # Geometry and physical properties #
        self.layer.assignLayerLengths(self.length)
        self.mass = self.calculateMass()
        self.minBendRadius = self.calculateMinBendRadius()

        # Find strength values #
        self.strengthUpperBound = self.unhelixStrength(output=False)
        self.strengthLowerBound = self.straightRunStrength(output=False)



    def tetherDetails(self, recursive=False, output=True, breakdownList=None):
        """Prints the details of the tether.

        Args:
            recursive (bool, optional): Prints the details of each layer recursively if desired. Defaults to False.
        """
        if(output):
            print("--- %s Round Tether Details ---" % self.name)
            print("  Tether Length: %f m" % self.length)
            print("  Tether Diameter: %f mm" % self.layer.outerDiameter)
            print("  Tether Radius: %f mm" % self.layer.outerRadius)
            print("  Tether Mass: %f g" % self.mass)
            print("  Strength Bounds: {} N, {} N".format(self.strengthLowerBound, self.strengthUpperBound))
            print("  Min Bend Radius: {} mm".format(self.minBendRadius))

        if recursive:
            self.layer.layerDetails(recursive=True, tabNum=1, output=output, breakdownList=breakdownList)

    def illustrate(self, borderless=False, figdim=5):
        """Illustrates a given tether design.

        Args:
            borderless (bool, optional): Whether to omit layer borders. Defaults to False.
            figdim (int, optional): Dimension of each side of the figure (inches). Defaults to 5.
        """

        self.layer.illustrate(borderless=borderless, figdim=figdim)


    def calculateMinBendRadius(self):
        """Calculates the minimum bend radius of the tether in mm.
        Based on rules of thumb from:
         https://www.thefoa.org/tech/ref/install/bend_radius.html#:~:text=The%20normal%20recommendation%20for%20fiber,10%20times%20the%20cable%20diameter
         https://www.anixter.com/en_us/resources/literature/wire-wisdom/minimum-bend-radius.html

        """

        tetherRule = 8 * self.diameter

        elecPaths = self.findWires("electrical")
        fiberPaths = self.findWires("optical")

        maxWireOD = 0
        maxFiberOD = 0

        for path in elecPaths:
            wire = self.getLayerAtPath(path)
            if wire.outerDiameter > maxWireOD:
                maxWireOD = wire.outerDiameter

        for path in fiberPaths:
            fiber = self.getLayerAtPath(path)
            if fiber.outerDiameter > maxFiberOD:
                maxFiberOD = fiber.outerDiameter

        wireRule = 12 * maxWireOD
        fiberRule = 20 * maxFiberOD # TODO add Curtis's math (ask for source)

        return max(tetherRule, wireRule, fiberRule)


    def calculateMass(self, breakdown=False, verbose=False):
        """Calculates the mass of the tether in grams. 

        Args:
            breakdown (bool, optional): Whether to return a breakdown list. Defaults to False.

        Returns:
            float: the mass in grams
            list: A list of layer data, each index consisting of a layer's [name, mass, length]
        """
        
        breakDownList = []
        mass = self.layer.calculateMass(self.length, breakDownList=breakDownList)

        if(verbose):
            print("--- %s Mass Breakdown ---" % self.name)
            print(" - Total Mass: %fg" % self.mass)
            for layerBreakdown in breakDownList:
                print("   - Layer Name: %s, Layer Path: %s" % (layerBreakdown[0], layerBreakdown[1]))
                print("     - Mass: %fg" % layerBreakdown[2])
                print("     - Material: %s" % layerBreakdown[3])
                print("     - Length: %fm" % layerBreakdown[4])

        if breakdown:
            return mass, breakDownList
        else:
            return mass
        

    def straightRunStrength(self, output=True, verbose=False):
        """Calculates the strength of the tether in newtons, assuming
        that all materials see load immediately. 

        Returns:
            float: The strength in newtons. 
        """
        
        tetherLayers = []
        self.layer.countLayers(list=tetherLayers)

        min_strain = np.inf
        eff_spring_const = 0

        layerDict = {}

        for layer in tetherLayers:
            
            # Grab yield stress and young's modulus #
            material_entry = databases.get_material_entry(layer.layerMaterial)
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
            print("  Tether Elongation: {} m".format(min_strain * self.length))

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


    def unhelixStrength(self, output=True, verbose=False):
        
        """Calculates the strength of the tether, allowing helixed members to
        straighten. (Not all members see load immeidately)

        Args:
            output (bool, optional): Whether to print analysis results. Defaults to True.
            verbose (bool, optional): Whether to include per-layer results. Defaults to False.

        Returns:
            _type_: _description_
        """

        tetherLayers = []
        self.layer.countLayers(list=tetherLayers)

        # Sort layers by ascending order of length #    
        tetherLayers.sort(key=lambda x: x.length)

        # the length the entire tether needs to stretch to in order to hit the yield strain of this member #
        yield_length_tether = [] 
        layerDict = {}

        for layer in tetherLayers:

            layerDict[layer.layerPath] = [layer.name, layer.layerPath]

            # Grab yield stress and young's modulus #
            material_entry = databases.get_material_entry(layer.layerMaterial)
            yield_stress = DB.get_material_property(material_entry, "stress")
            young_mod = DB.get_material_property(material_entry, "elastic_modulus")

            if young_mod == 0:
                continue
            
            # Calculate effective tether length, and update min strain #
            strain = yield_stress / young_mod
            yield_length = layer.length * (1 + strain)
            yield_length_tether.append(yield_length)

        # Find the breaking length of the tether (first layer to hit yield strain w.r.t. tether length) # 
        tether_break_len = min(yield_length_tether)

        break_force = 0

        forceCarryingLayers = []


        # Build a list of layers that will see stretch & contribute to strength #
        for layer in tetherLayers:
            if layer.length < tether_break_len:

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
                spring_const = (young_mod * csa * gpa_to_pa) / layer.length 

                layerDict[layer.layerPath].extend([spring_const * (tether_break_len - layer.length), (tether_break_len - layer.length), spring_const])

                break_force += spring_const * (tether_break_len - layer.length)

            else:
                layerDict[layer.layerPath].extend([0, 0, 0])

        if(output):
            print("--- Tether Unhelix-Allowed Strength Analysis ---")
            print("  Calculated Strength: {} N".format(break_force))
            print("  Tether Elongation: {} m".format(tether_break_len - self.length))
            
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


    def thermalReport(self):
        """Prints out properties useful to performing a thermal analysis on the tether design. 
        
        Properties printed include the name, path, material, thermal conductivity, dimensions, and position of each layer in the tether.
        
        Args: None
        
        Returns: None
        """
        
        breakdownList = []
        self.tetherDetails(recursive=True, output=False, breakdownList=breakdownList)
        
        # [self.name, self.layerPath, self.layerMaterial, self.innerRadius, self.outerRadius, self.layerThickness, (self.x, self.y)]
        
        for entry in breakdownList:
            material_entry = databases.get_material_entry(entry[2])
            thermal_cond = DB.get_material_property(material_entry, "thermal_cond")
            min_temp = DB.get_material_property(material_entry, "min_temp")
            max_temp = DB.get_material_property(material_entry, "max_temp")
            print("Layer: %s, Path: %s" % (entry[0], entry[1]))
            print("  - Material: %s, Thermal Conductivity: %f W/(m*K)" % (entry[2], thermal_cond))
            print("  - Min Temp: %f C, Max Temp: %f C" % (min_temp, max_temp))
            print("  - Dimensions:")
            print("     - Inner Radius: %f mm" % entry[3]) 
            print("     - Outer Radius Radius: %f mm"% entry[4]) 
            print("     - Thickness: %f mm"% entry[5]) 
            print("  - Centered at: (%f, %f)" % (entry[6][0], entry[6][1]))
            print("\n")



    def findWires(self, wireType):
        """Returns the path of all wires found in the tether

        Args:
            wireType (string): The type of wire to search for, must be electrical or optical

        Raises:
            ValueError: Invalid wire type passed

        Returns:
            list: A list of paths of all found wires
        """
        if wireType != "electrical" and wireType != "optical": 
            raise ValueError("Attempting to search for an invalid wire type: %s" % wireType)

        wireList = []

        self.layer.findWires(wireList, wireType)

        return wireList
    
    def buildLayerList(self, layerPaths):
        """Builds a list of layer objects in the tether from a set of paths

        Args:
            layerPaths (list): A list of paths to traverse

        Returns:
            list: A list of found layer objects
        """

        layerList = []

        for path in layerPaths:
            layer = self.getLayerAtPath(path)
            if layer is not None:
                layerList.append(layer)

        return layerList

    def assignLayerPaths(self):
        """Assigns layer paths for the tree within this tether. 
        """

        self.layer._assignLayerPath(1)

    def getLayerAtPath(self, path):
        """Traverses the tree and returns the layer at the end

        Args:
            path (string): The path within the tree to follow 
            (L2 to go 2 layers in, M2 to go to the second member, 
            L2M2L1 to go to the second layer-->first member-->first layer)

        Raises:
            ValueError: The passed path and the path of the found member differ

        Returns:
            Layer: The layer object found at the end of the path
        """
        regStr = r"[A-Z]\d+"
        opList = re.findall(regStr, path)

        object = self.layer

        for op in opList:
            moveDir = re.findall(r"[A-Z]", op)[0]
            num = int(re.findall(r"\d+", op)[0]) - 1

            if moveDir == "L":
                for i in range(0, num):
                    object = object.innerLayer
            elif moveDir == "M":
                object = object.memberList[num]

        if object.layerPath != path:
            raise ValueError("Something went wrong with pathing! Expected: \
                             %s but got: %s" % (path, object.layerPath))

        return object
    
    def printLayerPaths(self):
        """Print all of the layer paths within the tree
        """
        self.layer.printLayerPaths()
        
    def __eq__(self, other) : 
        """Allows comparison of object to others using the == operator. Checks for equality not identicality. 

        Args:
            other (RoundTetherDesign): Object to compare against

        Returns:
            bool: Whether ALL object attributes equal the other
        """
        return self.__dict__ == other.__dict__
    
    def adjustLayerPosition(self, path, x_delta, y_delta):
        """ Manually Adjust layer coordinates by x_delta and y_delta

        Args:
            path (string): The path of the target layer within the tree.
            x_delta (float): The adjustment to the x coordinate of the layer's position.
            y_delta (float): The adjustment to the y coordinate of the layer's position.
        """
        
        layer = self.getLayerAtPath(path)
        layer.x += x_delta
        layer.y += y_delta
        layer.polarCoords = cart_to_polar([layer.x, layer.y])
        layer.update_member_positions()

    def adjustLayerThickness(self, path, thickness):
        """ Manually adjust the thickness of a layer. Does not check for validity. 
        
        Args: 
            path (string): The path of the target layer within the tree.
            thickness (float): The desired thickness of the layer 
        """

        layer = self.getLayerAtPath(path)
        layer.thickness = thickness
        layer.outerDiameter = layer.innerDiameter + layer.thickness
        layer.outerRadius = layer.outerDiameter / 2

    def updateLayerSizeFromInner(self, path):
        """ Manually adjust the thicknes sof a path from its inner layer, using its set thickness

        Args:
            path (string): The path of the target layer within the tree
        """

        layer = self.getLayerAtPath(path)
        layer._setLayerDimensions()

    def addMemberToLayer(self, path, member, x, y):
        """ Manually adds a member to a given layer at the chosen coordinates
        
        Args:
            path (string): The layer which a member is added to. 
            member (Layer): The member to add.
            x (float): The x coordinate of the member. 
            y (float): The y coordinate of the member. 
        """

        # Grab layer and make a copy of the desired member #
        layer = self.getLayerAtPath(path)
        newMem = deepcopy(member)

        # Assign coordinates for the new member and update encapsulated layers # 
        newMem.x = x
        newMem.y = y
        newMem.polarCoords = cart_to_polar([newMem.x, newMem.y])
        
        newMem.update_member_positions()
        
        # Add it to the list and then make sure all paths are up-to-date # 
        layer.memberList.append(newMem)
        self.assignLayerPaths()
        self._calculateMechanicalProperties()



# Some helper functions #
def midpoint(p1, p2):
    """Finds the midpoint (cartesian) between points 1 and 2

    Args:
        p1 (list(float, float)): First point
        p2 (list(float, float)): Second point
    """

    # Calculate midpoint #
    x_val = (p1[0] + p2[0]) / 2
    y_val = (p1[1] + p2[1]) / 2
    return [x_val, y_val]


def checkOverlap(container, inner):
    """Checks whether there is any overlap between an inner layer and its container (meaning collision)

    Args:
        container (Layer): Containing (outer) layer
        inner (Layer): Inner layer that may overlap
    Returns:
        Bool: Whether overlap occurs
    """

    overlap = False
    innerPolarR = inner.polarCoords[0]
    outerPolarR = container.polarCoords[0]
    # Round the comparison to avoid errors in floating point approx. 
    if np.round(innerPolarR + inner.outerRadius, 9) > np.round(outerPolarR + container.outerRadius, 9):
        overlap = True


    return overlap


def check_tangent(p1, p2, r1, r2):
    """Checks whether two circles are at least tangent to each other, rounded to the 12th decimal place 
    (to avoid errors with floating point approximation)

    Args:
        p1 (list(float, float)): Center point of first circle
        p2 (list(float, float)): Cetner point of second circle
        r1 (int): Radius of first circle
        r2 (int): Radius of second circle

    Returns:
        Bool: Whether the circle are tangent
    """

    tangent = False
    pointDist = round(np.linalg.norm(np.array(p1) - np.array(p2)),12)
    rDist = round(r1 + r2, 12)

    if pointDist >= rDist:
        tangent = True

    return tangent


def cart_to_polar(point):
    """Converts a set of cartesian coordinates to polar

    Args:
        point (tuple(float, float)): Point in cartesian coordinates to convert

    Returns:
        tuple(float, float): Polar coordinates
    """

    rad = np.sqrt((point[0] ** 2) + (point[1] ** 2))
    theta = np.real(np.arctan2(point[1], point[0]))
    return [rad, theta]


def polar_to_cart(point):
    """Converts a set of polar coodinates to cartesian

    Args:
        point (tuple(float, float)): Point in polar coordinates to convert

    Returns:
        tuple(float, float): Cartesian coordinates
    """
    point[1] = np.real(point[1])
    x = point[0] * np.cos(point[1])
    y = point[0] * np.sin(point[1])
    return [x, y]

def get_angle_C(sideA, sideB, sideC):
    """Calculates angle C of a triangle using the lengths of each side. 

    Args:
        sideA (float): length of side A
        sideB (float): length of side B
        sideC (float): lenght of side C

    Returns:
        float: Angle of side C in radians
    """
    # Calculate angle C (theta) # 
    upperTerm = (sideC ** 2) - (sideA ** 2) - (sideB ** 2)
    lowerTerm = -2 * sideA * sideB
    return np.emath.arccos(upperTerm / lowerTerm)
