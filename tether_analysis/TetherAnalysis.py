""""
Copyright 2024, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government 
Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the 
California Institute of Technology. This software may be subject to U.S. export control laws. 
By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. 
User has the responsibility to obtain export licenses, or other export authority as may be required before exporting 
such information to foreign countries or providing access to foreign persons."

If you have questions regarding this, please contact the JPL Software Release Authority at x4-2458.
"""


import tether_analysis.TetherDesign as TD
import calculation_libraries.PowerAnalysis as PA
import databases.DatabaseControl as DB
import calculation_libraries.SpoolAnalysis as SP
import calculation_libraries.CommAnalysis as CA
import calculation_libraries.MechAnalysis as MA
import calculation_libraries.DielectricAnalysis as DA
from importlib import reload
import matplotlib.pyplot as plt
import math
import numpy as np

databases = DB.DatabaseControl()

def analysis_1(tether, annotate=True, annotation_spacing=None, temp=200, tether_volt=100, tether_power=10):
  
  pathList = tether.findWires("electrical")

  if len(pathList) % 2 == 0:
    send_wires = []
    return_wires = []
    
    for i in range(0,len(pathList)):
      if i % 2 == 0:
        send_wires.append(pathList[i])
      else:
        return_wires.append(pathList[i])
  else:
    send_wires = None
    return_wires = None

  
  print("Send wires: %s" % str(send_wires))
  print("Return wires: %s" % str(return_wires))
  
  

  if send_wires is not None:
    sols = PA.dc_power_transmission_analysis(tether, 
                                      tether_voltage=tether_volt, 
                                      desired_power=tether_power, 
                                      send_paths=send_wires, 
                                      return_paths=return_wires, 
                                      temp=temp)

  import calculation_libraries.MechAnalysis as MA
  strengthUpperBound = MA.unhelixStrength(tether, output=False)




  # fig, ax = plt.subplots(figsize=(5,5))
  fix, ax = plt.subplots(num=None, figsize=(7, 5), dpi=150)
  ax.set_anchor('W')
  xf = tether.diameter/2+.2
  tether.illustrate(ax=ax)

  ax.spines[['right', 'top']].set_visible(False)
  # arrowprops = dict(color='black', width=.5, headwidth=6, headlength=6)
  arrowprops = dict(color='black', arrowstyle='-')
  an = []

  specs = []
  specs.append("Diameter: %.1f mm" % tether.diameter)
  specs.append("Length: %.0f m" % tether.length)


  if send_wires is not None:
    tether_R = PA.dc_resistance(tether, send_wires, return_wires, temp=temp)
    specs.append("Design Temperature: %.0f C" % temp)
    specs.append("DC Resistance: %.0f Ohms" % tether_R)
    if len(sols) > 0:
      if tether_volt < 5000:
        specs.append("Voltage: %.0f V" % tether_volt)
      else:
        specs.append("Voltage: %.0f kV" % (tether_volt*.001))
      
      if tether_power < 10e3:
        specs.append("Delivered: %.0f W" % tether_power)
      else:
        specs.append("Delivered: %.0f kW" % (tether_power*.001))
      specs.append("Loss: %.0f W" % sols[0]['total_loss'])
      specs.append("Unit Loss: %.4f W/m" % (sols[0]['total_loss']/tether.length))
      specs.append("Efficiency: %.0f %%" % sols[0]['efficiency'])
      
      send = tether.getLayerAtPath(send_wires[0])
      ret = tether.getLayerAtPath(return_wires[0])
      strength, materials = DA.getDielectricStrength(send.innerLayer, ret.innerLayer, tether.layer)
    
      specs.append('Breakdown Voltage: %.02f kV' % (strength/1000))
    else:
      print("Warning, no power solutions found!!")
      
  specs.append("Total Mass: %.2f kg" % (tether.mass/1000))
  specs.append("Unit Mass: %.01f g/m" % (tether.mass/tether.length))
  # specs.append("Breaking Strength: %.0f N" % strengthUpperBound)

  if annotate:
    if annotation_spacing is None:
      annotation_spacing = tether.diameter / 15
    yf = xf
    for s in specs:
      ax.text(xf, yf, s)
      print(s)
      yf -= annotation_spacing


  ax.set_aspect('equal', adjustable='box')
  # ax.set_xlim(-2, 2) 
  # ax.set_ylim(-2, 2)
  plt.rcParams['pdf.fonttype'] = 42
  plt.rcParams['ps.fonttype'] = 42

  # tether.thermalReport()

  tether.tetherDetails(recursive=True)

  # spool_width= 250
  # SP.determine_spool_od(tether, 40, spool_width)


  if send_wires is not None:
    send = tether.getLayerAtPath(send_wires[0])
    ret = tether.getLayerAtPath(return_wires[0])
    ipts = DA.getBoundaries(send.innerLayer, ret.innerLayer, tether.layer)
    strength, materials = DA.getDielectricStrength(send.innerLayer, ret.innerLayer, tether.layer)
  
    print('Dielectric Breakdown Voltage: %f kV' % (strength/1000))
  
    imp = CA.characteristic_impedance(tether, pathList[0], pathList[1], temp=200, frequency=20e6)
    print('Characteristic Impedance (Not validated): %s' % str(imp))

  if len(tether.getLayerByName("FO Cladding")) > 0:
    fo = tether.getLayerByName("FO Cladding")[0].outerDiameter
    fo_r = tether.getLayerByName("FO")[0].polarCoords[0]
    lead = tether.getLayerByName("Fill")[0].memberHelixLead
    R = MA.calcHelixCurvature(lead, fo_r)
    print("Optical fiber lead: %.02f mm / rev" % lead)
    print("Optical fiber lead/diameter: %.02f" % (lead/tether.diameter))
    print("Optical fiber diamaeter: %.3f" % fo)
    print("Optical fiber bend radius: %.02f mm" % R)
    print("Optical curvature is %f X diameter" % (R/fo))

  st = tether.getLayerByName("Strength")
  if len(st) > 0:
    st = tether.getLayerByName("Strength")[0]
    d = (st.outerDiameter + st.innerDiameter)/2
    cir = d * math.pi
    lead = cir * np.tan(np.deg2rad(st.helixAngle))
    print("Strength layer lead: %.02f mm / rev" % lead)
    
  return ax
