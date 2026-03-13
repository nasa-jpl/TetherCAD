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

import databases.DatabaseControl as DB

databases = DB.DatabaseControl()

# Get all boundary points between two layers
def getBoundaries(l_start, l_end, layer):

  flip = False  

  if l_start.contains(l_end) or l_end.contains(l_start):
    """ nested copper layers """
    m = 0
    b = l_start.y
    
    # if l_start.contains(l_end):
    #   order = 1
    # else:
    #   order = -1
    
  else:
    """ Two discrete conductors """
    # put line in mx+b form
    dx = l_end.x - l_start.x
    if abs(dx) > 1e-3:
      m = (l_end.y - l_start.y) / dx
      b = l_start.y-m*l_start.x
    else:
    # Check for near-vertical line, flip axes
      dy = l_end.y - l_start.y
      m = (l_end.x - l_start.x) / dy
      b = l_start.x-m*l_start.y
      flip = True
    
    # # Set search order
    # if l_start.x < l_end.x:
    #   order = 1
    # else:
    #   order = -1
  
  # Compute all circle-line intersections
  layers = layer.enumerateLayers()
  
  if flip:
    ipts = []
    for l in layers:
     
      y = getIntersectionsY(m, b, l.x, l.y, l.outerRadius)
      if y[0] == y[1] or np.imag(y[0]) > 0:
        continue
      else:
        ipts.append({'layer': l, 'x':y[0]*m+b, 'y':y[0]})
        ipts.append({'layer': l, 'x':y[1]*m+b, 'y':y[1]})
        
    if l_start.y > l_end.y:
      sign = -1
    else:
      sign = 1
    def kfunc(val):
      return sign*val['y']
    
  else:
    ipts = []
    for l in layers:
     
      x = getIntersectionsX(m, b, l.x, l.y, l.outerRadius)
      if x[0] == x[1] or np.imag(x[0]) > 0:
        continue
      else:
        ipts.append({'layer': l, 'x':x[0], 'y':x[0]*m+b})
        ipts.append({'layer': l, 'x':x[1], 'y':x[1]*m+b})
    
    if l_start.x > l_end.x:
      sign = -1
    else:
      sign = 1    
    
    def kfunc(val):
      return sign*val['x']
    
  list.sort(ipts, key=kfunc)

  

  
  # Create series of edges through which we're calculating the dielectric strength
  
  cond = [l_start, l_end]
  entered = None
  # inside_start = False
  start_edge = None
  edges = []
  for pt in ipts:
    l = pt['layer']
    # print(pt)
    
    if entered is None and l in cond:
      #entering first conductor
      entered = l
      edges = []
      continue
    
    if entered is not None:
      if l == entered:
        # Now leaving outer surface of first conductor, start edge list
        edges = [pt]
      elif l in cond:
        # found the second conductor, exit
        edges.append(pt)
        break
      else:
        # found another insulator, add it to the edge list
        edges.append(pt)
      
    # # print(pt)
    # if pt['layer'] == l_start:
    #   if not inside_start:
    #     inside_start=True
    #   else:
    #     start_edge = pt
    # if pt['layer'] == l_end:
    #   edges.append(pt)
    #   break
    # if start_edge is not None:
    #   edges.append(pt)
  # print(ipts)
  # print("EDGES")
  # print(edges)
  return edges

# Returns the lengths and materials in the straight-line path between two conductors
def getDielectricStrength(l_start, l_end, bounding_layer):
  
  # Always start with the outer layer
  if l_end.contains(l_start):
    l = l_end
    l_end = l_start
    l_start = l
  
  edges = getBoundaries(l_start, l_end, bounding_layer)

  current_layer = l_start
  prev_edge = edges[0]
  materials = []
  for edge in edges:
    if edge['layer'] == current_layer: # We're leaving a layer
      p0 = np.array([prev_edge['x'], prev_edge['y']])
      p1 = np.array([edge['x'], edge['y']])
      dist = np.linalg.norm(p0-p1)
      materials.append({'length': dist, 'layer' : current_layer})
      current_layer = current_layer.parent
    else: # We're entering a layer
      p0 = np.array([prev_edge['x'], prev_edge['y']])
      p1 = np.array([edge['x'], edge['y']])
      dist = np.linalg.norm(p0-p1)
      materials.append({'length': dist, 'layer' : current_layer})  
      current_layer = edge['layer']
    prev_edge = edge
  
  strength = 0
  air = databases.get_material_entry('air')
  for mat in materials:
    mat['material_name'] = mat['layer'].material_entry['material_name']
    
    breakdown_strength = mat['layer'].material_entry['dielectric_strength'].values[0]
    if mat['layer'].fillRatio < 1:
      breakdown_strength = air['dielectric_strength'].values[0]
    
    mat['breakdown_voltage'] = breakdown_strength*mat['length']
    
    # print("Length: %f, Strength: %f" % (mat['length'], breakdown_strength))
    strength += mat['breakdown_voltage']
    
  
  return strength, materials


def getIntersectionsX(m, b, x0, y0, r):
  #y = mx + b
  #(x-x0)^2 + (y-y0)^2 = r^2
  #(x-x0)^2 + (m x + b -y0)^2 = r^2
  #x^2 - 2*x0*x + x0^2 + m^2 x^2 + 2*m*(b-y0)*x + (b-y0)^2 = r^2
  #(1+m^2) x^2 + (2*m*(b-y0) - 2*x0 )*x + (b-y0)^2 - r^2 + x0^2 = 0

  p = [(1+m**2), (2*m*(b-y0) - 2*x0), (b-y0)**2 - r**2 + x0**2]
  r = np.roots(p)
  
  # print(np.array(r)*m+b)
  
  return r

def getIntersectionsY(m, b, x0, y0, r):
  #x = my + b
  #(x-x0)^2 + (y-y0)^2 = r^2
  #(my + b-x0)^2 + (y-y0)^2 = r^2
  #y^2 - 2*y0*y + y0^2 + m^2 y^2 + 2*(b-x0)*y + (b-x0)^2 = r^2
  #(1+m^2) y^2 + (2*m*(b-x0) - 2*y0 )*y + (b-x0)^2 - r^2 + y0^2 = 0

  p = [(1+m**2), (2*m*(b-x0) - 2*y0), (b-x0)**2 - r**2 + y0**2]
  r = np.roots(p)
  
  return r


  
def drawDielectricLine(l_start, l_end, bounding_layer):
  import matplotlib.pyplot as plt
  
  bounding_layer.illustrate()
  
  ipts = getBoundaries(l_start, l_end, bounding_layer)
  plt.plot([ipts[0]['x'], ipts[-1]['x']], [ipts[0]['y'], ipts[-1]['y']])
  for pt in ipts:
    plt.plot(pt['x'], pt['y'], marker='o')
    
    
    
    
    
    
    
    
    

