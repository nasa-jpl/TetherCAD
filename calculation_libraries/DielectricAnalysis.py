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
    m = (l_end.y - l_start.y) / (l_end.x - l_start.x)
    b = l_start.y-m*l_start.x
    
    # # Set search order
    # if l_start.x < l_end.x:
    #   order = 1
    # else:
    #   order = -1
  
  # Compute all circle-line intersections
  layers = layer.enumerateLayers()
  
  ipts = []
  for l in layers:
    x = getIntersections(m, b, l.x, l.y, l.outerRadius)
    if x[0] == x[1] or np.imag(x[0]) > 0:
      continue
    else:
      ipts.append({'layer': l, 'x':x[0], 'y':x[0]*m+b})
      ipts.append({'layer': l, 'x':x[1], 'y':x[1]*m+b})
      # ipts.append([x[0], m*x[0]+b])
      # ipts.append([x[1], m*x[1]+b])

  
  def kfunc(val):
    return val['x']
  
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
        # Now leaving outer surface of first conductor, reset edge list
        edges = []
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
      p0 = np.array(prev_edge['x'], prev_edge['y'])
      p1 = np.array(edge['x'], edge['y'])
      dist = np.linalg.norm(p0-p1)
      materials.append({'length': dist, 'layer' : current_layer})
      current_layer = current_layer.parent
    else: # We're entering a layer
      p0 = np.array(prev_edge['x'], prev_edge['y'])
      p1 = np.array(edge['x'], edge['y'])
      dist = np.linalg.norm(p0-p1)
      materials.append({'length': dist, 'layer' : current_layer})  
      current_layer = edge['layer']
    prev_edge = edge
  
  strength = 0
  air = databases.get_material_entry('air')
  for mat in materials:
    if mat['length'] == 0:
      continue

    breakdown_strength = mat['layer'].material_entry['dielectric_strength'].values[0]
    if mat['layer'].fillRatio < 1:
      breakdown_strength = air['dielectric_strength'].values[0]
    
    # print("Length: %f, Strength: %f" % (mat['length'], breakdown_strength))
    strength += breakdown_strength*mat['length']
  
  
  return strength, materials


def getIntersections(m, b, x0, y0, r):
  #y = mx + b
  #(x-x0)^2 + (y-y0)^2 = r^2
  #(x-x0)^2 + (m x + b -y0)^2 = r^2
  #x^2 - 2*x0*x + x0^2 + m^2 x^2 + 2*(b-y0)*x + (b-y0)^2 = r^2
  #(1+m^2) x^2 + (2*m*(b-y0) - 2*x0 )*x + (b-y0)^2 - r^2 + x0^2 = 0

  p = [(1+m**2), (2*m*(b-y0) - 2*x0), (b-y0)**2 - r**2 + x0**2]
  r = np.roots(p)
  
  # print(np.array(r)*m+b)
  
  return r
  