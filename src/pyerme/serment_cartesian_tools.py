# Tools for building Cartesian reactor models

from detran import *
from serment import *

def build_model(nodes, nodal_map, global_bc) :

  # Ensure the map has valid entries and dimension
  nodal_map = np.array(nodal_map)
  assert(nodal_map.max() < len(nodes))
  assert(nodal_map.min() >= -2)
  D = len(nodal_map.shape)
  assert(D >= 1)
  assert(D <= 3)
  for node in nodes :
    assert(node.dimension() == D)
    assert(node.number_surfaces() == 2 * D)

  # Reshape to 1-D, 2-D, or 3-D
  number_nodes = sum(nodal_map.flatten() >= 0)
  nx = 1
  ny = 1
  nz = 1
  if D == 1 :
    nx = nodal_map.shape[0]
  elif D == 2 :
    nx = nodal_map.shape[1]
    ny = nodal_map.shape[0]
  else :
    nx = nodal_map.shape[2]
    ny = nodal_map.shape[1]
    nz = nodal_map.shape[0]
  nodal_map = nodal_map.reshape((nz, ny, nx))

  # Set the boundary conditions
  assert(len(global_bc) == 2 * D)
  bc = []
  for side in range(0, 2 * D) :
    if global_bc[side] == Node.VACUUM :
      bc.append(NeighborSurface(Node.VACUUM, 0))
    elif global_bc[side] == Node.REFLECT :
      bc.append(NeighborSurface(Node.REFLECT, 0))
    else :
      exit('ERROR: Invalid global boundary condition')

  # Create the node and neighbor lists
  node_list = NodeList.Create()
  for node in nodes :
    node_list.add_node(node)
  neighbors = vec2_neighbor()

  # Fill out the neighbor list
  count = 0
  for k in range(0, nz) :
    for j in range(0, ny) :
      for i in range(0, nx) :        
        # neighbors
        neighbor_index = [index(i-1, j  , k  , nx, ny), \
                          index(i+1, j  , k  , nx, ny), \
                          index(i  , j-1, k  , nx, ny), \
                          index(i  , j+1, k  , nx, ny), \
                          index(i  , j  , k-1, nx, ny), \
                          index(i  , j  , k+1, nx, ny)  ]  
        # this should be enough to allow jagged boundaries
        node_index = nodal_map[k][j][i]
        if node_index < 0 :
          break
        count += 1
        tmp_neigh = vec_neighbor(2*D, NeighborSurface(Node.VACUUM, 0))

        if i == 0 :
          tmp_neigh[CartesianNode.WEST] = bc[CartesianNode.WEST]
        else : # my west is their east
          tmp_neigh[CartesianNode.WEST] = NeighborSurface(neighbor_index[0], CartesianNode.EAST)
        if i == nx - 1:
          tmp_neigh[CartesianNode.EAST] = bc[CartesianNode.EAST]
        else : # my east is their west
          tmp_neigh[CartesianNode.EAST] = NeighborSurface(neighbor_index[1], CartesianNode.WEST)

        if (D > 1) :
          if j == 0 :
            tmp_neigh[CartesianNode.SOUTH] = bc[CartesianNode.SOUTH]
          else : # my south is their north
            tmp_neigh[CartesianNode.SOUTH] = NeighborSurface(neighbor_index[2], CartesianNode.NORTH)
          if j == ny - 1:
            tmp_neigh[CartesianNode.NORTH] = bc[CartesianNode.NORTH]
          else : # my north is their south
            tmp_neigh[CartesianNode.NORTH] = NeighborSurface(neighbor_index[3], CartesianNode.SOUTH)

        if (D > 2) :
          if k == 0 :
            tmp_neigh[CartesianNode.BOTTOM] = bc[CartesianNode.BOTTOM]
          else : # my bottom is their top
            tmp_neigh[CartesianNode.BOTTOM] = NeighborSurface(neighbor_index[4], CartesianNode.TOP)
          if k == nz - 1:
            tmp_neigh[CartesianNode.TOP] = bc[CartesianNode.TOP]
          else : # my top is their bottom
            tmp_neigh[CartesianNode.TOP] = NeighborSurface(neighbor_index[5], CartesianNode.BOTTOM)

        neighbors.push_back(tmp_neigh)

  assert(count == number_nodes)
  origins = compute_origins(nodes, nodal_map)  

  nodal_map_vec = vec_int(int(number_nodes), 0)
  nodal_map = nodal_map.reshape(number_nodes)
  for i in range(0, number_nodes) :
    nodal_map_vec[i] = int(nodal_map[i])

  node_list.set_nodal_map(nodal_map_vec, neighbors, origins)
  P = NodePartitioner()
  P.partition(node_list)

  return node_list

def compute_origins(nodes, nodal_map) :
  """ Ensure that nodal widths are consistent.  Then, compute 
      origins of each node.
  """
  nx = np.size(nodal_map, 2)
  ny = np.size(nodal_map, 1)
  nz = np.size(nodal_map, 0)

  number_nodes = sum(nodal_map.flatten() >= 0)
  origins = vec_point()
  n = 0
  Z = 0.0
  for k in range(0, nz) :
    Y = 0.0
    for j in range(0, ny) :
      X = 0.0 
      for i in range(0, nx) : 
        if nodal_map[k][j][i] < 0 : break
        # this makes sure each node in a plane (XY, YZ, or XZ) has 
        # identical measurement in the transverse direction
        node = as_cartesian_node(nodes[nodal_map[k][j][i]])
        if i == 0 : WX = node.width(0)
        if j == 0 : WY = node.width(1)
        if k == 0 : WZ = node.width(2)
        assert(node.width(0) == WX)
        assert(node.width(1) == WY)
        assert(node.width(2) == WZ)
        # add the origin 
        origins.push_back(Point(X, Y, Z))
        X += WX
        # increment the node index
        n += 1
      Y += WY
    Z += WZ
  print "ORIGIN", origins[0]
  return origins

def index(i, j, k, nx, ny) :
  return i + j * nx + k * nx * ny

def plot_slice(nodes, dims=[1.0, 1.0, 1.0], n=100) :
  """ Plot a 2-D slice of a geometry.
  """ 

  assert(nodes.is_finalized())

  D = nodes.node(0).dimension()

  if D == 1 :
    x = np.linspace(0, dims[0], n)
  
  elif D == 2 :

    plotter = PPMPlotter()
    plotter.initialize(n, n, "serment.ppm")
    
    x = np.linspace(0, dims[0], n)
    y = np.linspace(0, dims[1], n)
    X, Y = np.meshgrid(x, y)
    c = -1.0*np.ones(n**2)
    for j in range(0, n) :
      for i in range(0, n) :
        p1 = Point(X[j][i], Y[j][i], 0.0)   
        for k in range(0, nodes.number_global_nodes()) :
          p2 = p1 - nodes.origin(k)
          if nodes.node(k).color(p2) >= 0.0 :
            c[i+j*n] = nodes.node(k).color(p2)
            plotter.set_pixel(i, n-j-1, c[i+j*n])
            break
    plotter.write()
  
    
  


