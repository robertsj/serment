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

  # Create an index list that omits the negative entries
  indices = -np.ones(nz*ny*nx, 'i')
  nodal_map_vec = vec_int(int(number_nodes), 0)
  count = 0
  for k in range(0, nz) :
    for j in range(0, ny) :
      for i in range(0, nx) :
        if nodal_map[k][j][i] < 0 : continue
        indices[index(i,j,k,nx,ny)] = count
        nodal_map_vec[count] = int(nodal_map[k][j][i])
        count += 1

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
          continue
        count += 1
        tmp_neigh = vec_neighbor(2*D, NeighborSurface(Node.VACUUM, 0))

        if i == 0 :
          tmp_neigh[CartesianNode.WEST] = bc[CartesianNode.WEST]
        else : # my west is their east
          if indices[index(i-1, j, k, nx, ny)] < 0 : 
            tmp_neigh[CartesianNode.WEST] = bc[CartesianNode.WEST]
          else :
            tmp_neigh[CartesianNode.WEST] = NeighborSurface(int(indices[index(i-1, j, k, nx, ny)]), CartesianNode.EAST)
        
        if i == nx - 1:
          tmp_neigh[CartesianNode.EAST] = bc[CartesianNode.EAST]
        else : # my east is their west
          if indices[index(i+1, j, k, nx, ny)] < 0 :
            tmp_neigh[CartesianNode.EAST] = bc[CartesianNode.EAST]
          else :
            tmp_neigh[CartesianNode.EAST] = NeighborSurface(int(indices[index(i+1, j, k, nx, ny)]), CartesianNode.WEST)

        if (D > 1) :
          if j == 0 :
            tmp_neigh[CartesianNode.SOUTH] = bc[CartesianNode.SOUTH]
          else : # my south is their north
            if indices[index(i, j-1, k, nx, ny)] < 0 :
              tmp_neigh[CartesianNode.SOUTH] = bc[CartesianNode.SOUTH]
            else :
              tmp_neigh[CartesianNode.SOUTH] = NeighborSurface(int(indices[index(i, j-1, k, nx, ny)]), CartesianNode.NORTH)
          if j == ny - 1:
            tmp_neigh[CartesianNode.NORTH] = bc[CartesianNode.NORTH]
          else : # my north is their south
            if indices[index(i, j+1, k, nx, ny)] < 0 :
              tmp_neigh[CartesianNode.NORTH] = bc[CartesianNode.NORTH]
            else :
              tmp_neigh[CartesianNode.NORTH] = NeighborSurface(int(indices[index(i, j+1, k, nx, ny)]), CartesianNode.SOUTH)

        if (D > 2) :
          if k == 0 :
            tmp_neigh[CartesianNode.BOTTOM] = bc[CartesianNode.BOTTOM]
          else : # my bottom is their top
            if indices[index(i, j, k-1, nx, ny)] < 0 :
              tmp_neigh[CartesianNode.BOTTOM] = bc[CartesianNode.BOTTOM]
            else :
              tmp_neigh[CartesianNode.BOTTOM] = NeighborSurface(int(indices[index(i, j, k-1, nx, ny)]), CartesianNode.TOP)
          if k == nz - 1:
            tmp_neigh[CartesianNode.TOP] = bc[CartesianNode.TOP]
          else : # my top is their bottom
            if indices[index(i, j, k+1, nx, ny)] < 0 :
              tmp_neigh[CartesianNode.TOP] = bc[CartesianNode.TOP]
            else :
              tmp_neigh[CartesianNode.TOP] = NeighborSurface(int(indices[index(i, j, k+1, nx, ny)]), CartesianNode.BOTTOM)

        neighbors.push_back(tmp_neigh)

  assert(count == number_nodes)
  origins = compute_origins(nodes, nodal_map, nodal_map_vec)  
  print origins.size()
  assert(origins.size() == number_nodes)

  node_list.set_nodal_map(nodal_map_vec, neighbors, origins)
  P = NodePartitioner()
  P.partition(node_list)

  return node_list

def compute_origins(nodes, nodal_map, nodal_map_vec) :
  """ Ensure that nodal widths are consistent.  Then, compute 
      origins of each node.
  """
  nx = np.size(nodal_map, 2)
  ny = np.size(nodal_map, 1)
  nz = np.size(nodal_map, 0)
  
  number_nodes = sum(nodal_map.flatten() >= 0)
  origins = vec_point()

  WZ = -np.ones(nz)
  WY = -np.ones(ny)
  WX = -np.ones(nx)
  for k in range(0, nz) :
    for j in range(0, ny) :
      for i in range(0, nx) :
        if nodal_map[k][j][i] < 0 : continue
        W = as_cartesian_node(nodes[nodal_map[k][j][i]]).width(2)
        if (WZ[k] < 0) : WZ[k] = W
        assert(WZ[k] == W)
  for j in range(0, ny) :
    for k in range(0, nz) :
      for i in range(0, nx) :
        if nodal_map[k][j][i] < 0 : continue
        W = as_cartesian_node(nodes[nodal_map[k][j][i]]).width(1)
        if (WY[j] < 0) : WY[j] = W
        assert(WY[j] == W)
  for i in range(0, nx) :
    for j in range(0, ny) :
      for k in range(0, nz) :
        if nodal_map[k][j][i] < 0 : continue
        W = as_cartesian_node(nodes[nodal_map[k][j][i]]).width(0)
        if (WX[i] < 0) : WX[i] = W
        assert(WX[i] == W)
  c = 0
  Z = 0.0
  for k in range(0, nz) :
    assert(WZ[k] > 0.0)
    Y = 0.0
    for j in range(0, ny) :
      assert(WY[j] > 0.0)
      X = 0.0
      for i in range(0, nx) : 
        assert(WX[i] > 0.0)
        if nodal_map[k][j][i] >= 0 : 
          #print X, Y, Z, nodal_map[k][j][i], nodal_map_vec[c]
          origins.push_back(Point(X, Y, Z))
          c += 1
  
        X += WX[i]
      Y += WY[j]
    Z += WZ[k]

  return origins

def index(i, j, k, nx, ny) :
  return i + j * nx + k * nx * ny

def plot_slice(nodes, dims=[1.0, 1.0, 1.0], n=100) :
  """ Plot a 2-D slice of a geometry.
  
  Currently, this is limited to 2-D problems.  An extension
  to 3-D should be straightforward, especially if limited 
  to x, y, and z planes.
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
  
    
  


