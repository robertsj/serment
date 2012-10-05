//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_geometry.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran geometry.
 */
//---------------------------------------------------------------------------//

// Need to include this in order to get appropriate typemaps
//%include "detran_utilities.i"

%include "geometry/Mesh.hh"
//%include "Mesh1D.hh"
//%include "Mesh2D.hh"
//%include "Mesh3D.hh"
%include "geometry/PinCell.hh"
%include "geometry/Assembly.hh"
%include "geometry/Core.hh"
//
%include "geometry/MeshMOC.hh"
//%include "Segment.hh"
//%include "Track.hh"
%include "geometry/TrackDB.hh"
%include "geometry/Tracker.hh"

namespace detran_geometry
{

// Dummy interface to avoid Boost serialization.
class Mesh1D: public Mesh
{
public:  
  Mesh1D(detran_utilities::vec_int xfm, 
         detran_utilities::vec_dbl xcme, 
         detran_utilities::vec_int mat_map);
  Mesh1D(detran_utilities::vec_dbl xfme, 
         detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_int xfm, 
         detran_utilities::vec_dbl xcme, 
         detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_dbl xfme, 
         detran_utilities::vec_int mat_map);
};

// Dummy interface to avoid Boost serialization.
class Mesh2D: public Mesh
{
public:
  Mesh2D(detran_utilities::vec_int xfm, 
         detran_utilities::vec_int yfm, 
         detran_utilities::vec_dbl xcme, 
         detran_utilities::vec_dbl ycme, 
         detran_utilities::vec_int mat_map);
  Mesh2D(detran_utilities::vec_dbl xfme, detran_utilities::vec_dbl yfme, detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_int xfm, 
         detran_utilities::vec_int yfm, 
         detran_utilities::vec_dbl xcme, 
         detran_utilities::vec_dbl ycme, 
         detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_dbl xfme, detran_utilities::vec_dbl yfme, detran_utilities::vec_int mat_map);
};

// Dummy interface to avoid Boost serialization.
class Mesh3D: public Mesh
{
public:
  Mesh3D(detran_utilities::vec_int xfm,  detran_utilities::vec_int yfm,  detran_utilities::vec_int zfm,
         detran_utilities::vec_dbl xcme, detran_utilities::vec_dbl ycme, detran_utilities::vec_dbl zcme,
         detran_utilities::vec_int mat_map);
  Mesh3D(detran_utilities::vec_dbl xfme, detran_utilities::vec_dbl yfme, detran_utilities::vec_dbl zfme, detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_int xfm,  detran_utilities::vec_int yfm,  detran_utilities::vec_int zfm,
         detran_utilities::vec_dbl xcme, detran_utilities::vec_dbl ycme, detran_utilities::vec_dbl zcme,
         detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_dbl xfme, detran_utilities::vec_dbl yfme, detran_utilities::vec_dbl zfme, detran_utilities::vec_int mat_map);
};

class Segment
{
public:
  typedef detran_utilities::SP<Segment> SP_segment;
  typedef detran_utilities::size_t      size_t;
  Segment(const size_t r, double l);
  int region() const;
  double length() const;
  void scale(double v);
};

class Track
{
public:
  typedef detran_utilities::SP<Track>         SP_track;
  typedef Segment::SP_segment                 SP_segment;
  typedef detran_utilities::Point             Point;
  typedef std::vector<Segment>                vec_segment;
  typedef vec_segment::const_iterator         iterator;
  typedef vec_segment::const_reverse_iterator riterator;
  typedef detran_utilities::size_t            size_t;
  Track(Point r0, Point r1);
  static SP_track Create(Point r0, Point r1);
  Point enter() const;
  Point exit() const;
  int number_segments() const;
  void add_segment(Segment s);
  Segment& segment(size_t i);
  const Segment& segment(size_t i) const;
  iterator begin();
  riterator rbegin();
  double cos_phi() const;
  double sin_phi() const;
};

} // end namespace detran_geometry

%template(MeshSP)     detran_utilities::SP<detran_geometry::Mesh>;
%template(Mesh1DSP)   detran_utilities::SP<detran_geometry::Mesh1D>;
%template(Mesh2DSP)   detran_utilities::SP<detran_geometry::Mesh2D>;
%template(Mesh3DSP)   detran_utilities::SP<detran_geometry::Mesh3D>;
%template(PinCellSP)  detran_utilities::SP<detran_geometry::PinCell>;
%template(AssemblySP) detran_utilities::SP<detran_geometry::Assembly>;
%template(CoreSP)     detran_utilities::SP<detran_geometry::Core>;

%template(SegmentSP)  detran_utilities::SP<detran_geometry::Tracker>;
%template(TrackSP)    detran_utilities::SP<detran_geometry::Track>;
%template(TrackDBSP)  detran_utilities::SP<detran_geometry::TrackDB>;
%template(TrackerSP)  detran_utilities::SP<detran_geometry::Tracker>;


// Save old versions for later post process (plot) use
//----------------------------------------------------------------------------//
//class Mesh1D : public Mesh
//{
//public:
//  Mesh1D(std::vector<int>    xfm, 
//         std::vector<double> xcme, 
//         std::vector<int>    mat_map);
//  Mesh1D(std::vector<double> xfme, 
//         std::vector<int>    mat_map);
//  static SP<Mesh> Create(std::vector<int>    xfm, 
//                         std::vector<double> xcme, 
//                         std::vector<int>    mat_map);
//  static SP<Mesh> Create(std::vector<double> xfme, 
//                         std::vector<int>    mat_map); 
//  void add_coarse_mesh_map(std::string map_key, std::vector<int> mesh_map);
//  void add_mesh_map(std::string map_key, std::vector<int> mesh_map);
//  %pythoncode 
//  %{
//      
//  def check(self) :
//    """ Check for matplotlib and numpy.  Call at top of all plotting methods.
//    """  
//    try :
//      import matplotlib.pyplot as plt
//      import numpy as np
//      self.PLOT = True
//    except ImportError :
//      self.PLOT = False
//      
//  def plot_mesh_map(self, map_key) :
//    """ Plot the mesh map.
//    """
//    self.check()  
//    if self.PLOT :
//      import matplotlib.pyplot as plt
//      import numpy as np
//      # Get the map.
//      m_map = np.asarray(self.mesh_map(map_key))
//      x = self.mesh_axes()
//      plt.plot(x, m_map)
//      plt.title('mesh map')
//      plt.show()
//    else :
//      print "Warning: matplotlib or numpy not found; skipping plot."
//      
//  def plot_flux(self, f) :
//    """ Plot a flux or flux-shaped vector on the mesh."
//    """
//    import matplotlib.pyplot as plt
//    import numpy as np
//    x = self.mesh_axes()
//    plt.plot(x, f)
//    plt.title('flux map')
//    plt.show()
//  
//  def mesh_axes(self) :
//    """ Get the fine mesh esges defining cells for plotting.
//    """
//    import matplotlib.pyplot as plt
//    import numpy as np
//    x = np.zeros(self.number_cells_x())
//    x[0] = 0.5*self.dx(0)
//    for i in range(1, self.number_cells_x()) :
//      x[i] = x[i-1] + self.dx(i)
//    return x
// 
//  %}  
//  
//};
//
////----------------------------------------------------------------------------//
//class Mesh2D : public Mesh
//{
//public:
//  Mesh2D(std::vector<int>    xfm, 
//         std::vector<int>    yfm, 
//         std::vector<double> xcme, 
//         std::vector<double> ycme, 
//         std::vector<int>    mat_map);
//  Mesh2D(std::vector<double> xfme, 
//         std::vector<double> yfme, 
//         std::vector<int>    mat_map);
//  static SP<Mesh> Create(std::vector<int>    xfm, 
//                         std::vector<int>    yfm, 
//                         std::vector<double> xcme, 
//                         std::vector<double> ycme, 
//                         std::vector<int>    mat_map);
//  static SP<Mesh> Create(std::vector<double> xfme, 
//                         std::vector<double> yfme, 
//                         std::vector<int>    mat_map); 
//  void add_coarse_mesh_map(std::string map_key, std::vector<int> mesh_map);
//  void add_mesh_map(std::string map_key, std::vector<int> mesh_map);
//  
//  %pythoncode 
//  %{
//    
//  def check(self) :
//    """ Check for matplotlib and numpy.  Call at top of all plotting methods.
//    """  
//    try :
//      import matplotlib.pyplot as plt
//      import numpy as np
//      self.PLOT = True
//    except ImportError :
//      self.PLOT = False
//      
//  def plot_mesh_map(self, map_key, edges=False) :
//    """ Plot the mesh map.
//    """
//    self.check()  
//    if self.PLOT :
//    
//      # Imports
//      import matplotlib
//      import matplotlib.pyplot as plt
//      import numpy as np
//      
//      # Get the map.
//      m_map = np.asarray(self.mesh_map(map_key))
//      if len(m_map) == 0 :
//        print "Oops! Map has size zero!"
//        return
//            
//      # Reshape it.
//      m_map2 = m_map.reshape(self.number_cells_x(), self.number_cells_y())
//      
//      # Define random color map.
//      colormap = matplotlib.colors.ListedColormap ( np.random.rand ( 10000,3))
//      
//      # Get the mesh axes and then make a grid of them for plotting.
//      x, y = self.mesh_axes()
//      X, Y = np.meshgrid(x, y)
//      
//      if edges == False :
//        edges = None
//      else :
//        edges = "black"
//
//      # Create the plot a
//      plt.pcolor(X, Y, m_map2, cmap=colormap)
//      plt.show()
//      
//    else :
//      print "Warning: matplotlib or numpy not found; skipping plot."
//      
//  def plot_flux(self, f) :
//    """ Plot a flux or flux-shaped vector on the mesh."
//    """
//    import matplotlib.pyplot as plt
//    import numpy as np
//    # Reshape flux.
//    f = f.reshape(self.number_cells_x(), self.number_cells_y())
//    #print f
//    # Plot me. 
//    plt.imshow(f, interpolation='nearest', cmap=plt.cm.hot)
//    plt.title('flux map')
//    plt.colorbar()
//    plt.show()
//  
//  def mesh_axes(self) :
//    """ Get the fine mesh esges defining cells for plotting.
//    """
//    import matplotlib.pyplot as plt
//    import numpy as np
//    x = np.zeros(self.number_cells_x()+1)
//    y = np.zeros(self.number_cells_y()+1)
//    for i in range(0, self.number_cells_x()) :
//      x[i + 1] = x[i] + self.dx(i)
//    for j in range(0, self.number_cells_y()) :
//      y[j + 1] = y[j] + self.dy(j)
//    print y
//    return (x, y)
//         
//  %}
//};
  

//} // end namespace detran

//---------------------------------------------------------------------------//
//              end of detran_geometry.i
//---------------------------------------------------------------------------//





