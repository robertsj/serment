
// File: index.xml

// File: classerme_1_1AbsorptionOperator.xml
%feature("docstring") erme::AbsorptionOperator "

Converts a global moments vector into a global absorption rate.

C++ includes: AbsorptionOperator.hh ";

%feature("docstring")  erme::AbsorptionOperator::AbsorptionOperator "home robertsj Research serment source src erme AbsorptionOperator cc
erme::AbsorptionOperator::AbsorptionOperator(SP_nodelist nodes,
SP_indexer indexer, SP_server server)

Constructor.

Parameters:
-----------

nodes:  Pointer to node list

indexer:  Pointer to response indexer

server:  Pointer to response server ";

%feature("docstring")  erme::AbsorptionOperator::update "void
erme::AbsorptionOperator::update()

Update the vector data. ";


// File: classAbsorptionResponse.xml
%feature("docstring") AbsorptionResponse "

This class encapsulates a SermentVector and required RF updates.

C++ includes: AbsorptionResponse.hh ";

%feature("docstring")  AbsorptionResponse::AbsorptionResponse "AbsorptionResponse::AbsorptionResponse(GlobalInput &input,
ResponseFunctionServer *s) ";

%feature("docstring")  AbsorptionResponse::~AbsorptionResponse "AbsorptionResponse::~AbsorptionResponse() ";

%feature("docstring")  AbsorptionResponse::updateData "void
AbsorptionResponse::updateData(scalar k)

This function updates underlying matrix elements. ";


// File: classBase.xml
%feature("docstring") Base "

eigenvalue response matrix solver ";


// File: classerme__geometry_1_1CartesianNode.xml
%feature("docstring") erme_geometry::CartesianNode "

Base Cartesian node class.

C++ includes: CartesianNode.hh ";

%feature("docstring")  erme_geometry::CartesianNode::CartesianNode "home robertsj Research serment source src erme_geometry CartesianNode
cc home robertsj Research serment source src erme_geometry
CartesianNode cc erme_geometry::CartesianNode::CartesianNode(const
size_t dimension, std::string nodename, vec2_size_t so, vec_size_t po,
vec_size_t ao, vec_size_t eo, vec_dbl nodewidth) ";

%feature("docstring")  erme_geometry::CartesianNode::area "double
erme_geometry::CartesianNode::area(const size_t surface) const

Return the area of a node surface.

Parameters:
-----------

surface:  Surface index ";

%feature("docstring")  erme_geometry::CartesianNode::volume "double
erme_geometry::CartesianNode::volume() const

Return the volume of the node. ";

%feature("docstring")  erme_geometry::CartesianNode::color "virtual
double erme_geometry::CartesianNode::color(Point point)

Default color. Am I in the box or not? ";


// File: classerme__geometry_1_1CartesianNodeDetran.xml
%feature("docstring") erme_geometry::CartesianNodeDetran "

Specialization of CartesianNode for use with Detran.

The node essentially defines the local problem to be solved by Detran.
The required input consists of the parameters database, the material
definitions, and the mesh. Any other needed objects are created later
in the response generator.

C++ includes: CartesianNodeDetran.hh ";

%feature("docstring")
erme_geometry::CartesianNodeDetran::CartesianNodeDetran "erme_geometry::CartesianNodeDetran::CartesianNodeDetran(const size_t
dimension, std::string nodename, vec2_size_t so, vec_size_t po,
vec_size_t ao, vec_size_t eo, vec_dbl nodewidth, SP_db nodedb,
SP_material nodematerial, SP_mesh nodemesh)

Constructor.

Parameters:
-----------

dimension:  Dimension of the node

name:  Name

so:  Spatial orders [surface][axis]

po:  Polar orders [surface]

ao:  Azimuth orders [surface]

eo:  Energy orders [surface] ";

%feature("docstring")  erme_geometry::CartesianNodeDetran::color "double erme_geometry::CartesianNodeDetran::color(Point point)

Default color. Am I in the box or not? ";

%feature("docstring")  erme_geometry::CartesianNodeDetran::db "SP_db
erme_geometry::CartesianNodeDetran::db() const ";

%feature("docstring")  erme_geometry::CartesianNodeDetran::material "SP_material erme_geometry::CartesianNodeDetran::material() const ";

%feature("docstring")  erme_geometry::CartesianNodeDetran::mesh "SP_mesh erme_geometry::CartesianNodeDetran::mesh() const ";


// File: classerme__geometry_1_1CartesianNodeDummy.xml
%feature("docstring") erme_geometry::CartesianNodeDummy "

Dummy Cartesian node for testing

C++ includes: DummyNode.hh ";

%feature("docstring")
erme_geometry::CartesianNodeDummy::CartesianNodeDummy "erme_geometry::CartesianNodeDummy::CartesianNodeDummy(const size_t
dim, const size_t so=0, const size_t po=0, const size_t ao=0, const
size_t eo=0) ";

%feature("docstring")  erme_geometry::CartesianNodeDummy::color "double erme_geometry::CartesianNodeDummy::color(Point point)

Default color. Am I in the box or not? ";


// File: classserment__comm_1_1Comm.xml
%feature("docstring") serment_comm::Comm "

Parallel communication interface.

This is an easy API for using MPI (or some other parallel library).

Todo Consider implementing Comm as a singleton so that better
encapsulation can be used

C++ includes: Comm.hh ";


// File: structserment__comm_1_1Comm__Traits.xml
%feature("docstring") serment_comm::Comm_Traits "

This struct and its specializations are used to implement the type-
safe default message tags in Comm. Any other type-determined property
needed in Comm would also go here.

C++ includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01char_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< char * > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01char_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< char > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01double_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< double * > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01double_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< double > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01float_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< float * > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01float_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< float > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01int_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< int * > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01int_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< int > " C++ includes:
Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01long_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< long * > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01long_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< long > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01long_01double_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< long double * > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01long_01double_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< long double > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01short_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< short * > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01short_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< short > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01unsigned_01char_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< unsigned char * > "
C++ includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01unsigned_01char_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< unsigned char > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01unsigned_01int_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< unsigned int * > "
C++ includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01unsigned_01int_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< unsigned int > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01unsigned_01long_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< unsigned long * > "
C++ includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01unsigned_01long_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< unsigned long > " C++
includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01unsigned_01short_01_5_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< unsigned short * > "
C++ includes: Comm_Traits.hh ";


// File: structserment__comm_1_1Comm__Traits_3_01unsigned_01short_01_4.xml
%feature("docstring") serment_comm::Comm_Traits< unsigned short > "
C++ includes: Comm_Traits.hh ";


// File: classerme_1_1Connect.xml
%feature("docstring") erme::Connect "

Defines the geometric relationship between nodes.

In the simplest case of one unknown per node surface, the connectivity
matrix is simply an adjacency matrix in which a 1 is placed at any row
and column location representing the surface of one node adjoining the
surface of another node.

Vacuum or reflective boundary conditions are built into the
connectivity matrix. Reflective conditions require a polarity shift
for odd moments; this applies to some combinations of space and angle
basis functions. Rotational symmetry can be built into the node
connections.

C++ includes: Connect.hh ";

%feature("docstring")  erme::Connect::Connect "home robertsj Research
serment source src erme Connect cc erme::Connect::Connect(SP_nodelist
nodes, SP_indexer indexer)

Constructor.

Parameters:
-----------

nodes:  Pointer to list of nodes

indexer:  Pointer to response indexer ";


// File: classDataBase.xml
%feature("docstring") DataBase "

To be completed.

C++ includes: DataBase.hh ";

%feature("docstring")  DataBase::DataBase "DataBase::DataBase() ";

%feature("docstring")  DataBase::~DataBase "DataBase::~DataBase() ";


// File: structerme__response_1_1ResponseDatabase_1_1DBResponse.xml
%feature("docstring") erme_response::ResponseDatabase::DBResponse "C++ includes: ResponseDatabase.hh ";


// File: classDiff2dElement.xml
%feature("docstring") Diff2dElement "

A class for representing the domain of an individual diff2d problem.

Diff2dElement is a class that represents geometry and other domain-
specific aspects of a problem. A Diff2dProblem consists of one or
several such elements. This arrangement helps streamline response
function generation.

C++ includes: Diff2dElement.hh ";

%feature("docstring")  Diff2dElement::Diff2dElement "Diff2dElement::Diff2dElement() ";

%feature("docstring")  Diff2dElement::~Diff2dElement "Diff2dElement::~Diff2dElement() ";

%feature("docstring")  Diff2dElement::destroy "void
Diff2dElement::destroy() ";


// File: classDiff2dInput.xml
%feature("docstring") Diff2dInput "

A class for handling and containing data for diff2d problems.

Diff2dInput is a class that performs functions necessary for both
stand- alone problems using diff2d, a simple 2-D diffusion code, and
problems where diff2d is used to generate response functions.

C++ includes: Diff2dInput.hh ";

%feature("docstring")  Diff2dInput::readInput "void
Diff2dInput::readInput(char *file)

readInput

This function reads the local problem input and initializes all data
subsequently used to make the operater matrices and solve the problem.

The input file structure, and consequently, the reader, should be made
more general in future refinements. However, the basic data is all
that is needed for now. ";

%feature("docstring")  Diff2dInput::skipStuff "void
Diff2dInput::skipStuff(ifstream &in)

This function skips over blocks of comments and blank lines. ";


// File: classDiff2dOutput.xml
%feature("docstring") Diff2dOutput "

A class for outputing fluxes and other things.

Diff2dOutput is blah blah blah.

C++ includes: Diff2dOutput.hh ";

%feature("docstring")  Diff2dOutput::Diff2dOutput "Diff2dOutput::Diff2dOutput(Diff2dInput &inp) ";

%feature("docstring")  Diff2dOutput::~Diff2dOutput "Diff2dOutput::~Diff2dOutput() ";

%feature("docstring")  Diff2dOutput::doOutput "void
Diff2dOutput::doOutput(Diff2dProblem &prob, Diff2dSolver &sol) ";

%feature("docstring")  Diff2dOutput::closeSilo "void
Diff2dOutput::closeSilo() ";


// File: classDiff2dProblem.xml
%feature("docstring") Diff2dProblem "

A class containing operators (i.e. matrices) for diff2d problems.

To be completed.

C++ includes: Diff2dProblem.hh ";

%feature("docstring")  Diff2dProblem::Diff2dProblem "Diff2dProblem::Diff2dProblem() ";

%feature("docstring")  Diff2dProblem::Diff2dProblem "Diff2dProblem::Diff2dProblem(Diff2dInput &input, integer elid) ";

%feature("docstring")  Diff2dProblem::~Diff2dProblem "Diff2dProblem::~Diff2dProblem() ";

%feature("docstring")  Diff2dProblem::updateRHS "void
Diff2dProblem::updateRHS(scalar k, int O, int G, int S)

Update the RHS with new order, group, and side. ";

%feature("docstring")  Diff2dProblem::getKeff "scalar
Diff2dProblem::getKeff() ";

%feature("docstring")  Diff2dProblem::destroy "void
Diff2dProblem::destroy()

Destroy the K matrices ... finish me! ";

%feature("docstring")  Diff2dProblem::setK "void
Diff2dProblem::setK()

This function actually sets the K matrices. ";

%feature("docstring")  Diff2dProblem::setRHS "void
Diff2dProblem::setRHS()

This function sets things on the RHS affected by updates.

The purpose of a separate fixedSource evaluator is to separate source
terms updated during rf-generation from K and other quantities that
are static for a given geometry/material definition. This only
includes the fixedSource---all other quantities are static. By
separating, we can avoid the cost of re-constructing things. For now,
I'm setting it up only for response function generation, i.e. I'm not
accounting for stand alone problems with a boundary source. ";

%feature("docstring")  Diff2dProblem::setValues "void
Diff2dProblem::setValues(scalar A[], integer idx[], integer k, integer
lenk, integer g, scalar dv, scalar cix_i, scalar ciy_j, integer m_C)

Set values of the matrix. ";

%feature("docstring")  Diff2dProblem::resizeR "void
Diff2dProblem::resizeR(integer m, scalar kk) ";

%feature("docstring")  Diff2dProblem::clearR "void
Diff2dProblem::clearR(integer m) ";


// File: classDiff2dSolver.xml
%feature("docstring") Diff2dSolver "

A class for solving diff2d problems.

Diff2dSolver is blah blah blah.

C++ includes: Diff2dSolver.hh ";

%feature("docstring")  Diff2dSolver::Diff2dSolver "Diff2dSolver::Diff2dSolver() ";

%feature("docstring")  Diff2dSolver::Diff2dSolver "Diff2dSolver::Diff2dSolver(Diff2dInput &input, Diff2dProblem &problem,
int el) ";

%feature("docstring")  Diff2dSolver::~Diff2dSolver "Diff2dSolver::~Diff2dSolver() ";

%feature("docstring")  Diff2dSolver::solve "void
Diff2dSolver::solve()

This function drives the solution for a single problem. ";

%feature("docstring")  Diff2dSolver::compRespFct "void
Diff2dSolver::compRespFct()

This routine computes the response functions.

The current responses are computed. The absorption, leakage, and
fission responses are computed. Add some math here, maybe. Remember,
odd orders are reversed for the incident side and for its right hand
neighbor.

For the current responses, note that the incident basis functions are
always oriented with respect to the incident face. That means for the
top and left edges, the basis functions are reversed when the right
hand side is constructed. That is accounted for here where the alpha
terms are indexed in reverse for those edges. Similarly, the right and
bottom current responses are constructed in reverse. This is because
the responses must be expanded with respect to their face \"looking
outward\" into the next element. Since top and left are already
reversed, we expand as normal, but for right and bottom, we must
explicitly account for that. ";

%feature("docstring")  Diff2dSolver::destroy "void
Diff2dSolver::destroy()

Destroy the K matrices ... finish me! ";


// File: classlinear__algebra_1_1EigenSolver.xml
%feature("docstring") linear_algebra::EigenSolver "C++ includes:
EigenSolver.hh ";

%feature("docstring")  linear_algebra::EigenSolver::EigenSolver "home
robertsj Research serment source src linear_algebra EigenSolver cc
linear_algebra::EigenSolver::EigenSolver(SP_matrix A, SP_matrix
B=SP_matrix(0))

Constructor.

Parameters:
-----------

A:  Pointer to left hand operator

P:  Pointer to right hand operator (possibly null) ";

%feature("docstring")  linear_algebra::EigenSolver::solve "double
linear_algebra::EigenSolver::solve(SP_vector x)

Solve the eigenvalue problem.

Parameters:
-----------

x:  Eigenvector

Dominant eigenvalue ";


// File: classerme__solver_1_1EigenvalueUpdate.xml
%feature("docstring") erme_solver::EigenvalueUpdate "C++ includes:
EigenvalueUpdate.hh ";

%feature("docstring")  erme_solver::EigenvalueUpdate::EigenvalueUpdate
"erme_solver::EigenvalueUpdate::EigenvalueUpdate()

Constructor. ";

%feature("docstring")
erme_solver::EigenvalueUpdate::~EigenvalueUpdate "virtual
erme_solver::EigenvalueUpdate::~EigenvalueUpdate()

Virtual destructor. ";

%feature("docstring")  erme_solver::EigenvalueUpdate::compute "virtual double erme_solver::EigenvalueUpdate::compute(const double
keff, SP_vector J)

Computes and updated keff, possibly based on previous history.

Parameters:
-----------

keff:  latest eigenvalue estimate

J:  latest boundary unknowns

Default does nothing ";


// File: classerme_1_1FissionOperator.xml
%feature("docstring") erme::FissionOperator "

Converts a global moments vector into a global fission rate.

C++ includes: FissionOperator.hh ";

%feature("docstring")  erme::FissionOperator::FissionOperator "home
robertsj Research serment source src erme FissionOperator cc
erme::FissionOperator::FissionOperator(SP_nodelist nodes, SP_indexer
indexer, SP_server server)

Constructor.

Parameters:
-----------

nodes:  Pointer to node list

indexer:  Pointer to response indexer

server:  Pointer to response server ";

%feature("docstring")  erme::FissionOperator::update "void
erme::FissionOperator::update()

Update the vector data. ";


// File: classFissionResponse.xml
%feature("docstring") FissionResponse "

This class encapsulates a SermentVector and required RF updates.

C++ includes: FissionResponse.hh ";

%feature("docstring")  FissionResponse::FissionResponse "FissionResponse::FissionResponse(GlobalInput &input,
ResponseFunctionServer *s) ";

%feature("docstring")  FissionResponse::~FissionResponse "FissionResponse::~FissionResponse() ";

%feature("docstring")  FissionResponse::updateData "void
FissionResponse::updateData(scalar k)

This function updates underlying matrix elements. ";


// File: classGlobalSolver.xml
%feature("docstring") GlobalSolver "

A base class for all global solvers.

The eigenvalue

C++ includes: GlobalSolver.hh ";

%feature("docstring")  GlobalSolver::GlobalSolver "GlobalSolver::GlobalSolver(SP_globalproblem problem, SP_globalinput
input)

Constructor. ";

%feature("docstring")  GlobalSolver::~GlobalSolver "GlobalSolver::~GlobalSolver()

Pure virtual destructor. ";

%feature("docstring")  GlobalSolver::solve "virtual void
GlobalSolver::solve()=0

Solve the eigenvalue response matrix equations. ";


// File: classerme__solver_1_1GlobalSolverBase.xml
%feature("docstring") erme_solver::GlobalSolverBase "C++ includes:
GlobalSolverBase.hh ";

%feature("docstring")  erme_solver::GlobalSolverBase::GlobalSolverBase
"home robertsj Research serment source src erme_solver
GlobalSolverBase cc
erme_solver::GlobalSolverBase::GlobalSolverBase(SP_db db, SP_indexer
indexer, SP_server server, SP_state state, SP_R R, SP_M M, SP_F F,
SP_A A, SP_L L)

Constructor.

Parameters:
-----------

db:  Pointer to parameter database

indexer:  Pointer to response indexer

server:  Pointer to response server

state:  Pointer to state vector

R:  Pointer to response matrix

M:  Pointer to connectivity matrix

F:  Pointer to fission operator

A:  Pointer to absorption operator

L:  Pointer to leakage operator ";

%feature("docstring")
erme_solver::GlobalSolverBase::~GlobalSolverBase "erme_solver::GlobalSolverBase::~GlobalSolverBase()=0

Pure virtual destructor. ";

%feature("docstring")  erme_solver::GlobalSolverBase::solve "virtual
void erme_solver::GlobalSolverBase::solve()=0

Solve. ";


// File: classerme__solver_1_1GlobalSolverNewton.xml
%feature("docstring") erme_solver::GlobalSolverNewton "

Solves the problem using Picard (fixed point) iteration.

C++ includes: GlobalSolverNewton.hh ";

%feature("docstring")
erme_solver::GlobalSolverNewton::~GlobalSolverNewton "virtual
erme_solver::GlobalSolverNewton::~GlobalSolverNewton()

Virtual destructor. ";

%feature("docstring")  erme_solver::GlobalSolverNewton::solve "void
erme_solver::GlobalSolverNewton::solve()

Solve. ";


// File: classerme__solver_1_1GlobalSolverPicard.xml
%feature("docstring") erme_solver::GlobalSolverPicard "

Solves the problem using Picard (fixed point) iteration.

The eigenvalue response matrix equations can be cast in the form of an
inner current eigenvalue equation \\\\[
\\\\mathbf{MR}(k^{n})\\\\mathbf{J} = \\\\lambda\\\\mathbf{J} \\\\, ,
\\\\] with the associated $ k $ eigenvalue update \\\\[ k^{n+1} =
\\\\frac{ \\\\mathbf{F}(k^{n}) } { \\\\mathbf{A}(k^{n}) +
mathbf{L}(k^{n}) } \\\\] which is a mathematical statement of gains-
to-losses. These coupled equations represent a fixed-point iteration
in the nonlinear variable $ k $. Steffensen's method is easily
implemented by extrapolating from three successive $ k $ values and is
available as an optional update function.

C++ includes: GlobalSolverPicard.hh ";

%feature("docstring")
erme_solver::GlobalSolverPicard::GlobalSolverPicard "home robertsj
Research serment source src erme_solver GlobalSolverPicard cc home
robertsj Research serment source src erme_solver GlobalSolverPicard cc
erme_solver::GlobalSolverPicard::GlobalSolverPicard(SP_db db,
SP_indexer indexer, SP_server server, SP_state state, SP_R R, SP_M M,
SP_F F, SP_A A, SP_L L)

Constructor.

Parameters:
-----------

server:  Pointer to response server

state:  Pointer to state vector

R:  Pointer to response matrix

M:  Pointer to connectivity matrix

F:  Pointer to fission operator

A:  Pointer to absorption operator

L:  Pointer to leakage operator ";

%feature("docstring")
erme_solver::GlobalSolverPicard::~GlobalSolverPicard "virtual
erme_solver::GlobalSolverPicard::~GlobalSolverPicard()

Virtual destructor. ";

%feature("docstring")  erme_solver::GlobalSolverPicard::solve "void
erme_solver::GlobalSolverPicard::solve()

Solve. ";


// File: classerme__geometry_1_1HexagonalNode.xml
%feature("docstring") erme_geometry::HexagonalNode "

Base Hexagonal node class.

C++ includes: HexagonalNode.hh ";

%feature("docstring")  erme_geometry::HexagonalNode::HexagonalNode "erme_geometry::HexagonalNode::HexagonalNode(const size_t dimension,
std::string nodename=\"hexagonal_node\", vec2_size_t so, vec_size_t
po, vec_size_t ao, vec_size_t eo, const double face_width, const
double height=1.0) ";

%feature("docstring")  erme_geometry::HexagonalNode::area "double
erme_geometry::HexagonalNode::area(const size_t surface) const

Return the area of a node surface.

Parameters:
-----------

surface:  Surface index ";

%feature("docstring")  erme_geometry::HexagonalNode::volume "double
erme_geometry::HexagonalNode::volume() const

Return the volume of the node. ";

%feature("docstring")  erme_geometry::HexagonalNode::color "virtual
double erme_geometry::HexagonalNode::color(Point point)=0

Return the color associated with the spatial coordinate. ";


// File: classInnerIterBase.xml
%feature("docstring") InnerIterBase "

Base class for inner iterations within an outer power-like iteration.

C++ includes: InnerIterBase.hh ";

%feature("docstring")  InnerIterBase::InnerIterBase "InnerIterBase::InnerIterBase(SP_M M, SP_R R)

Constructor. ";

%feature("docstring")  InnerIterBase::~InnerIterBase "virtual
InnerIterBase::~InnerIterBase() ";

%feature("docstring")  InnerIterBase::solve "virtual scalar
InnerIterBase::solve(int max_iters, double tol, SP_vector J_in,
SP_vector J)=0

Perform inner iterations.

Parameters:
-----------

max_iterations:  Maximum number of iterations, nominally meaning
applications of MR on a vector.

tolerance:  Generic tolerance on convergence.

J_0:  Initial guess.

J:  Solved eigenvector.

Dominant eigenvalue. ";


// File: classInnerIterPower.xml
%feature("docstring") InnerIterPower "

Performs inner iterations via power method within outer iteration.

This is our own implementation of the power method for inners within
power-like outers. The SLEPc interface also provides the power method,
and their implementation is to be preferred.

The power (iteration) method is a standard procedure for finding the
largest eigenvalue of an operator.

C++ includes: InnerIterPower.hh ";

%feature("docstring")  InnerIterPower::InnerIterPower "InnerIterPower::InnerIterPower(SP_M M, SP_R R)

Constructor. ";

%feature("docstring")  InnerIterPower::~InnerIterPower "InnerIterPower::~InnerIterPower()

Destructor. ";

%feature("docstring")  InnerIterPower::solve "scalar
InnerIterPower::solve(int max_iters, double tol, SP_vector J_in,
SP_vector J)

Perform inner iterations.

Parameters:
-----------

max_iterations:  Maximum number of iterations, nominally meaning
applications of MR on a vector.

tolerance:  Generic tolerance on convergence.

X_0:  Initial guess. ";


// File: classInnerIterSLEPc.xml
%feature("docstring") InnerIterSLEPc "

Performs inner iterations using SLEPc within outer iteration.

This is an interface class for use of SLEPc eigensolvers for inner
iterations. SLEPc offers several algorithms which can be further
supplemented by additional libraries (e.g. ARPACK).

C++ includes: InnerIterSLEPc.hh ";

%feature("docstring")  InnerIterSLEPc::InnerIterSLEPc "InnerIterSLEPc::InnerIterSLEPc(SP_M M, SP_R R)

Constructor. ";

%feature("docstring")  InnerIterSLEPc::~InnerIterSLEPc "InnerIterSLEPc::~InnerIterSLEPc()

Destructor.

Frees the SLEPc items. ";

%feature("docstring")  InnerIterSLEPc::solve "scalar
InnerIterSLEPc::solve(int max_iters, double tol, SP_vector J_in,
SP_vector J)

Perform inner iterations.

Parameters:
-----------

max_iterations:  Maximum number of iterations, nominally meaning
applications of MR on a vector.

tolerance:  Generic tolerance on convergence.

X_0:  Initial guess. ";


// File: classInputXML.xml
%feature("docstring") InputXML "C++ includes: InputXML.hh ";

%feature("docstring")  InputXML::InputXML "InputXML::InputXML(string
S) ";

%feature("docstring")  InputXML::~InputXML "InputXML::~InputXML() ";

%feature("docstring")  InputXML::readInput "bool
InputXML::readInput(char *filename)

This opens, verifies, and processes an xml input file.

The opening and verification are meant to be General enough for use in
all input scenarios needing a \"load\" and \"verify\" sequence. The
processing is of course program-specific, and will be left as a
virtual method in this abstract class.

Parameters:
-----------

filename:  the xml input file ";


// File: classInvItShell.xml
%feature("docstring") InvItShell "C++ includes: InvItShell.hh ";

%feature("docstring")  InvItShell::InvItShell "InvItShell::InvItShell(integer a, integer b, void *ctx, GlobalProblem
*pr) ";

%feature("docstring")  InvItShell::~InvItShell "InvItShell::~InvItShell() ";

%feature("docstring")  InvItShell::myMatVec "void
InvItShell::myMatVec(Mat &M, Vec &x, Vec &y)

This performs the action of (M*R-lambda*I) for use in inverse it. ";

%feature("docstring")  InvItShell::updateEigs "void
InvItShell::updateEigs(scalar kval, scalar lval) ";


// File: classerme__solver_1_1Jacobian.xml
%feature("docstring") erme_solver::Jacobian "

This performs the action of the Jacobian ( $ \\\\mathbf{f}' $ )on a
vector $ \\\\mathbf{f} $.

The resulting vector is $ \\\\mathbf{f}'\\\\mathbf{f} $, or \"fpf\"
where \"p\" denotes the prime.

Note, the vector $ \\\\vec{f} $ does not necessarily begin as the
residual but rather as a scaled residual. The PETSc manual does not
seem to state this explicitly, but a little trial and error confirmed
this---so don't fret when debugging!

The Jacobian is defined as \\\\[ \\\\mathbf{f'(x)} = \\\\left
[\\\\begin{array}{ccc} (\\\\mathbf{M}\\\\mathbf{R}-\\\\lambda
\\\\mathbf{I}) & \\\\mathbf{M}\\\\mathbf{R_k}\\\\mathbf{J_-} &
\\\\mathbf{J_-} \\\\\\\\ (\\\\mathbf{F}-k\\\\mathbf{L}) &
(\\\\mathbf{F_k}-k\\\\mathbf{L_k}-\\\\mathbf{L}) \\\\mathbf{J_-} & 0
\\\\\\\\ \\\\mathbf{J^T_-} & 0 & 0 \\\\end{array} \\\\right ] \\\\, .
\\\\label{eq:jacobian} \\\\] For $ \\\\mathbf{R}(k) $ of size $
m\\\\times m $, the Jacobian is of size $ (m+2)\\\\times(m+2) $.

Most of the Jacobian is known  a priori, and so most of the action can
be computed directly. Only the first $ m-1 $ rows of the $ (m-1) $th
column require approximations via finite differences.

C++ includes: Jacobian.hh ";


// File: classJacobianShell.xml
%feature("docstring") JacobianShell "C++ includes: JacobianShell.hh
";

%feature("docstring")  JacobianShell::JacobianShell "JacobianShell::JacobianShell(integer a, integer b, void *ctx,
SermentVector *unk, SermentVector *res, GlobalProblem *pr, SNES &sn)
";

%feature("docstring")  JacobianShell::~JacobianShell "JacobianShell::~JacobianShell() ";

%feature("docstring")  JacobianShell::myMatVec "void
JacobianShell::myMatVec(Mat &M, Vec &f, Vec &fpf)

This performs the action of the Jacobian ( $ \\\\mathbf{f}' $ )on a
vector $ \\\\mathbf{f} $.

The resulting vector is $ \\\\mathbf{f}'\\\\mathbf{f} $, or \"fpf\"
where \"p\" denotes the prime.

Note, the vector $ \\\\vec{f} $ does not necessarily begin as the
residual but rather as a scaled residual. The PETSc manual does not
seem to state this explicitly, but a little trial and error confirmed
this---so don't fret when debugging!

The Jacobian is defined as \\\\[ \\\\mathbf{f'(x)} = \\\\left
[\\\\begin{array}{ccc} (\\\\mathbf{M}\\\\mathbf{R}-\\\\lambda
\\\\mathbf{I}) & \\\\mathbf{M}\\\\mathbf{R_k}\\\\mathbf{J_-} &
\\\\mathbf{J_-} \\\\\\\\ (\\\\mathbf{F}-k\\\\mathbf{L}) &
(\\\\mathbf{F_k}-k\\\\mathbf{L_k}-\\\\mathbf{L}) \\\\mathbf{J_-} & 0
\\\\\\\\ \\\\mathbf{J^T_-} & 0 & 0 \\\\end{array} \\\\right ] \\\\, .
\\\\label{eq:jacobian} \\\\] For $ \\\\mathbf{R}(k) $ of size $
m\\\\times m $, the Jacobian is of size $ (m+2)\\\\times(m+2) $.

Most of the Jacobian is known  a priori, and so most of the action can
be computed directly. Only the first $ m-1 $ rows of the $ (m-1) $th
column require approximations via finite differences. ";


// File: classerme_1_1LeakageOperator.xml
%feature("docstring") erme::LeakageOperator "

Leakage operator.

C++ includes: LeakageOperator.hh ";

%feature("docstring")  erme::LeakageOperator::LeakageOperator "home
robertsj Research serment source src erme LeakageOperator cc
erme::LeakageOperator::LeakageOperator(SP_nodelist nodes, SP_indexer
indexer, SP_server server)

Constructor.

Parameters:
-----------

nodes:  Pointer to node list

indexer:  Pointer to response indexer

server:  Pointer to response server ";

%feature("docstring")  erme::LeakageOperator::update "void
erme::LeakageOperator::update()

Update the response matrix data. ";

%feature("docstring")  erme::LeakageOperator::leakage "double
erme::LeakageOperator::leakage(linear_algebra::Vector &x)

Compute the net global leakage given a global moments vector.

Parameters:
-----------

x:  Global moments vector ";

%feature("docstring")  erme::LeakageOperator::display_leakage "void
erme::LeakageOperator::display_leakage()

Display the global leakage vector. ";


// File: classLeakageResponse.xml
%feature("docstring") LeakageResponse "

This is a matrix operator for computing the leakage response.

Because it is a matrix quantity similar to the response matrix, we
build from the same common base and use a BCRS matrix format. Note,
the block size is the number of edges/faces. For example, a 2d
diffusion element has a leakage response function comprised of four
vectors, one for each side. The inner product of these vectors with
the incident current vector (incident on that element) produce the net
leakage from the element.

C++ includes: LeakageResponse.hh ";

%feature("docstring")  LeakageResponse::LeakageResponse "LeakageResponse::LeakageResponse(GlobalInput &input,
ResponseFunctionServer *s, integer *idx, integer mmnum) ";

%feature("docstring")  LeakageResponse::~LeakageResponse "LeakageResponse::~LeakageResponse() ";

%feature("docstring")  LeakageResponse::updateData "void
LeakageResponse::updateData(scalar k)

This function updates underlying matrix elements. ";

%feature("docstring")  LeakageResponse::computeLeakage "scalar
LeakageResponse::computeLeakage(SermentVector &J)

This computes the leakage using the newest response function data. ";

%feature("docstring")  LeakageResponse::getLeakageVec "void
LeakageResponse::getLeakageVec(SermentVector &leakV) ";


// File: classLegendrePoly.xml
%feature("docstring") LegendrePoly "

This class is a complete Legendre polynomial basis set for a given
element.

LegendrePoly computes and stores the (discrete) Legendre polynomials
needed to produce incident boundary conditions and to expand outgoing
partial currents as response functions. The polynomial vectors are
produced once at initialization and stored for subsequent application.

C++ includes: LegendrePoly.hh ";

%feature("docstring")  LegendrePoly::LegendrePoly "LegendrePoly::LegendrePoly() ";

%feature("docstring")  LegendrePoly::~LegendrePoly "LegendrePoly::~LegendrePoly() ";

%feature("docstring")  LegendrePoly::buildMe "void
LegendrePoly::buildMe(integer nx, integer ny, integer maxOrder)

This creates the discrete Legendre polynomials for all edges. ";

%feature("docstring")  LegendrePoly::makeWeights "void
LegendrePoly::makeWeights()

This function creates the normalization weights. ";

%feature("docstring")  LegendrePoly::expandCur "void
LegendrePoly::expandCur(Vec &partialCurrent, scalar legCoefs[])

This function expands a partial current and gives back coefficients.

The user gives the partial current Petsc Vec and an array to be filled
with the coefficients. Note, the array *must* be the correct size,
i.e. order+1, as there is currently no check. Example use: scalar
coefsX[order+1], LP.expandCur( PcurX, coefsX ); ";

%feature("docstring")  LegendrePoly::factorial "scalar
LegendrePoly::factorial(integer x)

This is a simple factorial function. ";

%feature("docstring")  LegendrePoly::destroy "void
LegendrePoly::destroy()

This deallocates the Petsc Vecs explicitly. ";


// File: classlinear__algebra_1_1LinearSolver.xml
%feature("docstring") linear_algebra::LinearSolver "C++ includes:
LinearSolver.hh ";

%feature("docstring")  linear_algebra::LinearSolver::LinearSolver "home robertsj Research serment source src linear_algebra LinearSolver
cc linear_algebra::LinearSolver::LinearSolver(SP_matrix A, SP_matrix
P)

Constructor.

Parameters:
-----------

A:  Pointer to linear system matrix

P:  Pointer to preconditioning matrix, possibly equal to A ";

%feature("docstring")  linear_algebra::LinearSolver::solve "void
linear_algebra::LinearSolver::solve(SP_vector b, SP_vector x)

Solve the linear system.

Parameters:
-----------

b:  right hand side

x:  solution ";


// File: classLocalProblem.xml
%feature("docstring") LocalProblem "

This class encapsulates the entire local problem routine.

to be completed

C++ includes: LocalProblem.hh ";

%feature("docstring")  LocalProblem::LocalProblem "LocalProblem::LocalProblem(string in) ";

%feature("docstring")  LocalProblem::~LocalProblem "LocalProblem::~LocalProblem() ";

%feature("docstring")  LocalProblem::getResponseFunctions "virtual
ResponseFunction** LocalProblem::getResponseFunctions(scalar keff)=0
";

%feature("docstring")  LocalProblem::localTime "scalar
LocalProblem::localTime() ";


// File: classLocalProblemDiff2d.xml
%feature("docstring") LocalProblemDiff2d "

This class encapsulates the entire Diff2d local problem routine.

to be completed

C++ includes: LocalProblemDiff2d.hh ";

%feature("docstring")  LocalProblemDiff2d::LocalProblemDiff2d "LocalProblemDiff2d::LocalProblemDiff2d(string in) ";

%feature("docstring")  LocalProblemDiff2d::~LocalProblemDiff2d "LocalProblemDiff2d::~LocalProblemDiff2d() ";

%feature("docstring")  LocalProblemDiff2d::getResponseFunctions "ResponseFunction ** LocalProblemDiff2d::getResponseFunctions(scalar
keff)

This function returns an array of ResponseFunctionDiffusion objects.

to be completed ";


// File: classerme__utils_1_1ManagerERME.xml
%feature("docstring") erme_utils::ManagerERME "

Manages solution of eigenvalue response matrix problem.

C++ includes: ManagerERME.hh ";

%feature("docstring")  erme_utils::ManagerERME::ManagerERME "home
robertsj Research serment source src erme_utils ManagerERME cc
erme_utils::ManagerERME::ManagerERME(int argc, char *argv[], SP_db db)

Constructor. ";

%feature("docstring")  erme_utils::ManagerERME::build_erme "void
erme_utils::ManagerERME::build_erme(SP_nodelist nodes)

Construct the problem.

Parameters:
-----------

nodes:  List of nodes ";

%feature("docstring")  erme_utils::ManagerERME::solve "void
erme_utils::ManagerERME::solve()

Solve the problem.

The solver is built on-the-fly, which allows one to quickly change
parameters via the input database for successive solves. ";

%feature("docstring")  erme_utils::ManagerERME::indexer "SP_indexer
erme_utils::ManagerERME::indexer() const

Return the indexer. ";

%feature("docstring")  erme_utils::ManagerERME::finalize "void
erme_utils::ManagerERME::finalize()

Close libraries, etc. ";


// File: classlinear__algebra_1_1Matrix.xml
%feature("docstring") linear_algebra::Matrix "C++ includes: Matrix.hh
";

%feature("docstring")  linear_algebra::Matrix::Matrix "home robertsj
Research serment source src linear_algebra Matrix cc
linear_algebra::Matrix::Matrix(const size_type m, const size_type n,
const vec_int &number_nonzeros, const vec_int
&number_nonzeros_offdiagonal=vec_int(0)) ";


// File: classlinear__algebra_1_1MatrixBase.xml
%feature("docstring") linear_algebra::MatrixBase "C++ includes:
MatrixBase.hh ";

%feature("docstring")  linear_algebra::MatrixBase::MatrixBase "linear_algebra::MatrixBase::MatrixBase(const size_type m, const
size_type n)

Constructor.

Parameters:
-----------

m:  local number of rows

n:  local number of columns

number_nonzeros:

number_nonzeros_offdiagonal:  ";

%feature("docstring")  linear_algebra::MatrixBase::~MatrixBase "linear_algebra::MatrixBase::~MatrixBase()=0

Virtual destructor. ";

%feature("docstring")  linear_algebra::MatrixBase::insert_values "void linear_algebra::MatrixBase::insert_values(const size_type
number_rows, const int *rows, const size_type number_columns, const
int *columns, const double *values)

Insert values.

Parameters:
-----------

number_rows:  Number of column indices

rows:  Indices of global rows

number_columns:  Number of column indices

columns:  Indices of global columns

values:  Logically 2D array of values to insert ";

%feature("docstring")  linear_algebra::MatrixBase::assemble "void
linear_algebra::MatrixBase::assemble()

Assemble the matrix. ";

%feature("docstring")  linear_algebra::MatrixBase::multiply "void
linear_algebra::MatrixBase::multiply(Vector &x, Vector &y)

Matrix-vector multiplication.

Parameters:
-----------

x:  Input vector

y:  Output vector ";

%feature("docstring")  linear_algebra::MatrixBase::multiply_transpose
"void linear_algebra::MatrixBase::multiply_transpose(Vector &x,
Vector &y)

Matrix-vector multiplication using matrix transpose.

Parameters:
-----------

x:  Input vector

y:  Output vector ";

%feature("docstring")  linear_algebra::MatrixBase::A "Mat
linear_algebra::MatrixBase::A()

Get PETSc Mat object. ";

%feature("docstring")  linear_algebra::MatrixBase::number_global_rows
"size_type linear_algebra::MatrixBase::number_global_rows() const ";

%feature("docstring")
linear_algebra::MatrixBase::number_global_columns "size_type
linear_algebra::MatrixBase::number_global_columns() const ";

%feature("docstring")  linear_algebra::MatrixBase::number_local_rows "size_type linear_algebra::MatrixBase::number_local_rows() const ";

%feature("docstring")
linear_algebra::MatrixBase::number_local_columns "size_type
linear_algebra::MatrixBase::number_local_columns() const ";

%feature("docstring")  linear_algebra::MatrixBase::lower_bound "int
linear_algebra::MatrixBase::lower_bound() const

Return the lower bound. ";

%feature("docstring")  linear_algebra::MatrixBase::upper_bound "int
linear_algebra::MatrixBase::upper_bound() const

Return the upper bound. ";

%feature("docstring")  linear_algebra::MatrixBase::is_assembled "bool
linear_algebra::MatrixBase::is_assembled() const ";

%feature("docstring")  linear_algebra::MatrixBase::display "void
linear_algebra::MatrixBase::display(const int output=0, const
std::string name=\"matrix.out\") const

Display the matrix to screen (or to output)

Parameters:
-----------

output:  Flag indicating (stdout=0, ascii=1, binary=2)

name:  File name for ascii or binary file ";


// File: classlinear__algebra_1_1MatrixShell.xml
%feature("docstring") linear_algebra::MatrixShell "

Shell matrix class for user storage and operator scheme.

To use this class, a user must inherit it and implement the relevant
shell methods, e.g. shell_multiply. Unfortunately, inside these
methods the user must work entirely with the PETSc API.

C++ includes: MatrixShell.hh ";

%feature("docstring")  linear_algebra::MatrixShell::MatrixShell "home
robertsj Research serment source src linear_algebra MatrixShell cc
linear_algebra::MatrixShell::MatrixShell(const size_type m, const
size_type n, void *context) ";

%feature("docstring")  linear_algebra::MatrixShell::insert_values "virtual void linear_algebra::MatrixShell::insert_values(const
size_type number_rows, const int *rows, const size_type
number_columns, const int *columns, const double *values)

Insert values.

Parameters:
-----------

number_rows:  Number of column indices

rows:  Indices of global rows

number_columns:  Number of column indices

columns:  Indices of global columns

values:  Logically 2D array of values to insert ";

%feature("docstring")  linear_algebra::MatrixShell::assemble "virtual
void linear_algebra::MatrixShell::assemble()

Assemble the matrix. ";


// File: structserment__comm_1_1MPI__Traits.xml
%feature("docstring") serment_comm::MPI_Traits "

Provide a generic way to get MPI_Datatype arguments for MPI function
calls.

This struct provides a generic programming--common way to get
MPI_Datatype arguments for MPI function calls. The static function,
element_type(), returns an argument of type MPI_Datatype that matches
a C++ datatype with an MPI_Datatype.

C++ includes: MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01char_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< char > " C++ includes:
MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01double_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< double > " C++
includes: MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01float_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< float > " C++
includes: MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01int_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< int > " C++ includes:
MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01long_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< long > " C++ includes:
MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01long_01double_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< long double > " C++
includes: MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01short_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< short > " C++
includes: MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01unsigned_01char_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< unsigned char > " C++
includes: MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01unsigned_01int_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< unsigned int > " C++
includes: MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01unsigned_01long_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< unsigned long > " C++
includes: MPI_Traits.hh ";


// File: structserment__comm_1_1MPI__Traits_3_01unsigned_01short_01_4.xml
%feature("docstring") serment_comm::MPI_Traits< unsigned short > " C++
includes: MPI_Traits.hh ";


// File: classlinear__algebra_1_1MyMatrixShell.xml
%feature("docstring") linear_algebra::MyMatrixShell "C++ includes:
matrix_shell_fixture.hh ";

%feature("docstring")  linear_algebra::MyMatrixShell::MyMatrixShell "linear_algebra::MyMatrixShell::MyMatrixShell(const size_type m, const
size_type n) ";

%feature("docstring")  linear_algebra::MyMatrixShell::shell_multiply "PetscErrorCode linear_algebra::MyMatrixShell::shell_multiply(Vec x,
Vec y)

Matrix-vector multiplication.

Parameters:
-----------

x:  Input vector

y:  Output vector ";

%feature("docstring")
linear_algebra::MyMatrixShell::shell_multiply_transpose "PetscErrorCode
linear_algebra::MyMatrixShell::shell_multiply_transpose(Vec x, Vec y)

Matrix-vector multiplication using matrix transpose.

Parameters:
-----------

x:  Input vector

y:  Output vector ";


// File: classerme__geometry_1_1NeighborSurface.xml
%feature("docstring") erme_geometry::NeighborSurface "

Contains a neighbor index and the surface shared.

C++ includes: NeighborSurface.hh ";

%feature("docstring")  erme_geometry::NeighborSurface::NeighborSurface
"erme_geometry::NeighborSurface::NeighborSurface(const int n=0, const
size_type s=0) ";

%feature("docstring")  erme_geometry::NeighborSurface::neighbor "int
erme_geometry::NeighborSurface::neighbor() const ";

%feature("docstring")  erme_geometry::NeighborSurface::surface "size_type erme_geometry::NeighborSurface::surface() const ";


// File: classNewton.xml
%feature("docstring") Newton "

This class solves GlobalProblems via use of Newton's method.

to be completed

C++ includes: Newton.hh ";

%feature("docstring")  Newton::Newton "Newton::Newton(SP_globalproblem problem, SP_globalinput input) ";

%feature("docstring")  Newton::~Newton "Newton::~Newton() ";

%feature("docstring")  Newton::solve "void Newton::solve()

This is the Newton solver. ";

%feature("docstring")  Newton::Residual "PetscErrorCode
Newton::Residual(SNES snes, Vec X, Vec F, void *ptr)

This function computes the nonlinear residual.

The nonlinear residual is defined as \\\\[ \\\\mathbf{f(x)} = \\\\left
[\\\\begin{array}{c} (\\\\mathbf{M}\\\\mathbf{R}(k)-\\\\lambda
\\\\mathbf{I}) \\\\mathbf{J_-} \\\\\\\\
\\\\mathbf{F}(k)\\\\mathbf{J_-} - (k\\\\mathbf{L}(k)\\\\mathbf{J_-} )
\\\\\\\\ \\\\frac{1}{2} \\\\mathbf{J^T_-} \\\\mathbf{J_-} -
\\\\frac{1}{2} \\\\end{array} \\\\right ] = \\\\mathbf{0} \\\\, ,
\\\\] which is the same as used in the

Todo Residual should be a class so that the temp vectors don't have to
be rebuilt every time ";


// File: classerme__geometry_1_1Node.xml
%feature("docstring") erme_geometry::Node "

Base node class.

The response matrix method decomposes a system into independent
computational nodes. This class represents an abstract node, from
which nodes with specific geometries can be derived.

A Node represents everything needed by the underlying response
generator, be it Detran or something else. Note, the maximum response
order per variable is specified. However, the actual order realized
could depend on dynamic order schemes. Hence, while a unique node can
be specified once, the actual expansion orders used for multiple
instances of the node may vary.

Geometrically, nodes are rather flexible. In 1D, they can represent
slab slices, cylindrical shells, or spherical shells. In 2D, any
closed shape can be well-defined. In 3D, arbitrary right cylinders can
be defined. For 3D, bottom and top surfaces (in the xy plane) are
enforced to simplify polarity shifts in the polar angle upon
reflection at global boundaries. The bottom and top surface must be
the last two surfaces, respectively, in any 3D node.

Because visualization of a node (or the set of connected nodes) is
extremely useful, a Node can implement a color function that returns
some useful information about the node. Because all nodes are assigned
an origin, the implementation can determine whether a queried point is
within the node, and if so, return a color. By default, the south-
west-bottom corner of the global domain is the global origin, and the
s-w-b corner of nodes should be the local origin---different nodes
might need to be flexible with that notion.

C++ includes: Node.hh ";

%feature("docstring")  erme_geometry::Node::Node "erme_geometry::Node::Node(const size_t dimension, const size_t
number_surfaces, std::string name, vec2_size_t so, vec_size_t po,
vec_size_t ao, vec_size_t eo)

Constructor.

Parameters:
-----------

dimension:  Dimension of the node

number_surfaces:  Number of surfaces

name:  Unique Name

so:  Spatial orders [surface][axis]

po:  Polar orders [surface]

ao:  Azimuth orders [surface]

eo:  Energy orders [surface] ";

%feature("docstring")  erme_geometry::Node::~Node "erme_geometry::Node::~Node()=0

Pure virtual destructor. ";

%feature("docstring")  erme_geometry::Node::area "virtual double
erme_geometry::Node::area(const size_t surface) const =0

Return the area of a node surface.

Parameters:
-----------

surface:  Surface index ";

%feature("docstring")  erme_geometry::Node::volume "virtual double
erme_geometry::Node::volume() const =0

Return the volume of the node. ";

%feature("docstring")  erme_geometry::Node::color "virtual double
erme_geometry::Node::color(Point point)=0

Return a color. ";

%feature("docstring")  erme_geometry::Node::dimension "Node::size_t
erme_geometry::Node::dimension() const

Node dimension. ";

%feature("docstring")  erme_geometry::Node::number_surfaces "Node::size_t erme_geometry::Node::number_surfaces() const

Number nodal surfaces. ";

%feature("docstring")  erme_geometry::Node::name "std::string
erme_geometry::Node::name() const

Nodal name. ";

%feature("docstring")  erme_geometry::Node::spatial_order "Node::size_t erme_geometry::Node::spatial_order(const size_t s, const
size_t d) const

Spatial order for a surface and possible dimension. ";

%feature("docstring")  erme_geometry::Node::polar_order "Node::size_t
erme_geometry::Node::polar_order(const size_t s) const

Polar order for a surface. ";

%feature("docstring")  erme_geometry::Node::azimuthal_order "Node::size_t erme_geometry::Node::azimuthal_order(const size_t s)
const

Azimuthal order for a surface. ";

%feature("docstring")  erme_geometry::Node::energy_order "Node::size_t erme_geometry::Node::energy_order(const size_t s) const

Energy order for a surface. ";

%feature("docstring")  erme_geometry::Node::display "void
erme_geometry::Node::display() const

Pretty print of key characteristics. ";


// File: classerme__geometry_1_1NodeFactory.xml
%feature("docstring") erme_geometry::NodeFactory "

Base class for constructing various node types.

C++ includes: NodeFactory.hh ";

%feature("docstring")  erme_geometry::NodeFactory::NodeFactory "erme_geometry::NodeFactory::NodeFactory()

Constructor. ";

%feature("docstring")  erme_geometry::NodeFactory::~NodeFactory "virtual erme_geometry::NodeFactory::~NodeFactory()

Virtual destructor. ";

%feature("docstring")  erme_geometry::NodeFactory::create_node "virtual SP_node erme_geometry::NodeFactory::create_node(SP_db db,
SP_material material=SP_material(), SP_mesh mesh=SP_mesh())=0

Create a node.

Everything needed to build the node must exist in the database.

Parameters:
-----------

db:  Parameter database.

db:  ";


// File: classerme__geometry_1_1NodeFactoryDetran.xml
%feature("docstring") erme_geometry::NodeFactoryDetran "

Build Detran-based nodes.

C++ includes: NodeFactoryDetran.hh ";

%feature("docstring")
erme_geometry::NodeFactoryDetran::NodeFactoryDetran "erme_geometry::NodeFactoryDetran::NodeFactoryDetran()

Constructor. ";

%feature("docstring")
erme_geometry::NodeFactoryDetran::~NodeFactoryDetran "virtual
erme_geometry::NodeFactoryDetran::~NodeFactoryDetran()

Virtual destructor. ";

%feature("docstring")  erme_geometry::NodeFactoryDetran::create_node "home robertsj Research serment source src erme_geometry
NodeFactoryDetran cc Node::SP_node
erme_geometry::NodeFactoryDetran::create_node(SP_db db, SP_material
material, SP_mesh mesh)

Create a node.

Parameters:
-----------

db:  Parameter database. ";


// File: classerme__geometry_1_1NodeList.xml
%feature("docstring") erme_geometry::NodeList "

Container for nodes and indices of their neighbors.

The node list maintains a vector of unique nodes, their placement in
the global domain, and their neighbors. For now, a \"unique\" node is
defined by both its underlying transport problem *and* its requested
expansion orders. Hence, if the same node is used multiple times in a
problem with varying orders (e.g. for scoping adaptivity), then those
are unique.

A neighbor index must be supplied for each nodal surface. If the
surface is a global boundary, then either Node::REFLECT or
Node::VACUUM must be specified. The surfaces are in a Node-specific
ordering. For example, Cartesian nodes are indexed as follows: left,
right, bottom, top, south, north.

Once the user adds all nodes, finalize() must be called.

C++ includes: NodeList.hh ";

%feature("docstring")  erme_geometry::NodeList::NodeList "home
robertsj Research serment source src erme_geometry NodeList cc
erme_geometry::NodeList::NodeList()

Constructor. ";

%feature("docstring")  erme_geometry::NodeList::set_bounds "void
erme_geometry::NodeList::set_bounds(const size_t lb, const size_t ub)

Set the local node array bounds.

The entire vector of nodes lives on all processes. The partitioner
must set the bounds corresponding to a particular node. The unique
nodes on each process are also identified.

Parameters:
-----------

lb:  Lower bound

ub:  Upper bound ";

%feature("docstring")  erme_geometry::NodeList::add_node "void
erme_geometry::NodeList::add_node(SP_node n)

Add a unique node.

Parameters:
-----------

n:   Node to be added ";

%feature("docstring")  erme_geometry::NodeList::set_nodal_map "void
erme_geometry::NodeList::set_nodal_map(const vec_int &nodes, const
vec2_neighbor &neighbors)

Set the map of nodes with their neighbors.

Parameters:
-----------

nodal_indices:  Indices into the vector of unique nodes

neighbor:  Neighbor data for each node ";

%feature("docstring")  erme_geometry::NodeList::node "NodeList::SP_node erme_geometry::NodeList::node(const int n) const

Get a node.

Parameters:
-----------

n:  Index into vector of unique nodes via the global index ";

%feature("docstring")  erme_geometry::NodeList::unique_node "NodeList::SP_node erme_geometry::NodeList::unique_node(const int n)
const

Get a node.

Parameters:
-----------

n:  Index into vector of unique nodes via the cardinal index ";

%feature("docstring")  erme_geometry::NodeList::lower_bound "NodeList::size_t erme_geometry::NodeList::lower_bound() const

Get a node. Returns null pointer if not found.

Parameters:
-----------

name:  Name of node Return local lower bound ";

%feature("docstring")  erme_geometry::NodeList::upper_bound "NodeList::size_t erme_geometry::NodeList::upper_bound() const

Return local upper bound. ";

%feature("docstring")  erme_geometry::NodeList::number_global_nodes "NodeList::size_t erme_geometry::NodeList::number_global_nodes() const

Number of global nodes. ";

%feature("docstring")  erme_geometry::NodeList::number_local_nodes "NodeList::size_t erme_geometry::NodeList::number_local_nodes() const

Number of local nodes. ";

%feature("docstring")  erme_geometry::NodeList::number_global_surfaces
"NodeList::size_t erme_geometry::NodeList::number_global_surfaces()
const

Number of local surfaces. ";

%feature("docstring")  erme_geometry::NodeList::number_local_surfaces
"NodeList::size_t erme_geometry::NodeList::number_local_surfaces()
const

Number of local surfaces. ";

%feature("docstring")
erme_geometry::NodeList::number_unique_global_nodes "NodeList::size_t
erme_geometry::NodeList::number_unique_global_nodes() const

Number of unique global nodes. ";

%feature("docstring")
erme_geometry::NodeList::number_unique_local_nodes "NodeList::size_t
erme_geometry::NodeList::number_unique_local_nodes() const

Number of unique local nodes. ";

%feature("docstring")  erme_geometry::NodeList::neighbor "const
NeighborSurface & erme_geometry::NodeList::neighbor(const size_t
node_g, const size_t s) const

Get the global neighbor index for a node surface.

Parameters:
-----------

node_g:  Global node index

s:   Node surface ";

%feature("docstring")
erme_geometry::NodeList::global_index_from_local "NodeList::size_t
erme_geometry::NodeList::global_index_from_local(const size_t node_l)
const

Get the global index of a local node.

This indexes the nodal map that defines the geometry of the problem.

Parameters:
-----------

node_l:  Local node index ";

%feature("docstring")
erme_geometry::NodeList::local_index_from_global "int
erme_geometry::NodeList::local_index_from_global(const size_t node_g)
const

Get the local index of a global node.

Global index translated to the local portion of the nodal map. Returns
a negative value if the global index is not within the local range.

Parameters:
-----------

node_g:  Global node index ";

%feature("docstring")
erme_geometry::NodeList::unique_global_index_from_global "NodeList::size_t
erme_geometry::NodeList::unique_global_index_from_global(const size_t
node_g) const

Get the unique index of a global node.

This simply returns the nodal map entry. Hence, the entries in the map
must correspond to the nodes actually defined by the user.

Parameters:
-----------

n:  Global node index ";

%feature("docstring")
erme_geometry::NodeList::unique_local_index_from_unique_global "int
erme_geometry::NodeList::unique_local_index_from_unique_global(const
size_t node_ug) const

Get the unique local index from the unique global index.

Starting with the complete nodal map, the portion in this local group
is sorted and unique elements found. This indexer searches for the
location within that sorted segment for the given global identifier.
If not found, returns negative.

Parameters:
-----------

node_ug:  Unique global node index ";

%feature("docstring")
erme_geometry::NodeList::unique_global_index_from_unique_local "NodeList::size_t
erme_geometry::NodeList::unique_global_index_from_unique_local(const
size_t node_ul) const

Get the unique global index from the unique local index.

Parameters:
-----------

node_ul:  Unique local index ";

%feature("docstring")  erme_geometry::NodeList::is_finalized "bool
erme_geometry::NodeList::is_finalized() const

Have all nodes been added? ";

%feature("docstring")  erme_geometry::NodeList::display "void
erme_geometry::NodeList::display() const

Display all the nodes in the list. ";


// File: classerme__geometry_1_1NodePartitioner.xml
%feature("docstring") erme_geometry::NodePartitioner "

Partition nodes one level 1 communicator.

This is a very light weight partitioning that simply broadcasts the
list of nodes and assigns array bounds for each receiving process.

C++ includes: NodePartitioner.hh ";

%feature("docstring")  erme_geometry::NodePartitioner::NodePartitioner
"home robertsj Research serment source src erme_geometry
NodePartitioner cc erme_geometry::NodePartitioner::NodePartitioner()

Constructor. ";

%feature("docstring")  erme_geometry::NodePartitioner::partition "void erme_geometry::NodePartitioner::partition(SP_nodelist &nodes)

Partition a list of nodes.

Parameters:
-----------

nodes:   Node list ";


// File: classerme__response_1_1NodeResponse.xml
%feature("docstring") erme_response::NodeResponse "

Container for nodal response functions.

NodeResponse is a very simple container for responses. For the
simplest problems, the only responses required are for the boundary
function (usually the partial current or angular flux), fission,
absorption, and leakage. Suppose we have a node with $ S $ surfaces,
and on each surface, we have a partial current expanded in $ M $
terms. The total boundary vector for this node has a size of $ N = SM
$. The total size of the boundary response data is then $ N^2 $, since
each of the $ N $ moments contributes to outgoing currents that are
also expanded in $ N $ terms. The fission and absorption operators
represent vectors with which the node incident boundary function is
folded to yield total fission and absorption rates; that vector has a
size of $ N $. The leakage operator yields the total leakage from each
of $ S $ when operated on the node boundary function; hence, the node
leakage data has a size of $ SN $.

To be efficient, we want these responses stored contiguously if
possible. If we view the node boundary responses as an $ N \\\\times N
$ block, then one column is produced per incident response. Hence, we
are best served using column-oriented storage for this data. Moreover,
this facilitates the probable case of response server processes
sending these columns back to a response driver process, after which
the driver participates in the global solve.

In this initial implementation, we'll use STL vectors. For the
boundary responses, columns will be stored contiguously.

C++ includes: NodeResponse.hh ";

%feature("docstring")  erme_response::NodeResponse::NodeResponse "home robertsj Research serment source src erme_response NodeResponse
cc erme_response::NodeResponse::NodeResponse(const size_t N, const
size_t number_surfaces)

Constructor.

Parameters:
-----------

moments_size:  Size of node moments vector for each surface ";

%feature("docstring")  erme_response::NodeResponse::boundary_response
"const double & erme_response::NodeResponse::boundary_response(const
size_t out, const size_t in) const

Const access to boundary response. ";

%feature("docstring")  erme_response::NodeResponse::boundary_response
"double & erme_response::NodeResponse::boundary_response(const size_t
out, const size_t in)

Mutable access to boundary response. ";

%feature("docstring")  erme_response::NodeResponse::fission_response "const double & erme_response::NodeResponse::fission_response(const
size_t in) const

Const access to fission response. ";

%feature("docstring")  erme_response::NodeResponse::fission_response "double & erme_response::NodeResponse::fission_response(const size_t
in)

Mutable access to fission response. ";

%feature("docstring")
erme_response::NodeResponse::absorption_response "const double &
erme_response::NodeResponse::absorption_response(const size_t in)
const

Const access to absorption response. ";

%feature("docstring")
erme_response::NodeResponse::absorption_response "double &
erme_response::NodeResponse::absorption_response(const size_t in)

Mutable access to absorption response. ";

%feature("docstring")  erme_response::NodeResponse::leakage_response "const double & erme_response::NodeResponse::leakage_response(const
size_t surface, const size_t in) const

Const access to leakage response. ";

%feature("docstring")  erme_response::NodeResponse::leakage_response "double & erme_response::NodeResponse::leakage_response(const size_t
surface, const size_t in)

Mutable access to leakage response. ";

%feature("docstring")  erme_response::NodeResponse::size "size_t
erme_response::NodeResponse::size() const

Return moment size. ";

%feature("docstring")  erme_response::NodeResponse::number_surfaces "size_t erme_response::NodeResponse::number_surfaces() const

Return number of surfaces. ";

%feature("docstring")  erme_response::NodeResponse::display "void
erme_response::NodeResponse::display() const

Display the response data. ";


// File: classerme__solver_1_1NonlinearResidual.xml
%feature("docstring") erme_solver::NonlinearResidual "C++ includes:
NonlinearResidual.hh ";

%feature("docstring")
erme_solver::NonlinearResidual::NonlinearResidual "erme_solver::NonlinearResidual::NonlinearResidual(SP_R R, SP_M M, SP_F
F, SP_A A, SP_L L)

Constructor.

Parameters:
-----------

R:  Pointer to response matrix

M:  Pointer to connectivity matrix

F:  Pointer to fission operator

A:  Pointer to absorption operator

L:  Pointer to leakage operator ";

%feature("docstring")  erme_solver::NonlinearResidual::compute_norm "double erme_solver::NonlinearResidual::compute_norm(Vector &x)

Computes the $ L_2 $ norm of the nonlinear residual.

The nonlinear residual is defined as \\\\[ \\\\mathbf{f(x)} = \\\\left
[\\\\begin{array}{c} (\\\\mathbf{M}\\\\mathbf{R}(k)-\\\\lambda
\\\\mathbf{I}) \\\\mathbf{J_-} \\\\\\\\
\\\\mathbf{F}(k)\\\\mathbf{J_-} - (k\\\\mathbf{L}(k)\\\\mathbf{J_-} )
\\\\\\\\ \\\\frac{1}{2} \\\\mathbf{J^T_-} \\\\mathbf{J_-} -
\\\\frac{1}{2} \\\\end{array} \\\\right ] = \\\\mathbf{0} \\\\, ,
\\\\] which is the same as used in the Newton-based schemes. The $ L_2
$ norm is then $ \\\\sqrt{ \\\\mathbf{f(x)}^T \\\\mathbf{f(x)} } $.

Parameters:
-----------

x:  vector of boundary unknowns with $ k $ and $ \\\\lambda $ ";

%feature("docstring")  erme_solver::NonlinearResidual::compute_norm "double erme_solver::NonlinearResidual::compute_norm(Vector &x, const
double k, const double l)

Computes the $ L_2 $ norm of the nonlinear residual.

This version offers an interface for Picard iteration. ";


// File: classerme__solver_1_1OperatorMR.xml
%feature("docstring") erme_solver::OperatorMR "

Performs the action of M*R.

C++ includes: OperatorMR.hh ";

%feature("docstring")  erme_solver::OperatorMR::OperatorMR "erme_solver::OperatorMR::OperatorMR(SP_R R, SP_M M)

Constructor.

Parameters:
-----------

R:  Response matrix

M:  Connectivity matrix ";


// File: classPCAppxJacobian.xml
%feature("docstring") PCAppxJacobian "

An approximate Jacobian for preconditioning.

For effective JFNK, a preconditioner is needed for the linear solves.
A straightforward and effective approach is to use an approximate
Jacobian. The approach here is to construct a partial Jacobian from
the initial guess, and to use it (or its factorization) for the rest
of the problem.

The approximate Jacobian is defined as \\\\[ \\\\mathbf{f'(x)} =
\\\\left [\\\\begin{array}{ccc} (\\\\mathbf{M}\\\\mathbf{R}-\\\\lambda
\\\\mathbf{I}) & 0 & \\\\mathbf{J_-} \\\\\\\\
(\\\\mathbf{F}-k\\\\mathbf{L}) & -\\\\mathbf{L} \\\\mathbf{J_-} & 0
\\\\\\\\ \\\\mathbf{J^T_-} & 0 & 0 \\\\end{array} \\\\right ] \\\\, ,
\\\\label{eq:jacobian} \\\\] where $ \\\\mathbf{J_-} $, $ k $, and $
\\\\lambda $ are the initial guess, likely found from a single crude
power iteration.

Notice this approximate Jacobian is almost complete, lacking only the
finite differences for the derivatives with respect to $ k $. No
studies have been performed to assess how much dropping these terms
affects the convergence of linear solves (or how much time is saved by
dropping the terms). Studies have shown this matrix is an effective
preconditioner when coupled with incomplete factorization.

C++ includes: PCAppxJacobian.hh ";

%feature("docstring")  PCAppxJacobian::PCAppxJacobian "PCAppxJacobian::PCAppxJacobian(integer a, integer b, integer c,
integer const nnz[], SermentVector *unk, GlobalProblem *pr)

Constructs an approximate Jacobian for preconditioning.

Parameters:
-----------

a:

b:

c:

nzz:

unk:

pr:  ";

%feature("docstring")  PCAppxJacobian::~PCAppxJacobian "PCAppxJacobian::~PCAppxJacobian()

Default destructor. ";

%feature("docstring")  PCAppxJacobian::myMatVec "void
PCAppxJacobian::myMatVec(Mat &M, Vec &f, Vec &fpf) ";


// File: classPowerIter.xml
%feature("docstring") PowerIter "

This class solves a GlobalProblem via power iteration.

The power (iteration) method is a standard procedure for finding the
largest eigenvalue of an operator. Here, the method is probably more
accurately called a Picard iteration, since the operator is strictly
nonlinear (whereas power iteration implies a linear operator).

Acceleration via Steffensen's method is available. In this case,
Aitken's $\\\\delta^2$ process is used to accelerate $k$ after two
outer iterations. The Aitken-extrapolants are fed back into response
function evaluation, and, in the best case, convergence is second
order.

The convergence criterion is currently limited to a nonlinear residual
evaluation that makes comparison to Newton methods easier, though
other criteria (perhaps on $ k $ alone) could be used.

The outer PowerIter is templated with an InnerIter. This can be the
built-in power iteration (as meant in the traditional sense) or a
SLEPc wrapper class allowing power and several Krylov iteration
schemes.

C++ includes: PowerIter.hh ";

%feature("docstring")  PowerIter::PowerIter "PowerIter< Inner
>::PowerIter(SP_globalproblem problem, SP_globalinput input)

Constructor. ";

%feature("docstring")  PowerIter::~PowerIter "PowerIter< Inner
>::~PowerIter() ";

%feature("docstring")  PowerIter::solve "void PowerIter< Inner
>::solve()

Solve the ERME via power iteration.

The inner eigenvalue problem is solved \\\\[ \\\\Big (
\\\\mathbf{M}\\\\mathbf{R}(k^{(m)}) - \\\\lambda \\\\mathbf{I}
J_{-}^{(n,m)} \\\\Big ) = 0 \\\\] within an outer iteration $ m $. The
inner iterations can be solved via the power method or one of SLEPc's
solvers.

The outer iteration is then defined by the eigenvalue update \\\\[
k^{(m+1)} = \\\\frac{ \\\\mathbf{F}(k^{(m)})J_{-}^{(n,m)} } {
\\\\mathbf{L}(k^{(m)})J_{-}^{(n,m)} } \\\\] for fission operator $
\\\\mathbf{F} $ and loss (absorption plus leakage) operator $
\\\\mathbf{L} $, and where $ n $ is the number of inner iterations
required for convergence on $ J_{-} $. ";


// File: classerme__response_1_1ResponseDatabase.xml
%feature("docstring") erme_response::ResponseDatabase "

Provides precomputed responses stored on disk.

Provides precomputed responses from a database.

This first implementation uses a singleton pattern. The reason for
this is so that just one instance of the database is created for all
db-derived responses. This model may break down when several local
groups are used.

We take the following approach. At the construction of the response
server, the database is read in completely and then broadcasted to all
nodes. Interpolation and expansion happens using data in memory.

The routines for reading the data can alternatively be used to read
data in as needed. We'll leave that for later.

C++ includes: ResponseDatabase.hh ";

%feature("docstring")
erme_response::ResponseDatabase::~ResponseDatabase "erme_response::ResponseDatabase::~ResponseDatabase()

Destructor. ";

%feature("docstring")  erme_response::ResponseDatabase::get "void
erme_response::ResponseDatabase::get(std::string nodename, SP_response
response, ResponseIndex index, const double keff)

Get the response for a node, index, and keff.

Parameters:
-----------

nodename:  name of node to look up in database

response:  response container to fill

index:  index of incident response

keff:  eigenvalue for requested response ";

%feature("docstring")  erme_response::ResponseDatabase::filename "std::string erme_response::ResponseDatabase::filename() const ";


// File: classResponseFunction.xml
%feature("docstring") ResponseFunction "

This base class is the foundation for holding and giving rf data.

Each instance of ResponseFunction (or more specifically, it's
subclasses) contains all the response function data for one element.
To pass RF data of all elements, a pointer array of ResponseFunction
objects is passed. Since we're using only pointers, this is pretty
memory-efficient, I think.

C++ includes: ResponseFunction.hh ";

%feature("docstring")  ResponseFunction::ResponseFunction "ResponseFunction::ResponseFunction(scalar *C, scalar *L, scalar *F,
scalar *A, scalar k) ";

%feature("docstring")  ResponseFunction::~ResponseFunction "ResponseFunction::~ResponseFunction() ";

%feature("docstring")  ResponseFunction::getCurrentResponse "scalar *
ResponseFunction::getCurrentResponse()

This function returns a pointer to current response data. ";

%feature("docstring")  ResponseFunction::getLeakageResponse "scalar *
ResponseFunction::getLeakageResponse()

This function returns a pointer to leakage response data. ";

%feature("docstring")  ResponseFunction::getFissionResponse "scalar *
ResponseFunction::getFissionResponse()

This function returns a pointer to fission response data. ";

%feature("docstring")  ResponseFunction::getAbsorptionResponse "scalar * ResponseFunction::getAbsorptionResponse()

This function returns a pointer to absorption response data. ";


// File: classResponseFunctionDiffusion.xml
%feature("docstring") ResponseFunctionDiffusion "

This base class is the foundation for holding and giving rf data.

Each instance of ResponseFunctionDiffusion contains all the response
function data for one element in a diffusion-based local problem. To
pass RF data of all elements, a pointer array of
ResponseFunctionDiffusion objects is passed. Since we're using only
pointers, this is pretty memory-efficient, I think.

C++ includes: ResponseFunctionDiffusion.hh ";

%feature("docstring")
ResponseFunctionDiffusion::ResponseFunctionDiffusion "ResponseFunctionDiffusion::ResponseFunctionDiffusion(scalar *C, scalar
*L, scalar *F, scalar *A, scalar k) ";

%feature("docstring")
ResponseFunctionDiffusion::~ResponseFunctionDiffusion "ResponseFunctionDiffusion::~ResponseFunctionDiffusion() ";


// File: classResponseFunctionServer.xml
%feature("docstring") ResponseFunctionServer "

To be completed.

C++ includes: ResponseFunctionServer.hh ";

%feature("docstring")  ResponseFunctionServer::ResponseFunctionServer
"ResponseFunctionServer::ResponseFunctionServer(GlobalInput &input)
";

%feature("docstring")  ResponseFunctionServer::~ResponseFunctionServer
"ResponseFunctionServer::~ResponseFunctionServer() ";

%feature("docstring")  ResponseFunctionServer::updateResponseFunctions
"ResponseFunction **
ResponseFunctionServer::updateResponseFunctions(scalar k)

This function updates the response functions via the server. ";

%feature("docstring")  ResponseFunctionServer::serverTime "scalar
ResponseFunctionServer::serverTime() ";

%feature("docstring")  ResponseFunctionServer::resetTime "scalar
ResponseFunctionServer::resetTime() ";


// File: structerme__response_1_1ResponseIndex.xml
%feature("docstring") erme_response::ResponseIndex "

Convenience container for indices.

Parameters:
-----------

n:  Unique (global) node index

s:  Node surface index

e:  Energy moment order

p:  Polar moment order

a:  Azimuthal moment order

s0:  First spatial dimension order

s1:  Second spatial dimension order

eo:  Is the combined function even=false or odd=true?

ni:  Moment index within the local node

C++ includes: ResponseIndex.hh ";

%feature("docstring")  erme_response::ResponseIndex::ResponseIndex "erme_response::ResponseIndex::ResponseIndex(size_t n=0, size_t s=0,
size_t e=0, size_t p=0, size_t a=0, size_t s0=0, size_t s1=0, bool
eo=false, size_t ni=0)

Constructor. ";


// File: classerme__response_1_1ResponseIndexer.xml
%feature("docstring") erme_response::ResponseIndexer "

Indexes a node response vector.

Each unique Node has assigned maximum orders in each phase space
variable on each surface of the node. Each variable is expanded in a
set of basis functions, and so the combined response is a tensor
product of such functions. Hence, the total order of a particular
response might be larger than the maximum allowed for any particular
variable. The user can limit the order of cross terms by setting the
order reduction level. The current options are: 0 - no reduction
(default) 1 - spatial order is limited by the largest of the two
spatial orders (applicable to 3D only) 2 - angular order is limited by
the largest of azimuthal and polar order (not applicable to 1D) 3 -
combination of 1 and 2 4 - maximum space-angle order is the maximum of
the available space and angle orders (not applicable to 1D)

Note, the indexer is constructed by all processes and applies to the
entire node list. However, a number of index functions help provide
indices into the global, local, and nodal moments.

The indexer is responsible for giving easy access into a moments
vector or operator row/column. This can be done via a global, local,
or nodal view, and within those, a unique view, since responses in
general can be repeated.

Relevant database entries: dimension

erme_order_reduction

C++ includes: ResponseIndexer.hh ";

%feature("docstring")  erme_response::ResponseIndexer::ResponseIndexer
"home robertsj Research serment source src erme_response
ResponseIndexer cc
erme_response::ResponseIndexer::ResponseIndexer(SP_db db, SP_nodelist
nodes)

Constructor.

Parameters:
-----------

db:  Pointer to parameter database

nodes:  Pointer to node list ";

%feature("docstring")  erme_response::ResponseIndexer::number_nodes "ResponseIndexer::size_t erme_response::ResponseIndexer::number_nodes()
const

Total number of nodes in the problem. ";

%feature("docstring")
erme_response::ResponseIndexer::number_node_moments "ResponseIndexer::size_t
erme_response::ResponseIndexer::number_node_moments(const size_t
node_ug) const

Number of moments associated with a node.

Parameters:
-----------

node_ug:  Unique global node index ";

%feature("docstring")
erme_response::ResponseIndexer::number_surface_moments "ResponseIndexer::size_t
erme_response::ResponseIndexer::number_surface_moments(const size_t
node_ug, const size_t surface_n) const

Number of moments on a node surface.

Parameters:
-----------

node_g:  Unique global node index

surface:  Node surface index ";

%feature("docstring")
erme_response::ResponseIndexer::number_unique_moments "ResponseIndexer::size_t
erme_response::ResponseIndexer::number_unique_moments() const

Return the number of unique moments of local nodes. ";

%feature("docstring")
erme_response::ResponseIndexer::number_local_moments "ResponseIndexer::size_t
erme_response::ResponseIndexer::number_local_moments() const

Return the number of moments of all local nodes. ";

%feature("docstring")
erme_response::ResponseIndexer::number_global_moments "ResponseIndexer::size_t
erme_response::ResponseIndexer::number_global_moments() const

Return the number of moments of all nodes. ";

%feature("docstring")  erme_response::ResponseIndexer::response_index
"ResponseIndex erme_response::ResponseIndexer::response_index(const
size_t node_g, const size_t surface_n, const size_t index_s) const

Get moment indices from cardinal index within node.

Parameters:
-----------

node_g:  Global index of node

surface_n:  Surface index of node

index_s:  Moment index on surface of node ";

%feature("docstring")
erme_response::ResponseIndexer::response_index_from_unique_local "ResponseIndex
erme_response::ResponseIndexer::response_index_from_unique_local(const
size_t index_ul) const

Get moment indices from unique cardinal index.

Parameters:
-----------

index_ul:  Moment index within unique local moments ";

%feature("docstring")
erme_response::ResponseIndexer::response_index_from_local "ResponseIndex
erme_response::ResponseIndexer::response_index_from_local(const size_t
index_l) const

Get moment indices from local cardinal index.

Parameters:
-----------

index_ul:  Moment index within local moments ";

%feature("docstring")
erme_response::ResponseIndexer::nodal_index_to_local "ResponseIndexer::size_t
erme_response::ResponseIndexer::nodal_index_to_local(const size_t
node_l, const size_t index_n) const

Get local moment index from a cardinal index within node.

Parameters:
-----------

node_l:  Local node index

index_n:  Cardinal moment index within node ";

%feature("docstring")
erme_response::ResponseIndexer::global_index_to_local "int
erme_response::ResponseIndexer::global_index_to_local(const size_t
index_g) const

Get local moment index from global moment index.

Note, this returns -1 if the global index doesn't represent a local
index.

Parameters:
-----------

index_g:  Global node index ";

%feature("docstring")
erme_response::ResponseIndexer::nodal_index_to_global "ResponseIndexer::size_t
erme_response::ResponseIndexer::nodal_index_to_global(const size_t
node_g, const size_t index_n) const

Get global moment index from a nodal moment index.

Parameters:
-----------

node_g:  Global node index

index_n:  Moment index within node ";

%feature("docstring")
erme_response::ResponseIndexer::local_index_to_global "ResponseIndexer::size_t
erme_response::ResponseIndexer::local_index_to_global(const size_t
index_l) const

Get global moment index from a local moment index.

Parameters:
-----------

index_l:  Local moment index ";

%feature("docstring")
erme_response::ResponseIndexer::local_index_to_unique "ResponseIndexer::size_t
erme_response::ResponseIndexer::local_index_to_unique(const size_t
index_l) const

Get the unique local index from cardinal local index.

Parameters:
-----------

index_l:  Local moment index ";

%feature("docstring")  erme_response::ResponseIndexer::display "void
erme_response::ResponseIndexer::display() const

Display the indices in a nice format. ";


// File: classResponseMatrix.xml
%feature("docstring") ResponseMatrix "

This base class is part of the foundation of a whole response matrix.

This class contains items specific to response matrices. However,
since we wish to leverage the various matrix classes (full versus
shell), we need multiple inheritance.

C++ includes: ResponseMatrix.hh ";

%feature("docstring")  ResponseMatrix::ResponseMatrix "ResponseMatrix::ResponseMatrix(GlobalInput &input,
ResponseFunctionServer *s) ";

%feature("docstring")  ResponseMatrix::~ResponseMatrix "ResponseMatrix::~ResponseMatrix() ";


// File: classerme_1_1ResponseMatrix.xml
%feature("docstring") erme::ResponseMatrix "

Response matrix operator.

C++ includes: ResponseMatrix.hh ";

%feature("docstring")  erme::ResponseMatrix::ResponseMatrix "home
robertsj Research serment source src erme ResponseMatrix cc
ResponseMatrix::ResponseMatrix(SP_nodelist nodes, SP_indexer indexer,
SP_server server)

Constructor.

Parameters:
-----------

nodes:  Pointer to node list

indexer:  Pointer to response indexer

server:  Pointer to response server ";

%feature("docstring")  erme::ResponseMatrix::update "void
ResponseMatrix::update()

Update the response matrix data. ";


// File: classResponseMatrixFull.xml
%feature("docstring") ResponseMatrixFull "

This base class is part of the foundation of a whole response matrix.

This class contains items specific to response matrices. However,
since we wish to leverage the various matrix classes (full versus
shell), we need multiple inheritance.

C++ includes: ResponseMatrixFull.hh ";

%feature("docstring")  ResponseMatrixFull::ResponseMatrixFull "ResponseMatrixFull::ResponseMatrixFull(GlobalInput &input,
ResponseFunctionServer *s) ";

%feature("docstring")  ResponseMatrixFull::~ResponseMatrixFull "ResponseMatrixFull::~ResponseMatrixFull() ";

%feature("docstring")  ResponseMatrixFull::updateData "void
ResponseMatrixFull::updateData(scalar k)

This function updates underlying matrix elements. ";


// File: classResponseOperator.xml
%feature("docstring") ResponseOperator "

This base class has general traits of a response operator.

Such operators include the response matrix, leakage response matric,
etc.

C++ includes: ResponseOperator.hh ";

%feature("docstring")  ResponseOperator::ResponseOperator "ResponseOperator::ResponseOperator(GlobalInput &input,
ResponseFunctionServer *s) ";

%feature("docstring")  ResponseOperator::~ResponseOperator "ResponseOperator::~ResponseOperator() ";

%feature("docstring")  ResponseOperator::updateData "virtual void
ResponseOperator::updateData(scalar k)=0 ";


// File: classerme_1_1ResponseOperator.xml
%feature("docstring") erme::ResponseOperator "

Base class for response operators.

C++ includes: ResponseOperator.hh ";

%feature("docstring")  erme::ResponseOperator::ResponseOperator "erme::ResponseOperator::ResponseOperator(SP_nodelist nodes, SP_indexer
indexer, SP_server server)

Constructor.

Parameters:
-----------

indexer:  Pointer to response indexer

server:  Pointer to response server ";

%feature("docstring")  erme::ResponseOperator::~ResponseOperator "virtual erme::ResponseOperator::~ResponseOperator()

Virtual Destructor. ";

%feature("docstring")  erme::ResponseOperator::update "virtual void
erme::ResponseOperator::update()=0

Update responses.

This assumes that the response server is updated. ";


// File: classerme__response_1_1ResponseServer.xml
%feature("docstring") erme_response::ResponseServer "

Serve nodal responses to clients.

A ResponseServer lives on a local communicator. A server is in charge
of one or more Node objects. The nodal responses are produced by a
ResponseSource that solves the local problems. There is one source for
each unique node.

C++ includes: ResponseServer.hh ";

%feature("docstring")  erme_response::ResponseServer::ResponseServer "home robertsj Research serment source src erme_response ResponseServer
cc home robertsj Research serment source src erme_response
ResponseServer cc
erme_response::ResponseServer::ResponseServer(SP_nodelist nodes,
SP_indexer indexer, std::string dbname=\"\", size_t dborder=1)

Constructor.

Parameters:
-----------

nodes:  Pointer to node list

indexer:  Pointer to indexer

dbname:  Filename of response database (optional)

dborder:  Interpolation order for database (optional) ";

%feature("docstring")  erme_response::ResponseServer::update "void
erme_response::ResponseServer::update(const double keff)

Update the eigenvalue and compute the new responses. ";

%feature("docstring")  erme_response::ResponseServer::response "ResponseServer::SP_response
erme_response::ResponseServer::response(size_t node)

Return a nodal response.

Parameters:
-----------

node:  Local node index ";


// File: classerme__response_1_1ResponseSource.xml
%feature("docstring") erme_response::ResponseSource "

Abstract response source.

A ResponseSource provides its ResponseServer with responses for use in
the global solve. Each ResponseSource is unique for a given Node. The
ResponseSource represents the interface between Serment and local
solvers such as Detran. Each concrete Node implementation must have a
corresponding concrete ResponseSource and ResponseSourceFactory::build
specialization.

C++ includes: ResponseSource.hh ";

%feature("docstring")  erme_response::ResponseSource::ResponseSource "erme_response::ResponseSource::ResponseSource(SP_node node)

Constructor.

Parameters:
-----------

node:  Pointer to node object for which this source generates
responses ";

%feature("docstring")  erme_response::ResponseSource::~ResponseSource
"virtual erme_response::ResponseSource::~ResponseSource()

Virtual destructor. ";

%feature("docstring")  erme_response::ResponseSource::update "void
erme_response::ResponseSource::update(const double keff)

Update the k-eigenvalue. ";

%feature("docstring")  erme_response::ResponseSource::compute "virtual void erme_response::ResponseSource::compute(SP_response
response, ResponseIndex index)=0

Compute a response for the requested incident index.

The client passes the response to be updated and the index of the
corresponding incident condition. The client is then responsible for
moving data from the sources to the server.

Parameters:
-----------

response:  Pointer to response object to be updated

index:  Response indices ";


// File: classerme__response_1_1ResponseSourceDatabase.xml
%feature("docstring") erme_response::ResponseSourceDatabase "C++
includes: ResponseSourceDatabase.hh ";

%feature("docstring")
erme_response::ResponseSourceDatabase::ResponseSourceDatabase "erme_response::ResponseSourceDatabase::ResponseSourceDatabase(SP_node
node)

Constructor.

Parameters:
-----------

node:  Pointer to a response database node ";

%feature("docstring")
erme_response::ResponseSourceDatabase::~ResponseSourceDatabase "erme_response::ResponseSourceDatabase::~ResponseSourceDatabase()

Virtual destructor. ";

%feature("docstring")  erme_response::ResponseSourceDatabase::compute
"void erme_response::ResponseSourceDatabase::compute(SP_response
response, ResponseIndex index)

Compute a response for the requested incident index. ";


// File: classerme__response_1_1ResponseSourceDetran.xml
%feature("docstring") erme_response::ResponseSourceDetran "

Compute responses using Detran.

C++ includes: ResponseSourceDetran.hh ";

%feature("docstring")
erme_response::ResponseSourceDetran::ResponseSourceDetran "erme_response::ResponseSourceDetran::ResponseSourceDetran(SP_node
node)

Constructor. ";

%feature("docstring")  erme_response::ResponseSourceDetran::compute "void erme_response::ResponseSourceDetran::compute(SP_response
response, ResponseIndex index)

Compute a response for the requested incident index.

The client passes the response to be updated and the index of the
corresponding incident condition. The client is then responsible for
moving data from the sources to the server.

Parameters:
-----------

response:  Pointer to response object to be updated

index:  Response indices ";


// File: classerme__response_1_1ResponseSourceDummy.xml
%feature("docstring") erme_response::ResponseSourceDummy "

Fake response source for testing purposes.

C++ includes: ResponseSourceDummy.hh ";

%feature("docstring")
erme_response::ResponseSourceDummy::ResponseSourceDummy "erme_response::ResponseSourceDummy::ResponseSourceDummy(SP_node node)

Constructor. ";

%feature("docstring")  erme_response::ResponseSourceDummy::compute "void erme_response::ResponseSourceDummy::compute(SP_response response,
ResponseIndex index)

Compute a response for the requested incident index.

The client passes the response to be updated and the index of the
corresponding incident condition. The client is then responsible for
moving data from the sources to the server.

Parameters:
-----------

response:  Pointer to response object to be updated

index:  Response indices ";


// File: classerme__response_1_1ResponseSourceFactory.xml
%feature("docstring") erme_response::ResponseSourceFactory "

Constructs response.

C++ includes: ResponseSourceFactory.hh ";

%feature("docstring")  erme_response::ResponseSourceFactory::build "SP_source erme_response::ResponseSourceFactory::build(SP_NODE node)

Build a response source.

This is a factory method that calls an implemention for each type of
node to be constructed.

Parameters:
-----------

node:  Smart pointer to node ";


// File: classSermentMatrix.xml
%feature("docstring") SermentMatrix "

An abstract matrix class that encapsulates an external matrix class.

SermentMatrix is the base class used for all matrix quantities used in
SERMENT's linear operations. The base (i.e. abstract) class and its
subclasses support one common contract, namely to provide action of a
vector. Since the ultimate implementation depends on the underlying
matrix format, this is defined as a virtual method with appropriate
signatures.

C++ includes: SermentMatrix.hh ";

%feature("docstring")  SermentMatrix::SermentMatrix "SermentMatrix::SermentMatrix(integer a, integer b) ";

%feature("docstring")  SermentMatrix::~SermentMatrix "SermentMatrix::~SermentMatrix() ";

%feature("docstring")  SermentMatrix::matVec "virtual void
SermentMatrix::matVec(SermentVector &x, SermentVector &y)=0 ";

%feature("docstring")  SermentMatrix::matVec "virtual void
SermentMatrix::matVec(SermentVector x)=0 ";

%feature("docstring")  SermentMatrix::releaseMe "void
SermentMatrix::releaseMe()

Deallocate the memory for the Petsc matrix explicitly. ";


// File: classSermentMatrixBCRS.xml
%feature("docstring") SermentMatrixBCRS "

This concrete matrix class encapsulates a Petsc block CSR matrix.

SermentMatrixCRS builds off SermentMatrixFull. It implements methods
specific to a complete matrix construction of a compressed row storage
matrix.

C++ includes: SermentMatrixBCRS.hh ";

%feature("docstring")  SermentMatrixBCRS::SermentMatrixBCRS "SermentMatrixBCRS::SermentMatrixBCRS(integer a, integer b, integer c,
integer d) ";

%feature("docstring")  SermentMatrixBCRS::~SermentMatrixBCRS "SermentMatrixBCRS::~SermentMatrixBCRS() ";

%feature("docstring")  SermentMatrixBCRS::insertVals "void
SermentMatrixBCRS::insertVals(scalar values[], integer numrow, integer
idxrow[], integer numcol, integer idxcol[])

This function insert values into a block CSR matrix.

This function adds values to the matrix given a 2-D array of values,
the corresponding row and column counts, and indices for placing those
values. Note, the values are placed in block format, but they
themselves are still in a 1d array. For example, if m=n=2, and we want
the following blocks 1 2 | 3 4 5 6 | 7 8 - - | - - - 9 10 | 11 12 13
14 | 15 16 the vector v would have {1, 5, ..., 16} sequentially, i.e.
in column-major form. Likely, the BCSR matrices to be computed in
Serment, namely the response matrix, will be constructed one block at
a time, since the current response blocks are stored in separate 1d
arrays for each element. ";


// File: classSermentMatrixCRS.xml
%feature("docstring") SermentMatrixCRS "

This concrete matrix class encapsulates a Petsc CSR matrix.

SermentMatrixCRS builds off SermentMatrixFull. It implements methods
specific to a complete matrix construction of a compressed row storage
matrix.

C++ includes: SermentMatrixCRS.hh ";

%feature("docstring")  SermentMatrixCRS::SermentMatrixCRS "SermentMatrixCRS::SermentMatrixCRS(integer a, integer b, integer c) ";

%feature("docstring")  SermentMatrixCRS::SermentMatrixCRS "SermentMatrixCRS::SermentMatrixCRS(integer a, integer b, integer c,
const integer nnz[]) ";

%feature("docstring")  SermentMatrixCRS::~SermentMatrixCRS "SermentMatrixCRS::~SermentMatrixCRS() ";

%feature("docstring")  SermentMatrixCRS::insertVals "void
SermentMatrixCRS::insertVals(scalar values[], integer numrow, integer
idxrow[], integer numcol, integer idxcol[])

Insert values.

This function adds values to the matrix given a 2-D array of values,
the corresponding row and column counts, and indices for placing those
values. ";

%feature("docstring")  SermentMatrixCRS::insertVal "void
SermentMatrixCRS::insertVal(scalar value, integer row, integer col)

Insert value.

This function adds a single values to the matrix given an index. ";


// File: classSermentMatrixFull.xml
%feature("docstring") SermentMatrixFull "

An abstract matrix class encapsulating a full (constructed) matrix.

SermentMatrixFull is an abstract class that builds off SermentMatrix.
It specifies additional methods specific to a complete matrix
construction, e.g. the case of a sequential compressed row storage
matrix.

C++ includes: SermentMatrixFull.hh ";

%feature("docstring")  SermentMatrixFull::SermentMatrixFull "SermentMatrixFull::SermentMatrixFull(integer a, integer b, integer c)
";

%feature("docstring")  SermentMatrixFull::~SermentMatrixFull "SermentMatrixFull::~SermentMatrixFull() ";

%feature("docstring")  SermentMatrixFull::matVec "void
SermentMatrixFull::matVec(SermentVector &x, SermentVector &y)

Matrix-vector multiplication, Mx-->y.

This function performs matrix-vector multiplcication of the form y=Mx.
The actual multiplication is done by the Petsc MatMult function, where
x and y cannot be the same. ";

%feature("docstring")  SermentMatrixFull::matVec "void
SermentMatrixFull::matVec(SermentVector x)

Matrix-vector multiplication, Mx-->x.

This function performs matrix-vector multiplcication of the form x=Mx.
The actual multiplication is done by the Petsc MatMult function.
Because Petsc's function does not allow x and y to be the same, we
initiate a local Petsc with vector into which the results are placed.
The contents of x.V are then replaced with the output. ";

%feature("docstring")  SermentMatrixFull::insertVals "virtual void
SermentMatrixFull::insertVals(scalar values[], integer numrow, integer
idxrow[], integer numcol, integer idxcol[])=0 ";

%feature("docstring")  SermentMatrixFull::checkReady "void
SermentMatrixFull::checkReady()

Method check whether SermentVector is ready for operations.

Because the underlying external library (here Petsc) has special
operations to ready vectors for operations, they are inserted here as
a check to be called before any operation. If isReady is true, nothing
happens, but if it is false, the vector is assembled. ";

%feature("docstring")  SermentMatrixFull::viewMe "void
SermentMatrixFull::viewMe()

This function displays the matrix contents. ";


// File: classSermentMatrixShell.xml
%feature("docstring") SermentMatrixShell "

This is an abstract class encapsulating a shell matrix (action only).

SermentMatrixShell is an abstract class that builds off SermentMatrix.
It specifies additional methods specific for the internal mechanism of
matrix-vector operations

C++ includes: SermentMatrixShell.hh ";

%feature("docstring")  SermentMatrixShell::SermentMatrixShell "SermentMatrixShell::SermentMatrixShell(integer a, integer b, void
*ctx) ";

%feature("docstring")  SermentMatrixShell::~SermentMatrixShell "SermentMatrixShell::~SermentMatrixShell() ";

%feature("docstring")  SermentMatrixShell::matVec "void
SermentMatrixShell::matVec(SermentVector &x, SermentVector &y)

Matrix-vector multiplication, Mx-->y.

This function performs matrix-vector multiplcication of the form y=Mx.
The actual multiplication is done by the MyMatVec function, allowing
for a completely matrix-free approach. ";

%feature("docstring")  SermentMatrixShell::matVec "void
SermentMatrixShell::matVec(SermentVector x) ";

%feature("docstring")  SermentMatrixShell::myMatVec "virtual void
SermentMatrixShell::myMatVec(Mat &A, Vec &x, Vec &y)=0 ";


// File: classSermentVector.xml
%feature("docstring") SermentVector "

A Vector class that encapsulates an external (PETSc) vector class.

SermentVector is a class used for all vector quantities that are
operated on by matrices (i.e. SermentMatrices). SermentVector
encapsulates an external package's vector to hide external-specific
code within the main body of Serment. Currently, PETSc is used for all
linear algebra. Also, only a sequential code is being implemented
currently. However, because all vector (and matrix) construction is
hidded in these encapsulations, a move to parallel will not change the
main code.

C++ includes: SermentVector.hh ";

%feature("docstring")  SermentVector::SermentVector "SermentVector::SermentVector(integer m) ";

%feature("docstring")  SermentVector::~SermentVector "SermentVector::~SermentVector() ";

%feature("docstring")  SermentVector::insertVal "void
SermentVector::insertVal(integer row, scalar value)

Method to insert single indexed value into the vector.

Note, this is not the ideal way to construct Petsc vectors and is
included primarily for testing. ";

%feature("docstring")  SermentVector::insertVals "void
SermentVector::insertVals(integer ni, const integer ix[], const scalar
y[])

Method to insert several indexed values into the matrix.

For standard compressed row storage matrices, this is the preferred
method for inserting values. Arguments: ni - number of elements to add
ix - indices where to add y - values to add The Petsc documentation is
as follows: x - vector to insert in ni - number of elements to add ix
- indices where to add y - array of values iora - either INSERT_VALUES
or ADD_VALUES, where ADD_VALUES adds values to any existing entries,
and INSERT_VALUES replaces existing entries with new values ";

%feature("docstring")  SermentVector::vecSet "void
SermentVector::vecSet(scalar a)

Set all entries to single value. ";

%feature("docstring")  SermentVector::vecAVPY "void
SermentVector::vecAVPY(scalar a, SermentVector &Y)

V = a*V + Y. ";

%feature("docstring")  SermentVector::vecAYPV "void
SermentVector::vecAYPV(scalar a, SermentVector &Y)

V = a*Y + V. ";

%feature("docstring")  SermentVector::vecDot "scalar
SermentVector::vecDot(SermentVector &Y)

val = y'x ";

%feature("docstring")  SermentVector::vecPointMult "void
SermentVector::vecPointMult(SermentVector &Y, SermentVector &Z)

Z = X.*Y. ";

%feature("docstring")  SermentVector::vecScale "void
SermentVector::vecScale(scalar a)

x = a*x

Parameters:
-----------

a:  scalar value by which the vector x is scaled ";

%feature("docstring")  SermentVector::vecCopy "void
SermentVector::vecCopy(SermentVector &Y)

Copy me to Y.

Parameters:
-----------

a:  scalar value by which the vector x is scaled ";

%feature("docstring")  SermentVector::Length "integer
SermentVector::Length()

Return my length. ";

%feature("docstring")  SermentVector::releaseMe "void
SermentVector::releaseMe()

Method to release memory for the Petsc Vec.

Ideally, this method would be part of the destructor, but since it
would be called after PetscFinalize---which seems to do deallocations
of its own---it's best to use explicit deallocation when a
SermentVector is no longer needed. ";

%feature("docstring")  SermentVector::checkReady "void
SermentVector::checkReady()

Method check whether SermentVector is ready for operations.

Because the underlying external library (here Petsc) has special
operations to ready vectors for operations, they are inserted here as
a check to be called before any operation. If isReady is true, nothing
happens, but if it is false, the vector is assembled. ";

%feature("docstring")  SermentVector::viewMe "void
SermentVector::viewMe()

View the vector contents. ";


// File: classState.xml
%feature("docstring") State "

Represents the problem state vector.

A solution for the eigenvalue response matrix equations consists of a
global boundary vector, which contains moments for each surface of
each cell, and the k-eigenvalue.

C++ includes: StateERME.hh ";


// File: classerme_1_1StateERME.xml
%feature("docstring") erme::StateERME "C++ includes: StateERME.hh ";

%feature("docstring")  erme::StateERME::StateERME "home robertsj
Research serment source src erme StateERME cc
erme::StateERME::StateERME(const size_t size)

Constructor.

Parameters:
-----------

size:  Size of local state vector ";

%feature("docstring")  erme::StateERME::set_k "void
erme::StateERME::set_k(const double k_val) ";

%feature("docstring")  erme::StateERME::set_lambda "void
erme::StateERME::set_lambda(const double lambda_val) ";

%feature("docstring")  erme::StateERME::k "double
erme::StateERME::k() const ";

%feature("docstring")  erme::StateERME::lambda "double
erme::StateERME::lambda() const ";

%feature("docstring")  erme::StateERME::local_size "StateERME::size_t
erme::StateERME::local_size() const ";

%feature("docstring")  erme::StateERME::global_size "StateERME::size_t erme::StateERME::global_size() const ";

%feature("docstring")  erme::StateERME::moments "const
StateERME::Vector & erme::StateERME::moments() const

Const reference to moments vector. ";

%feature("docstring")  erme::StateERME::moments "StateERME::Vector &
erme::StateERME::moments()

Mutable reference to moments vector. ";


// File: classerme__solver_1_1SteffensenUpdate.xml
%feature("docstring") erme_solver::SteffensenUpdate "C++ includes:
SteffensenUpdate.hh ";

%feature("docstring")  erme_solver::SteffensenUpdate::SteffensenUpdate
"erme_solver::SteffensenUpdate::SteffensenUpdate()

Constructor. ";

%feature("docstring")
erme_solver::SteffensenUpdate::~SteffensenUpdate "virtual
erme_solver::SteffensenUpdate::~SteffensenUpdate()

Virtual destructor. ";

%feature("docstring")  erme_solver::SteffensenUpdate::compute "double
erme_solver::SteffensenUpdate::compute(const double keff, SP_vector J)

Computes and updated keff, possibly based on previous history.

Parameters:
-----------

keff:  latest eigenvalue estimate

J:  latest boundary unknowns ";


// File: classerme_1_1TwoGroupKernel.xml
%feature("docstring") erme::TwoGroupKernel "

Simple two-dimensional, two-group response function generator.

Simple diff2d alternative that borrows some of its code.

This class produces response functions for homogeneous cells based on
a semi-analytical model augmented with fitting parameters. For the
tests used to guide fitting, the response functions are all within
about 10% relative error, and in absolute terms, the residuals were
usually below 0.001.

This aims to be a simple, fast solver limited to two group,
homogeneous nodes without upscatter.

C++ includes: TwoGroupKernel.hh ";

%feature("docstring")  erme::TwoGroupKernel::TwoGroupKernel "home
robertsj Research serment source src erme local TwoGroupKernel cc home
robertsj Research serment source src erme local TwoGroupKernel cc
erme::TwoGroupKernel::TwoGroupKernel()

Constructor. ";

%feature("docstring")  erme::TwoGroupKernel::~TwoGroupKernel "erme::TwoGroupKernel::~TwoGroupKernel()

Destructor. ";

%feature("docstring")  erme::TwoGroupKernel::response "void
erme::TwoGroupKernel::response(SP_responsefunction responses, Vec_Dbl
&data)

Get the responses.

Responses are returned as a function of the following: + Delta --
Assembly dimension + k -- k-effective + D1 -- Group 1 diffusion
coefficient + D2 -- Group 2 diffusion coefficient + R1 -- Group 1
removal cross-section + A2 -- Group 2 absorption cross-section + F1 --
Group 1 fission cross-section times nu + F2 -- Group 2 fission cross-
section times nu + S12 -- Group 1 to 2 scattering cross-section These
data are held sequentially in the data vector.

Parameters:
-----------

responses:  Response function smart pointer to be filled.

data:  Two group data. ";


// File: classerme_1_1TwoGroupKernelSolver.xml
%feature("docstring") erme::TwoGroupKernelSolver "C++ includes:
TwoGroupKernelSolver.hh ";

%feature("docstring")
erme::TwoGroupKernelSolver::TwoGroupKernelSolver "home robertsj
Research serment source src erme local TwoGroupKernelSolver cc
erme::TwoGroupKernelSolver::TwoGroupKernelSolver(int nh, int order)

Constructor.

Initializes the matrix, solvers, etc.

Parameters:
-----------

nh:  number of spatial mesh

order:  order of responses to generate ";

%feature("docstring")  erme::TwoGroupKernelSolver::update_matrix "void erme::TwoGroupKernelSolver::update_matrix(Vec_Dbl &data)

Update the matrix.

Parameters:
-----------

data:  group data for this state ";

%feature("docstring")  erme::TwoGroupKernelSolver::set_source "void
erme::TwoGroupKernelSolver::set_source(Vec_Dbl &data, int group, int
order) ";

%feature("docstring")  erme::TwoGroupKernelSolver::solve "void
erme::TwoGroupKernelSolver::solve()

Solve the system. ";


// File: classlinear__algebra_1_1Vector.xml
%feature("docstring") linear_algebra::Vector "

Lightweight wrapper for PETSc Vec.

Use of Vector and the corresponding Matrix class should, in theory,
eliminate a lot of PETSc code from the rest of Serment.

C++ includes: Vector.hh ";

%feature("docstring")  linear_algebra::Vector::Vector "linear_algebra::Vector::Vector(const size_type m, const double
val=0.0)

Constructor.

Parameters:
-----------

m:  Local number of rows

val:  Optional initial value ";

%feature("docstring")  linear_algebra::Vector::Vector "linear_algebra::Vector::Vector(const Vector &V)

Copy Constructor.

Parameters:
-----------

V:   Vector to copy ";

%feature("docstring")  linear_algebra::Vector::~Vector "linear_algebra::Vector::~Vector()

Destructor. ";

%feature("docstring")  linear_algebra::Vector::insert_values "void
linear_algebra::Vector::insert_values(const unsigned int number, const
int *rows, const double *values)

Insert values.

Parameters:
-----------

values:  Array of values to insert

number:  Number of values to insert

rows:  Indices of rows where values are inserted ";

%feature("docstring")  linear_algebra::Vector::assemble "void
linear_algebra::Vector::assemble()

Assemble the vector. ";

%feature("docstring")  linear_algebra::Vector::norm "double
linear_algebra::Vector::norm(const int type=L2)

Compute my norm. ";

%feature("docstring")  linear_algebra::Vector::norm_residual "double
linear_algebra::Vector::norm_residual(const Vector &x, const int
type=L2) ";

%feature("docstring")  linear_algebra::Vector::dot "double
linear_algebra::Vector::dot(Vector &x) ";

%feature("docstring")  linear_algebra::Vector::scale "void
linear_algebra::Vector::scale(const double factor) ";

%feature("docstring")  linear_algebra::Vector::set "void
linear_algebra::Vector::set(const double v) ";

%feature("docstring")  linear_algebra::Vector::add "void
linear_algebra::Vector::add(const Vector &x) ";

%feature("docstring")  linear_algebra::Vector::subtract "void
linear_algebra::Vector::subtract(const Vector &x) ";

%feature("docstring")  linear_algebra::Vector::multiply "void
linear_algebra::Vector::multiply(const Vector &x) ";

%feature("docstring")  linear_algebra::Vector::divide "void
linear_algebra::Vector::divide(const Vector &x) ";

%feature("docstring")  linear_algebra::Vector::copy "void
linear_algebra::Vector::copy(const Vector &x) ";

%feature("docstring")  linear_algebra::Vector::add_a_times_x "void
linear_algebra::Vector::add_a_times_x(const double a, const Vector &x)
";

%feature("docstring")  linear_algebra::Vector::V "Vec
linear_algebra::Vector::V() const

Return the PETSc vector. ";

%feature("docstring")  linear_algebra::Vector::global_size "int
linear_algebra::Vector::global_size() const

Return the global size. ";

%feature("docstring")  linear_algebra::Vector::local_size "int
linear_algebra::Vector::local_size() const

Return the local size. ";

%feature("docstring")  linear_algebra::Vector::lower_bound "int
linear_algebra::Vector::lower_bound() const

Return the lower bound. ";

%feature("docstring")  linear_algebra::Vector::upper_bound "int
linear_algebra::Vector::upper_bound() const

Return the upper bound. ";

%feature("docstring")  linear_algebra::Vector::is_assembled "bool
linear_algebra::Vector::is_assembled() const

Return assembled flag. ";

%feature("docstring")  linear_algebra::Vector::display "void
linear_algebra::Vector::display() const

View via standard output. ";


// File: classutil_1_1Vector__Lite.xml
%feature("docstring") util::Vector_Lite "

Array container that is a wrapper around a standard C array.

It adds iterator and arithemtic support, along with bounds checking
(via Nemesis' DBC).

An alternative to this class is boost::array (www.boost.org). However,
boost::array is an aggregate type, which has advantages (can use
initializers) and disadvantages (public data, cannot be a base class).
boost::array also doesn't do bounds checking.

Parameters:
-----------

T:  Type of each array element.

N:  Length of array.

C++ includes: Vector_Lite.hh ";

%feature("docstring")  util::Vector_Lite::Vector_Lite "util::Vector_Lite< T, N >::Vector_Lite(const T &u=T())

Constructor based on a scalar value.

Initializes all values to u. Note that this ctor also acts as the
default ctor.

Parameters:
-----------

u:  Scalar value. ";

%feature("docstring")  util::Vector_Lite::Vector_Lite "util::Vector_Lite< T, N >::Vector_Lite(const T &u0, const T &u1)

Constructor for N = 2.

Parameters:
-----------

u0:  1st element.

u1:  2nd element. ";

%feature("docstring")  util::Vector_Lite::Vector_Lite "util::Vector_Lite< T, N >::Vector_Lite(const T &u0, const T &u1, const
T &u2)

Constructor for N = 3.

Parameters:
-----------

u0:  1st element.

u1:  2nd element.

u2:  3rd element. ";

%feature("docstring")  util::Vector_Lite::Vector_Lite "util::Vector_Lite< T, N >::Vector_Lite(const T &u0, const T &u1, const
T &u2, const T &u3)

Constructor for N = 4.

Parameters:
-----------

u0:  1st element.

u1:  2nd element.

u2:  3rd element.

u3:  4th element. ";

%feature("docstring")  util::Vector_Lite::Vector_Lite "util::Vector_Lite< T, N >::Vector_Lite(const T &u0, const T &u1, const
T &u2, const T &u3, const T &u4)

Constructor for N = 5.

Parameters:
-----------

u0:  1st element.

u1:  2nd element.

u2:  3rd element.

u3:  4th element.

u4:  5th element. ";

%feature("docstring")  util::Vector_Lite::fill "void
util::Vector_Lite< T, N >::fill(const T u[N])

Fill in from C array.

Parameters:
-----------

u:  c-style pointer to array of size N. ";

%feature("docstring")  util::Vector_Lite::~Vector_Lite "util::Vector_Lite< T, N >::~Vector_Lite(void)

Destructor. ";

%feature("docstring")  util::Vector_Lite::valid_index "bool
util::Vector_Lite< T, N >::valid_index(const size_type i) const

Returns true if i is a valid array index. ";

%feature("docstring")  util::Vector_Lite::begin "iterator
util::Vector_Lite< T, N >::begin()

Iterator begin. ";

%feature("docstring")  util::Vector_Lite::begin "const_iterator
util::Vector_Lite< T, N >::begin() const

Const iterator begin. ";

%feature("docstring")  util::Vector_Lite::end "iterator
util::Vector_Lite< T, N >::end()

Iterator end. ";

%feature("docstring")  util::Vector_Lite::end "const_iterator
util::Vector_Lite< T, N >::end() const

Const iterator end. ";

%feature("docstring")  util::Vector_Lite::size "size_type
util::Vector_Lite< T, N >::size() const

Number of elements ( N); for STL support. ";

%feature("docstring")  util::Vector_Lite::max_size "size_type
util::Vector_Lite< T, N >::max_size() const

Max number of elements ( N); for STL support. ";

%feature("docstring")  util::Vector_Lite::empty "bool
util::Vector_Lite< T, N >::empty() const

True if N = 0; for STL support. ";


// File: structxmlFileInfo.xml
%feature("docstring") xmlFileInfo "C++ includes: InputXML.hh ";


// File: namespaceconstants.xml


// File: namespacedetran.xml


// File: namespacedetran__test.xml


// File: namespaceerme.xml


// File: namespaceerme__geometry.xml
%feature("docstring")  erme_geometry::cartesian_node_detran "Node::SP_node erme_geometry::cartesian_node_detran(const int dim)

The Detran Cartesian nodes are all 10 cm in extent with a
monoenergetic material. All applicable responses are first order. ";

%feature("docstring")  erme_geometry::cartesian_node_detran_list_2d "NodeList::SP_nodelist erme_geometry::cartesian_node_detran_list_2d()
";

%feature("docstring")  erme_geometry::cartesian_node_dummy_list_1d "NodeList::SP_nodelist erme_geometry::cartesian_node_dummy_list_1d() ";

%feature("docstring")  erme_geometry::cartesian_node_dummy_list_2d "NodeList::SP_nodelist erme_geometry::cartesian_node_dummy_list_2d(int
so=4, int ao=2, int po=2) ";

%feature("docstring")  erme_geometry::cartesian_node_dummy_list_3d "NodeList::SP_nodelist erme_geometry::cartesian_node_dummy_list_3d() ";

%feature("docstring")
erme_geometry::cartesian_node_dummy_list_2d_variable "NodeList::SP_nodelist
erme_geometry::cartesian_node_dummy_list_2d_variable(int so, int N) ";

%feature("docstring")
erme_geometry::cartesian_node_dummy_list_3d_variable "NodeList::SP_nodelist
erme_geometry::cartesian_node_dummy_list_3d_variable(int so, int N) ";


// File: namespaceerme__response.xml
%feature("docstring")  erme_response::interpolate_linear "double
erme_response::interpolate_linear(double x, double x0, double x1,
double r0, double r1)

Interpolate a function evaluated at two points using linear
interpolation. ";

%feature("docstring")  erme_response::interpolate_quadratic "double
erme_response::interpolate_quadratic(double x, double x0, double x1,
double x2, double r0, double r1, double r2)

Interpolate a function evaluated at three points using quadratic
interpolation.

This was generated using Maple's code generator with optimization. ";

%feature("docstring")  erme_response::interpolate_cubic "double
erme_response::interpolate_cubic(double x, double x0, double x1,
double x2, double x3, double r0, double r1, double r2, double r3)

Interpolate a function evaluated at four points using cubic
interpolation.

This was generated using Maple's code generator with optimization. ";

%feature("docstring")  erme_response::interpolate "double
erme_response::interpolate(double xi, detran_utilities::vec_dbl &x,
detran_utilities::vec_dbl &r)

Interpolation helper function.

Note, we invert the abscissa, since interpolating responses on the
inverse of k reduces the error significantly. This was suggested in
Zhang & Rahnema (2012). While they only use linear interpolation, the
1/k dependences works for quadratic and cubic as well.

Parameters:
-----------

xi:  Value of independent variable at which to evaluate function

x:  Abscissa

r:  Function evaluated at abscissa ";


// File: namespaceerme__solver.xml


// File: namespaceerme__utils.xml


// File: namespacelinear__algebra.xml
%feature("docstring")  linear_algebra::initialize "void
linear_algebra::initialize(int &argc, char **&argv)

Initialize a parallel job. ";

%feature("docstring")  linear_algebra::finalize "void
linear_algebra::finalize()

Finish a parallel job. ";

%feature("docstring")  linear_algebra::multiply_wrapper "PetscErrorCode linear_algebra::multiply_wrapper(Mat A, Vec x, Vec y)

PETSc calls this function when applying the operator. The function
then calls the shell multiply function, when must be implemented by
the derived class. ";


// File: namespacepyerme.xml


// File: namespacepyserment.xml


// File: namespaceserment__comm.xml


// File: namespacestd.xml


// File: namespaceutil.xml
%feature("docstring")  util::inner_product "T
util::inner_product(const Vector_Lite< T, N > &a, const Vector_Lite<
T, N > &b)

Computes the inner product between two vectors.

Parameters:
-----------

a:  1st vector.

b:  2nd vector. ";


// File: ____init_____8py.xml


// File: AbsorptionOperator_8cc.xml


// File: AbsorptionOperator_8hh.xml


// File: AbsorptionResponse_8cc.xml


// File: AbsorptionResponse_8hh.xml


// File: CartesianNode_8cc.xml


// File: CartesianNode_8hh.xml


// File: CartesianNodeDetran_8cc.xml


// File: CartesianNodeDetran_8hh.xml


// File: Comm_8hh.xml


// File: Comm__Traits_8hh.xml


// File: Connect_8cc.xml


// File: Connect_8hh.xml


// File: Constants_8hh.xml


// File: DataBase_8hh.xml


// File: Definitions_8hh.xml


// File: diff2d_8cc.xml
%feature("docstring")  std::main "int main(int argc, char *args[]) ";


// File: Diff2dElement_8hh.xml


// File: Diff2dInput_8cc.xml


// File: Diff2dInput_8hh.xml


// File: Diff2dOutput_8cc.xml


// File: Diff2dOutput_8hh.xml


// File: Diff2dProblem_8cc.xml


// File: Diff2dProblem_8hh.xml


// File: Diff2dSolver_8cc.xml


// File: Diff2dSolver_8hh.xml


// File: DummyNode_8cc.xml


// File: DummyNode_8hh.xml


// File: EigenSolver_8cc.xml


// File: EigenSolver_8hh.xml


// File: EigenvalueUpdate_8hh.xml


// File: FissionOperator_8cc.xml


// File: FissionOperator_8hh.xml


// File: FissionResponse_8cc.xml


// File: FissionResponse_8hh.xml


// File: GlobalSolver_8cc.xml


// File: GlobalSolver_8hh.xml


// File: GlobalSolverBase_8cc.xml


// File: GlobalSolverBase_8hh.xml


// File: GlobalSolverNewton_8cc.xml


// File: GlobalSolverNewton_8hh.xml


// File: GlobalSolverPicard_8cc.xml


// File: GlobalSolverPicard_8hh.xml


// File: GlobalSolverPicard_8i_8hh.xml


// File: HexagonalNode_8hh.xml


// File: InnerIterBase_8cc.xml


// File: InnerIterBase_8hh.xml


// File: InnerIterPower_8cc.xml


// File: InnerIterPower_8hh.xml


// File: InnerIterSLEPc_8cc.xml


// File: InnerIterSLEPc_8hh.xml
%feature("docstring")  apply_MR "PetscErrorCode apply_MR(Mat A, Vec
X_in, Vec X_out)

A matrix-vector multiplication wrapper function for MR.

This is needed because PETSc works with function pointers, which does
not include pointers to member functions.

Parameters:
-----------

A:  PETSc shell matrix

X_in:  Incoming PETSc vector

X_out:  Outgoing PETSc vector ";


// File: InputXML_8cc.xml


// File: InputXML_8hh.xml


// File: Interpolation_8hh.xml


// File: InvItMatVecWrap_8hh.xml
%feature("docstring")  InvItMatVecWrap "PetscErrorCode
InvItMatVecWrap(Mat M, Vec X, Vec Y)

This is a matrix-vector multiplication wrapper function.

It exists outside the SermentMatrixShell class so that it can be
passed to Petsc, as Petsc needs pointers to functions, and it doesn't
like pointers to member methods. ";


// File: InvItShell_8cc.xml
%feature("docstring")  InvItMatVecWrap "PetscErrorCode
InvItMatVecWrap(Mat M, Vec f, Vec fpf)

This is a matrix-vector multiplication wrapper function.

It exists outside the SermentMatrixShell class so that it can be
passed to Petsc, as Petsc needs pointers to functions, and it doesn't
like pointers to member methods. ";


// File: InvItShell_8hh.xml


// File: Jacobian_8hh.xml


// File: JacobianEmpty_8hh.xml
%feature("docstring")  JacobianEmpty "PetscErrorCode
JacobianEmpty(SNES snes, Vec X, Mat *J, Mat *B, MatStructure *flag,
void *ptr) ";


// File: JacobianShell_8cc.xml
%feature("docstring")  myMatVecWrap "PetscErrorCode myMatVecWrap(Mat
M, Vec f, Vec fpf)

This is a matrix-vector multiplication wrapper function for use with
JacobianShell.

It exists outside the SermentMatrixShell class so that it can be
passed to Petsc, as Petsc needs pointers to functions but doesn't like
pointers to member functions. Here, the context associated with the
instance of JacobianShell is retrieved. That context is simply a
pointer to the global ( Newton) solver in disguise, so we cast it into
the correct form and call JacobianShell's own member function for
multiplication. ";


// File: JacobianShell_8hh.xml


// File: LeakageOperator_8cc.xml


// File: LeakageOperator_8hh.xml


// File: LeakageResponse_8cc.xml


// File: LeakageResponse_8hh.xml


// File: LegendrePoly_8cc.xml


// File: LegendrePoly_8hh.xml


// File: LinAlg_8hh.xml


// File: LinearAlgebraSetup_8hh.xml


// File: LinearSolver_8cc.xml


// File: LinearSolver_8hh.xml


// File: LocalProblem_8hh.xml


// File: LocalProblemDiff2d_8cc.xml


// File: LocalProblemDiff2d_8hh.xml


// File: ManagerERME_8cc.xml


// File: ManagerERME_8hh.xml


// File: Matrix_8cc.xml


// File: Matrix_8hh.xml


// File: matrix__shell__fixture_8hh.xml


// File: matrix__test_8cc.xml
%feature("docstring")  main "int main(int argc, char *args[]) ";


// File: MatrixBase_8cc.xml


// File: MatrixBase_8hh.xml


// File: MatrixBase_8i_8hh.xml


// File: MatrixShell_8cc.xml


// File: MatrixShell_8hh.xml


// File: MPI_8cc.xml


// File: MPI_8hh.xml


// File: mpi__matrix__test_8cc.xml
%feature("docstring")  buildmatrix "void buildmatrix(Mat &A, integer
N)

This builds the 1-D finite difference matrix. ";

%feature("docstring")  buildvector "void buildvector(Vec &x, Vec &b,
integer N)

This builds the right hand side. ";

%feature("docstring")  main "int main(int argc, char *args[])

This is a simple parallel parallel linear solver demonstration.

The system solved is y''(x) = 1 for x = [0,1], subject to y(0)=y(1)=0.
It is implemented as the linear system Ax = b, where A is the 1-D
second-difference operator defined as (1/h)^2 * [ 2 -1 0 ..., -1 2 -1
....] where h = 1/N, and N is the number of divisions. The right hand
side is a vector of 1's.

Uses as a reference: ksp/examples/tutorials/ex2.c ";


// File: MPI__Traits_8hh.xml


// File: MyMatVecWrap_8hh.xml
%feature("docstring")  myMatVecWrap "PetscErrorCode myMatVecWrap(Mat
M, Vec X, Vec Y)

This is a matrix-vector multiplication wrapper function for use with
JacobianShell.

It exists outside the SermentMatrixShell class so that it can be
passed to Petsc, as Petsc needs pointers to functions but doesn't like
pointers to member functions. Here, the context associated with the
instance of JacobianShell is retrieved. That context is simply a
pointer to the global ( Newton) solver in disguise, so we cast it into
the correct form and call JacobianShell's own member function for
multiplication. ";


// File: NeighborSurface_8hh.xml


// File: Newton_8cc.xml


// File: Newton_8hh.xml


// File: Node_8cc.xml


// File: Node_8hh.xml


// File: Node_8i_8hh.xml


// File: node__fixture_8hh.xml


// File: NodeFactory_8hh.xml


// File: NodeFactoryDetran_8cc.xml


// File: NodeFactoryDetran_8hh.xml


// File: NodeList_8cc.xml


// File: NodeList_8hh.xml


// File: NodeList_8i_8hh.xml


// File: nodelist__fixture_8hh.xml


// File: NodePartitioner_8cc.xml


// File: NodePartitioner_8hh.xml


// File: NodeResponse_8cc.xml


// File: NodeResponse_8hh.xml


// File: NodeResponse_8i_8hh.xml


// File: NodeSerialization_8hh.xml


// File: NonlinearResidual_8hh.xml


// File: NonlinearResidual_8i_8hh.xml


// File: OperatorMR_8hh.xml


// File: PCAppxJacobian_8cc.xml


// File: PCAppxJacobian_8hh.xml


// File: PicardInnerBase_8hh.xml


// File: PowerIter_8cc.xml


// File: PowerIter_8hh.xml


// File: ResidualWrap_8hh.xml
%feature("docstring")  ResidualWrap "PetscErrorCode ResidualWrap(SNES
snes, Vec X, Vec F, void *ptr) ";


// File: ResponseDatabase_8cc.xml


// File: ResponseDatabase_8hh.xml


// File: ResponseDatabase_8i_8hh.xml


// File: ResponseFunction_8cc.xml


// File: ResponseFunction_8hh.xml


// File: ResponseFunctionDiffusion_8cc.xml


// File: ResponseFunctionDiffusion_8hh.xml


// File: ResponseFunctionServer_8cc.xml


// File: ResponseFunctionServer_8hh.xml


// File: ResponseIndex_8hh.xml


// File: ResponseIndexer_8cc.xml


// File: ResponseIndexer_8hh.xml


// File: ResponseIndexer_8i_8hh.xml


// File: response_2ResponseMatrix_8cc.xml


// File: ResponseMatrix_8cc.xml


// File: response_2ResponseMatrix_8hh.xml


// File: ResponseMatrix_8hh.xml


// File: ResponseMatrixFull_8cc.xml


// File: ResponseMatrixFull_8hh.xml


// File: response_2ResponseOperator_8hh.xml


// File: ResponseOperator_8hh.xml


// File: ResponseServer_8cc.xml


// File: ResponseServer_8hh.xml


// File: ResponseServer_8i_8hh.xml


// File: ResponseSource_8hh.xml


// File: ResponseSourceDatabase_8cc.xml


// File: ResponseSourceDatabase_8hh.xml


// File: ResponseSourceDetran_8hh.xml


// File: ResponseSourceDummy_8cc.xml


// File: ResponseSourceDummy_8hh.xml


// File: ResponseSourceFactory_8hh.xml


// File: ResponseSourceFactoryDatabase_8hh.xml


// File: ResponseSourceFactoryDetran_8hh.xml


// File: ResponseSourceFactoryDummy_8hh.xml


// File: Serial_8cc.xml


// File: Serial_8hh.xml


// File: SermentMatrix_8hh.xml


// File: SermentMatrixBCRS_8cc.xml


// File: SermentMatrixBCRS_8hh.xml


// File: SermentMatrixCRS_8cc.xml


// File: SermentMatrixCRS_8hh.xml


// File: SermentMatrixFull_8cc.xml


// File: SermentMatrixFull_8hh.xml


// File: SermentMatrixShell_8cc.xml


// File: SermentMatrixShell_8hh.xml


// File: SermentVector_8cc.xml


// File: SermentVector_8hh.xml


// File: StateERME_8cc.xml


// File: StateERME_8hh.xml


// File: StateERME_8i_8hh.xml


// File: SteffensenUpdate_8hh.xml


// File: test__AbsorptionOperator_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_AbsorptionOperator "int
test_AbsorptionOperator(int argc, char *argv[]) ";


// File: test__CartesianNodeDetran_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_CartesianNodeDetran "int
test_CartesianNodeDetran(int argc, char *argv[]) ";

%feature("docstring")  test_CartesianNodeDetran_serialize "int
test_CartesianNodeDetran_serialize(int argc, char *argv[]) ";


// File: test__Communicator_8cc.xml
%feature("docstring")  detran_test::main "int main(int argc, char
*argv[]) ";

%feature("docstring")  detran_test::test_Communicator "int
test_Communicator(int argc, char *argv[]) ";


// File: test__Connect_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Connect "int test_Connect(int argc, char
*argv[]) ";


// File: test__FissionOperator_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_FissionOperator "int
test_FissionOperator(int argc, char *argv[]) ";


// File: test__GlobalReduction_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_GlobalReduction "int
test_GlobalReduction(int argc, char *argv[]) ";


// File: test__GlobalSolverPicard_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_GlobalSolverPicard "int
test_GlobalSolverPicard(int argc, char *argv[]) ";


// File: test__Interpolation_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Interpolation "int test_Interpolation(int
argc, char *argv[]) ";


// File: test__LeakageOperator_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_LeakageOperator "int
test_LeakageOperator(int argc, char *argv[]) ";


// File: test__Matrix_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Matrix_actual "int test_Matrix_actual()
";

%feature("docstring")  test_Matrix "int test_Matrix(int argc, char
*argv[]) ";


// File: test__MatrixShell_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_MatrixShell_actual "int
test_MatrixShell_actual() ";

%feature("docstring")  test_MatrixShell "int test_MatrixShell(int
argc, char *argv[]) ";


// File: test__MR__scaling_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_MR_scaling "int test_MR_scaling(int argc,
char *argv[]) ";


// File: test__NodeList_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_NodeList "int test_NodeList(int argc,
char *argv[]) ";


// File: test__NodePartitioner_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_NodePartitioner "int
test_NodePartitioner(int argc, char *argv[]) ";


// File: test__NodeResponse_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_NodeResponse "int test_NodeResponse(int
argc, char *argv[]) ";


// File: test__PingPong_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_PingPong_sendreceive "int
test_PingPong_sendreceive(int argc, char *argv[]) ";

%feature("docstring")  test_PingPong_bandwidth "int
test_PingPong_bandwidth(int argc, char *argv[]) ";

%feature("docstring")  test_PingPong_latency "int
test_PingPong_latency(int argc, char *argv[]) ";


// File: test__ResponseDatabase_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_ResponseDatabase "int
test_ResponseDatabase(int argc, char *argv[]) ";


// File: test__ResponseIndexer_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_ResponseIndexer "int
test_ResponseIndexer(int argc, char *argv[]) ";

%feature("docstring")  test_ResponseIndexer_zeroth_order "int
test_ResponseIndexer_zeroth_order(int argc, char *argv[]) ";


// File: test__ResponseMatrix_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_ResponseMatrix "int
test_ResponseMatrix(int argc, char *argv[]) ";


// File: test__ResponseServer_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_ResponseServer "int
test_ResponseServer(int argc, char *argv[]) ";


// File: test__StateERME_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_StateERME "int test_StateERME(int argc,
char *argv[]) ";


// File: test__Vector_8cc.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";

%feature("docstring")  test_Vector_actual "int test_Vector_actual()
";

%feature("docstring")  test_Vector "int test_Vector(int argc, char
*argv[]) ";


// File: TwoGroupKernel_8cc.xml


// File: TwoGroupKernel_8hh.xml


// File: TwoGroupKernelSolver_8cc.xml


// File: TwoGroupKernelSolver_8hh.xml


// File: typedefs_8hh.xml


// File: Vector_8cc.xml


// File: Vector_8hh.xml


// File: Vector_8i_8hh.xml


// File: Vector__Lite_8hh.xml


// File: Vector__Lite_8i_8hh.xml


// File: vector__test_8cc.xml
%feature("docstring")  main "int main(int argc, char *args[]) ";


// File: xerme_8cc.xml
%feature("docstring")  print_welcome "void print_welcome() ";

%feature("docstring")  main "int main(int argc, char **argv) ";


// File: Serment.xml


// File: todo.xml


// File: dir_e01a5fc97b175d2a37c8aa7e41f81f49.xml


// File: dir_aec6c3cec0d5ab4596e7af8b48b8c357.xml


// File: dir_05c849ba6d5972d285e2dc0a8ce6c25d.xml


// File: dir_cf3570b8581645c72ea1c163c2803a80.xml


// File: dir_40607dea5cea8ff12f9efda32114847f.xml


// File: dir_ef4b33cf5aff78c628b5a529476b31c3.xml


// File: dir_3c19d26104445ff48ac374210665d76d.xml


// File: dir_99b865f0751ba1cbccfe2eca18e968ba.xml


// File: dir_528e51a2644c2cdfc562d45df0e1bd6e.xml


// File: dir_38bb8b0840704712ec725cd151201e4d.xml


// File: dir_d633216edf99b850b786ebd5be677684.xml


// File: dir_ba05b86b9ac25edd42943bdc53fac6e4.xml


// File: dir_349f70b507179a51ff8e4cca288d2e33.xml


// File: dir_0830de72ed239dc531f2ebd2d5415313.xml


// File: dir_7cfd8e74fcab8e41b93963e1f3012ad6.xml


// File: dir_b5d6fad520a19978fa347f97d81bb5c3.xml


// File: dir_2fc1cbdeaae777e160b5fc1226da74a4.xml


// File: dir_c7588ae18addb1e9248d05d3049ba339.xml


// File: dir_845728f9def611742969df60838a1361.xml


// File: dir_733b01c174046f1c06df7d9f3f9c99ac.xml


// File: dir_fbd3184be27cc6d2a4277ca94f8cd754.xml


// File: dir_8103e6d399404576baedf3dccc30d6f0.xml


// File: erme_2test_2test_AbsorptionOperator-example.xml


// File: erme_2test_2test_FissionOperator-example.xml


// File: erme_2test_2test_LeakageOperator-example.xml


// File: erme_2test_2test_ResponseMatrix-example.xml


// File: erme_2test_2test_StateERME_8cc-example.xml


// File: erme_geometry_2test_2test_NodeList_8cc-example.xml


// File: erme_response_2test_2test_NodeResponse_8cc-example.xml


// File: erme_response_2test_2test_ResponseIndexer_8cc-example.xml


// File: erme_response_2test_2test_ResponseServer_8cc-example.xml


// File: linear_algebra_2test_2test_Vector_8cc-example.xml


// File: utils_2test_2tstVector_Lite_8cc-example.xml


// File: indexpage.xml

