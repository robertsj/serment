#------------------------------------------------------------------------------#
# serment/src/erme_solvers/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding serment erme_solvers tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_GlobalSolverPicard				test_GlobalSolverPicard.cc)
TARGET_LINK_LIBRARIES(test_GlobalSolverPicard 		erme_solver)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_GlobalSolverPicard 			test_GlobalSolverPicard 0)

