#------------------------------------------------------------------------------#
# serment/src/erme_solvers/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding serment erme_solvers tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_NonlinearResidual               test_NonlinearResidual.cc)
TARGET_LINK_LIBRARIES(test_NonlinearResidual        erme_solver)
ADD_TEST(test_NonlinearResidual                     test_NonlinearResidual 0)

ADD_EXECUTABLE(test_GlobalSolverPicard				test_GlobalSolverPicard.cc)
TARGET_LINK_LIBRARIES(test_GlobalSolverPicard 		erme_solver)
ADD_TEST(test_GlobalSolverPicard                    test_GlobalSolverPicard 0)

ADD_EXECUTABLE(test_GlobalSolverNewton              test_GlobalSolverNewton.cc)
TARGET_LINK_LIBRARIES(test_GlobalSolverNewton       erme_solver)
ADD_TEST(test_GlobalSolverNewton                    test_GlobalSolverNewton 0)
