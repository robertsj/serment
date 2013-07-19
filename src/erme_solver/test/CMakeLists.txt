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

ADD_EXECUTABLE(test_Jacobian                        test_Jacobian.cc)
TARGET_LINK_LIBRARIES(test_Jacobian                 erme_solver)
ADD_TEST(test_Jacobian                              test_Jacobian 0)

ADD_EXECUTABLE(test_GlobalSolverPicard				test_GlobalSolverPicard.cc)
TARGET_LINK_LIBRARIES(test_GlobalSolverPicard 		erme_solver)
ADD_TEST(test_GlobalSolverPicard                    test_GlobalSolverPicard 0)

ADD_EXECUTABLE(test_GlobalSolverNewton              test_GlobalSolverNewton.cc)
TARGET_LINK_LIBRARIES(test_GlobalSolverNewton       erme_solver)
ADD_TEST(test_GlobalSolverNewton                    test_GlobalSolverNewton 0)

ADD_EXECUTABLE(test_EigenvalueUpdate                test_EigenvalueUpdate.cc)
TARGET_LINK_LIBRARIES(test_EigenvalueUpdate         erme_solver)
ADD_TEST(test_EigenvalueUpdate                      test_EigenvalueUpdate 0)
ADD_TEST(test_SteffenseneUpdate                     test_EigenvalueUpdate 1)
ADD_TEST(test_AnghelUpdateEXP                       test_EigenvalueUpdate 2)
ADD_TEST(test_AnghelUpdateLINLIN                    test_EigenvalueUpdate 2)
ADD_TEST(test_AnghelUpdateINVLIN                    test_EigenvalueUpdate 3)
ADD_TEST(test_AnghelUpdateINVINV                    test_EigenvalueUpdate 4)
