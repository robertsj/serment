#------------------------------------------------------------------------------#
# serment/src/linear_algebra/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding serment linear algebra tests." )

ADD_EXECUTABLE(test_Vector                  test_Vector.cc)
TARGET_LINK_LIBRARIES(test_Vector           linear_algebra )
ADD_TEST(test_Vector                        test_Vector 0)

ADD_EXECUTABLE(test_Matrix                  test_Matrix.cc)
TARGET_LINK_LIBRARIES(test_Matrix           linear_algebra)
ADD_TEST(test_Matrix                        test_Matrix 0)

ADD_EXECUTABLE(test_MatrixShell             test_MatrixShell.cc)
TARGET_LINK_LIBRARIES(test_MatrixShell      linear_algebra)
ADD_TEST(test_MatrixShell                   test_MatrixShell 0)

ADD_EXECUTABLE(test_NonlinearSolver         test_NonlinearSolver.cc)
TARGET_LINK_LIBRARIES(test_NonlinearSolver  linear_algebra)
ADD_TEST(test_NonlinearSolver               test_NonlinearSolver 0)