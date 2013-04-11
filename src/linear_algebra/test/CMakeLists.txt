#------------------------------------------------------------------------------#
# serment/src/linear_algebra/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding serment linear algebra tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_Vector                test_Vector.cc)
TARGET_LINK_LIBRARIES(test_Vector         linear_algebra )

ADD_EXECUTABLE(test_Matrix                test_Matrix.cc)
TARGET_LINK_LIBRARIES(test_Matrix         linear_algebra)

ADD_EXECUTABLE(test_MatrixShell           test_MatrixShell.cc)
TARGET_LINK_LIBRARIES(test_MatrixShell    linear_algebra)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_Vector test_Vector 0)
ADD_TEST(test_Matrix test_Matrix 0)
ADD_TEST(test_MatrixShell test_MatrixShell 0)