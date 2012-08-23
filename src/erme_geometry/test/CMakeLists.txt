#------------------------------------------------------------------------------#
# serment/src/erme_geometry/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding serment erme geometry tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_CartesianNodeDetran			test_CartesianNodeDetran.cc)
TARGET_LINK_LIBRARIES(test_CartesianNodeDetran 	erme_geometry)

ADD_EXECUTABLE(test_NodeList			test_NodeList.cc)
TARGET_LINK_LIBRARIES(test_NodeList 	erme_geometry)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_CartesianNodeDetran test_CartesianNodeDetran 0)

ADD_TEST(test_NodeList test_NodeList 0)
