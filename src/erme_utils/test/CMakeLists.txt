#------------------------------------------------------------------------------#
# serment/src/erme/test/CMakeLists.txt
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_Archive				    test_Archive.cc)
TARGET_LINK_LIBRARIES(test_Archive          erme_utils)
ADD_TEST(test_Archive                       test_Archive 0)
