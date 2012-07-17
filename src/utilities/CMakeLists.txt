#-------------------------------------------------------------------------------
# serment/src/utilities/CMakeLists.txt
#
# $Rev:: 103                                           $:Rev of last commit
# $Author:: j.alyn.roberts@gma                         $:Author of last commit
# $Date:: 2011-06-14 22:22:46 -0400 (Tue, 14 Jun 2011) $:Date of last commit  
#
#-------------------------------------------------------------------------------

# Source subdirectories
SET( UTILS_SRC 
     LegendrePoly.cc
     GenException.cc
     InputXML.cc
     DBC.cc
)

# Make the library
ADD_LIBRARY( utils 
             ${UTILS_SRC} 
)
TARGET_LINK_LIBRARIES( utils
                       ${LIBXML2_LIBRARIES}
)

