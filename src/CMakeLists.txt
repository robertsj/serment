# Includes
INCLUDE_DIRECTORIES(linear_algebra)
INCLUDE_DIRECTORIES(comm)
INCLUDE_DIRECTORIES(erme_geometry)
INCLUDE_DIRECTORIES(erme_response)
INCLUDE_DIRECTORIES(erme)
INCLUDE_DIRECTORIES(erme_solvers)

# Source subdirectories.
ADD_SUBDIRECTORY(linear_algebra)
ADD_SUBDIRECTORY(comm)
ADD_SUBDIRECTORY(erme_geometry)
ADD_SUBDIRECTORY(erme_response)
ADD_SUBDIRECTORY(erme)
ADD_SUBDIRECTORY(erme_solvers)

