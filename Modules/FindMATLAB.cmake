SET(MATLAB_INCLUDE_SEARCH_PATH
    /usr/local/MATLAB/R2020b/extern/include
)

SET(MATLAB_LIBRARY_SEARCH_PATH
    /usr/local/MATLAB/R2020b/bin/glnxa64
)

FIND_PATH(MATLAB_H    mat.h    ${MATLAB_INCLUDE_SEARCH_PATH})
FIND_LIBRARY(MATLAB_MAT_L    NAMES mat    PATHS ${MATLAB_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(MATLAB_MX_L    NAMES mx    PATHS ${MATLAB_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(MATLAB_MEX_L    NAMES mex    PATHS ${MATLAB_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(MATLAB_ENG_L    NAMES eng    PATHS ${MATLAB_LIBRARY_SEARCH_PATH})

SET(MATLAB_FOUND 1)

FOREACH(var MATLAB_H MATLAB_MAT_L MATLAB_MX_L MATLAB_MEX_L MATLAB_ENG_L)
    IF(NOT ${var})
    SET(MATLAB_FOUND 0)
    ENDIF(NOT ${var})
ENDFOREACH(var)

IF(MATLAB_FOUND)
  SET(MATLAB_INCLUDE_DIR  ${MATLAB_H})
  # MESSAGE(${MATLAB_INCLUDE_DIR}"+++++++++++")
  SET(MATLAB_LIBRARIES    ${MATLAB_MAT_L} ${MATLAB_MX_L} ${MATLAB_MEX_L} ${MATLAB_ENG_L})
ENDIF(MATLAB_FOUND)