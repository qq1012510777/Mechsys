#####################################################################################
# MechSys - A C++ library to simulate Mechanical Systems                            #
# Copyright (C) 2010 Sergio Galindo                                                 #
#                                                                                   #
# This file is part of MechSys.                                                     #
#                                                                                   #
# MechSys is free software; you can redistribute it and/or modify it under the      #
# terms of the GNU General Public License as published by the Free Software         #
# Foundation; either version 2 of the License, or (at your option) any later        #
# version.                                                                          #
#                                                                                   #
# MechSys is distributed in the hope that it will be useful, but WITHOUT ANY        #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.          #
#                                                                                   #
# You should have received a copy of the GNU General Public License along with      #
# MechSys; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, #
# Fifth Floor, Boston, MA 02110-1301, USA                                           #
#####################################################################################

SET(FADE2D_INCLUDE_SEARCH_PATH
  /usr/include
  /usr/local/include
  $ENV{MECHSYS_ROOT}/fadeRelease_v1.84/include_fade2d
  $ENV{HOME}/pkg/fadeRelease_v1.84/include_fade2d
  )

SET(FADE2D_LIBRARY_SEARCH_PATH
  /usr/lib
  /usr/local/lib
  $ENV{MECHSYS_ROOT}/fadeRelease_v1.84/lib_ubuntu18.04_x86_64
  $ENV{HOME}/pkg/fadeRelease_v1.84/lib_ubuntu18.04_x86_64
  )

FIND_PATH(FADE2D_H Fade_2D.h ${FADE2D_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(FADE2D_LIB NAMES fade2d PATHS ${FADE2D_LIBRARY_SEARCH_PATH})

SET(FADE2D_FOUND 1)
FOREACH(var FADE2D_H FADE2D_LIB)
  IF(NOT ${var})
	SET(FADE2D_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(FADE2D_FOUND)
  SET(FADE2D_INCLUDE_DIRS ${FADE2D_H})
  SET(FADE2D_LIBRARIES    ${FADE2D_LIB})
ENDIF(FADE2D_FOUND)
