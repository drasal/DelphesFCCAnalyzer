#############################################################################
#
# CMAKE uninstall module
#
#############################################################################

# Check that software was installed
IF( NOT EXISTS "@CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt" )
  MESSAGE( FATAL_ERROR "Cannot find install manifest: @CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt" )
ENDIF( NOT EXISTS "@CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt" )

# Read install_manifest.txt and uninstall all what was automatically installed
FILE( READ "@CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt" files )

STRING( REGEX REPLACE "\n" ";" files "${files}" )
FOREACH( file ${files} )
  MESSAGE( STATUS "Uninstalling $ENV{DESTDIR}${file}" )
  IF( IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}" )
    EXEC_PROGRAM(
      "@CMAKE_COMMAND@" ARGS "-E remove \"$ENV{DESTDIR}${file}\""
      OUTPUT_VARIABLE rm_out
      RETURN_VALUE rm_retval
      )
    IF( NOT "${rm_retval}" STREQUAL 0)
      MESSAGE( FATAL_ERROR "Problem when removing $ENV{DESTDIR}${file}" )
    ENDIF(NOT "${rm_retval}" STREQUAL 0)
  ELSE(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    MESSAGE(STATUS "File $ENV{DESTDIR}${file} does not exist.")
  ENDIF(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
ENDFOREACH( file )

# Remove symlink
MESSAGE( STATUS "Uninstalling symlink $ENV{HOME}/bin/analyzeDelphes" )
IF( IS_SYMLINK "$ENV{HOME}/bin/analyzeDelphes")
  EXEC_PROGRAM(
         "@CMAKE_COMMAND@" ARGS "-E remove \"$ENV{HOME}/bin/analyzeDelphes\""
         OUTPUT_VARIABLE rm_out
         RETURN_VALUE rm_retval
         )
  IF( NOT "${rm_retval}" STREQUAL 0)
      MESSAGE( FATAL_ERROR "Problem when removing symlink $ENV{HOME}/bin/analyzeDelphes" )
  ENDIF(NOT "${rm_retval}" STREQUAL 0)
ELSE()
MESSAGE(STATUS "Symlink $ENV{HOME}/bin/analyzeDelphes does not exist.")
ENDIF()
