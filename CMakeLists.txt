#############################################################################
#
# CMAKE file for building FCC analyzer software
#
#############################################################################

#----------------------------------------------------------------------------
# CMAKE compatibility issues
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

#----------------------------------------------------------------------------
# Project name and version
PROJECT(DelphesFCCAnalyzer)

#set( ${PROJECT_NAME}_MAJOR_VERSION 1 )
#set( ${PROJECT_NAME}_MINOR_VERSION 0 )

#----------------------------------------------------------------------------
# Project options
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )

#----------------------------------------------------------------------------
# Advanced settings
# Library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION
    "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}" )

# Output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# Set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" )
ENDIF()


#----------------------------------------------------------------------------
# Project dependencies

#
# Get boost include & lib directory
#
SET( BOOST_INCLUDE_DIR "$ENV{BOOST_INCLUDE}" )

# Find all the libraries and save them 
SET( BOOST_LIBS_NAMES "boost_system" 
                      "boost_filesystem" 
                      "boost_program_options" )
SET( BOOST_LIBS_DIR "$ENV{BOOST_LIB}" )

FOREACH( libname ${BOOST_LIBS_NAMES} )

    SET( BOOST_LIB_${libname} BOOST_LIB_${libname}-NOTFOUND )
    MARK_AS_ADVANCED( BOOST_LIB_${libname} )
    #MESSAGE( STATUS "${libname} >> ${BOOST_LIBS_DIR}")
    FIND_LIBRARY( BOOST_LIB_${libname} 
        NAMES ${libname} 
        PATHS ${BOOST_LIBS_DIR} 
        PATH_SUFFIXES lib )
    #MESSAGE( STATUS ">>> ${BOOST_LIB_${libname}}" )
    IF( NOT BOOST_LIB_${libname} )
        SET( BOOST_FINDLIB_FAILED TRUE )
        IF( NOT BOOST_FIND_QUIETLY )
            MESSAGE( STATUS "Check for BOOST: ${BOOST_LIB_DIR}" 
                            " -- failed to find BOOST ${libname} library!!" )
        ENDIF()
    ELSE()
        LIST( APPEND BOOST_LIBS ${BOOST_LIB_${libname}} )
    ENDIF()
ENDFOREACH()

#
# Get podio include & lib directory
#
SET( PODIO_INCLUDE_DIR "$ENV{PODIO}/include" )

# Find all the libraries and save them 
SET( PODIO_LIBS_NAMES "podio" 
                      "podioDict" ) 
#                      "podioDict" ) 
SET( PODIO_LIBS_DIR "$ENV{PODIO}/lib" )

FOREACH( libname ${PODIO_LIBS_NAMES} )

    SET( PODIO_LIB_${libname} PODIO_LIB_${libname}-NOTFOUND )
    MARK_AS_ADVANCED( PODIO_LIB_${libname} )
    #MESSAGE( STATUS "${libname} >> ${PODIO_LIBS_DIR}")
    FIND_LIBRARY( PODIO_LIB_${libname} 
        NAMES ${libname} 
        PATHS ${PODIO_LIBS_DIR} 
        PATH_SUFFIXES lib )
    MESSAGE( STATUS "Loading podio library: ${PODIO_LIB_${libname}}" )
    IF( NOT PODIO_LIB_${libname} )
        SET( PODIO_FINDLIB_FAILED TRUE )
        IF( NOT PODIO_FIND_QUIETLY )
            MESSAGE( STATUS "Check for PODIO: ${PODIO_LIB_DIR}" 
                            " -- failed to find PODIO ${libname} library!!" )
        ENDIF()
    ELSE()
       LIST( APPEND PODIO_LIBS ${PODIO_LIB_${libname}} )
    ENDIF()
ENDFOREACH()

#
# Get fcc-edm include & lib directory
#
SET( FCCEDM_INCLUDE_DIR "$ENV{FCCEDM}/include" )

# Find all the libraries and save them 
SET( FCCEDM_LIBS_NAMES "datamodel" 
                       "datamodelDict" )
SET( FCCEDM_LIBS_DIR "$ENV{FCCEDM}/lib" )

FOREACH( libname ${FCCEDM_LIBS_NAMES} )

    SET( FCCEDM_LIB_${libname} FCCEDM_LIB_${libname}-NOTFOUND )
    MARK_AS_ADVANCED( FCCEDM_LIB_${libname} )
    #MESSAGE( STATUS "${libname} >> ${FCCEDM_LIBS_DIR}")
    FIND_LIBRARY( FCCEDM_LIB_${libname} 
        NAMES ${libname} 
        PATHS ${FCCEDM_LIBS_DIR} 
        PATH_SUFFIXES lib )
    MESSAGE( STATUS "Loading fcc-edm library: ${FCCEDM_LIB_${libname}}" )
    IF( NOT FCCEDM_LIB_${libname} )
        SET( FCCEDM_FINDLIB_FAILED TRUE )
        IF( NOT FCCEDM_FIND_QUIETLY )
            MESSAGE( STATUS "Check for FCCEDM: ${FCCEDM_LIB_DIR}" 
                            " -- failed to find FCCEDM ${libname} library!!" )
        ENDIF()
    ELSE()
       LIST( APPEND FCCEDM_LIBS ${FCCEDM_LIB_${libname}} )
    ENDIF()
ENDFOREACH()

#
# Get include directory using root-config
#
SET( ROOT_HOME "$ENV{ROOT_CONFIG}" )
EXEC_PROGRAM( "${ROOT_HOME}/root-config" "."
    ARGS --incdir
    OUTPUT_VARIABLE out_var
    RETURN_VALUE out_ret )

# Check if everything went OK
IF( out_ret)
    IF( NOT ROOT_FIND_QUIETLY )
        MESSAGE( STATUS "Check for ROOT: ${ROOT_HOME}/root-config --incdir"
                        " -- failed to find ROOT include directory!!")
    ENDIF()
    SET( ROOT_FINDINCL_FAILED TRUE )
    MARK_AS_ADVANCED( ROOT_FINDINCL_FAILED )
ELSE()
    SET( ROOT_FINDINCL_FAILED FALSE )
    SET( ROOT_INCLUDE_DIRS ${out_var} )
    MARK_AS_ADVANCED( ROOT_FINDINCL_FAILED ROOT_INCLUDE_DIRS )
ENDIF()


#
# Get libraries using root-config
#
SET( ROOT_LIBS_NAMES )

EXEC_PROGRAM( "${ROOT_HOME}/root-config" "."
    ARGS --noauxlibs --evelibs
    OUTPUT_VARIABLE out_var
    RETURN_VALUE out_ret )

# Check if everything went OK 
IF( out_ret) 
    IF( NOT ROOT_FIND_QUIETLY ) 
        MESSAGE( STATUS "Check for ROOT: ${ROOT_HOME}/root-config --noauxlibs --evelibs" 
                        " -- failed to find ROOT libraries!!") 
    ENDIF()                      
    SET( ROOT_FINDLIBS_FAILED TRUE ) 
    MARK_AS_ADVANCED( ROOT_FINDLIB_FAILED ) 
ELSE() 
    # Each argument will be separated by semicolon 
    SEPARATE_ARGUMENTS( out_var ) 
     
    # Get library path and remove -L compiler flag 
    LIST( GET out_var 0 libdir ) 
    LIST( REMOVE_AT out_var 0 ) 
    STRING( REGEX REPLACE "^-.(.*)$" "\\1" ROOT_LIBS_DIR "${libdir}") 
    MARK_AS_ADVANCED( ROOT_LIBS_DIR ) 
     
    # Extract library names 
    FOREACH( lib ${out_var} ) 
        STRING( REGEX REPLACE "^-.(.*)$" "\\1" libname "${lib}") 
        LIST( APPEND ROOT_LIBS_NAMES ${libname} ) 
    ENDFOREACH() 
     
    # Find all the libraries and save them 
    FOREACH( libname ${ROOT_LIBS_NAMES} ) 
     
        SET( ROOT_LIB_${libname} ROOT_LIB_${libname}-NOTFOUND ) 
        MARK_AS_ADVANCED( ROOT_LIB_${libname} ) 
 
        FIND_LIBRARY( ROOT_LIB_${libname} 
            NAMES ${libname} 
            PATHS ${ROOT_LIBS_DIR} 
            PATH_SUFFIXES lib ) 
 
        IF( NOT ROOT_LIB_${libname} ) 
            SET( ROOT_FINDLIB_FAILED TRUE ) 
            IF( NOT ROOT_FIND_QUIETLY ) 
                MESSAGE( STATUS "Check for ROOT: ${ROOT_HOME}/root-config" 
                                " -- failed to find ROOT ${libname} library!!" ) 
            ENDIF() 
        ELSE() 
            #MESSAGE(STATUS ">>> ${ROOT_LIB_${libname}}")
            LIST( APPEND ROOT_LIBS ${ROOT_LIB_${libname}} ) 
        ENDIF() 
    ENDFOREACH() 
ENDIF() 

#----------------------------------------------------------------------------
# Include, source, dirs 
#
# Compiler flags
#add_definitions( "-Wall -Wno-long-long -std=c++11 -pedantic" )
SET ( CMAKE_CXX_COMPILER "g++" )
ADD_DEFINITIONS( "-Wall -std=c++11 -g -pedantic" )
#SET ( CMAKE_EXE_LINKER_FLAGS "-Wl" )

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include 
                     ${BOOST_INCLUDE_DIR}
                     ${PODIO_INCLUDE_DIR}
                     ${FCCEDM_INCLUDE_DIR}
                     ${ROOT_INCLUDE_DIRS} )
FILE( GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc )
FILE( GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h )

#----------------------------------------------------------------------------
# add the executable, and link it to the geant4 libraries
#
#SET(CMAKE_INSTALL_RPATH ${BOOST_LIBS_DIR} ${PYTHIA_LIBS_DIR} ${ROOT_LIBS_DIR})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

ADD_EXECUTABLE(analyzeDelphesFCC ${sources} ${headers} )
# explicitly say that the executable depends on custom target

TARGET_LINK_LIBRARIES(analyzeDelphesFCC ${BOOST_LIBS} ${PODIO_LIBS} ${FCCEDM_LIBS} ${ROOT_LIBS})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
INSTALL(TARGETS analyzeDelphesFCC DESTINATION bin)

IF(CMAKE_HOST_UNIX)
    INSTALL( CODE "MESSAGE(STATUS \"Creating symlink $ENV{HOME}/bin/analyzeDelphesFCC to ${PROJECT_SOURCE_DIR}/bin/analyzeDelphesFCC ... \" )
                   EXECUTE_PROCESS(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink ${PROJECT_SOURCE_DIR}/bin/analyzeDelphesFCC $ENV{HOME}/bin/analyzeDelphesFCC)" )
ENDIF()

#----------------------------------------------------------------------------
# Uninstall 
#
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeUninstall.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
