cmake_minimum_required( VERSION 2.8.8 )

project( graphilt )

# VERSION
set( ${PROJECT_NAME}_MAJOR_VERSION "0" )
set( ${PROJECT_NAME}_MINOR_VERSION "1" )
set( ${PROJECT_NAME}_PATCH_VERSION "0" )
set( ${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}_PATCH_VERSION}" )
set( PROJECT_VERSION ${${PROJECT_NAME}_VERSION} )
set( PROJECT_SOVERSION ${${PROJECT_NAME}_MAJOR_VERSION} )

# Always include the source and build directories in the include path
# to save doing so manually in every subdirectory.
SET( CMAKE_INCLUDE_CURRENT_DIR ON )
# CMAKE_MODULE_PATH
#set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/" )
# 3rd Party dir
set( 3RDPARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty" )

#### OPTIONS
option( BUILD_TESTS "Build the test suite" ON )
option( BUILD_TOOLS "Build standard tools" ON )
option( BUILD_DOC "Build documentation" ON )

#### LIBS
# Boost
# set( Boost_USE_STATIC_LIBS ON )
# set( Boost_USE_MULTITHREADED ON ) 
# find_package( Boost COMPONENTS thread filesystem system graph serialization REQUIRED )
# include_directories( ${Boost_INCLUDE_DIRS} )

# find_package( Threads REQUIRED )

# if( WIN32 AND __COMPILER_GNU )
#     # mingw-gcc fails to link boost::thread
#     add_definitions( -DBOOST_THREAD_USE_LIB )
# endif( WIN32 AND __COMPILER_GNU )

# Readline
#find_package( ReadLine REQUIRED )
# Custom readline OS X is too old ...
# set( ReadLine_DIR /usr/local/Cellar/readline/6.2.4 )
# set( ReadLine_INCLUDE_DIR ${ReadLine_DIR}/include )
# set( ReadLine_LIB_DIR ${ReadLine_DIR}/lib )
# include_directories( ${ReadLine_INCLUDE_DIR} )
# link_directories( ${ReadLine_LIB_DIR} )

#### GENERAL OPTIONS
set( LIB_SUFFIX "" CACHE STRING "Define suffix of directory name (32/64)" )
set( EXEC_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE PATH "Installation prefix for executables and object code libraries" FORCE )
set( BIN_INSTALL_DIR ${EXEC_INSTALL_PREFIX}/bin CACHE PATH "Installation prefix for user executables" FORCE )
set( LIB_INSTALL_DIR ${EXEC_INSTALL_PREFIX}/lib${LIB_SUFFIX} CACHE PATH  "Installation prefix for object code libraries" FORCE )
set( INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include CACHE PATH "Installation prefix for C++ header files" FORCE )
set( DOC_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/doc CACHE PATH "Installation prefix for documentation" FORCE )

# Configure output bin path
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib )

# Check of resources folder exists, if not create it
set( RES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/resources" )
if( NOT EXISTS ${RES_DIR} )
    file( MAKE_DIRECTORY ${RES_DIR} )
endif()

if( CMAKE_BUILD_TYPE STREQUAL "Release" )
    add_definitions( -DQT_NO_DEBUG )
    add_definitions( -DQT_NO_DEBUG_OUTPUT )
    add_definitions( -O3 )
elseif( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    add_definitions( -DQT_DEBUG )
    add_definitions( -DGZK_DEBUG )
endif()

#### COMPILER SPECIFICS
message( STATUS "Using ${CMAKE_CXX_COMPILER_ID}" )
add_definitions( -Wall -Wextra -Wno-unused-parameter -Wdeprecated-declarations )
# C++11
add_definitions(-std=c++11) #-std=c++11 std=gnu++11
add_definitions(-stdlib=libc++) #stdlib=libstdc++ for gcc 4.2.1 includes ..
set( LIBCPLUSPLUS "-stdlib=libc++" )

if( BUILD_SHARED_LIBS )
#    add_definitions( -force_load ) # load all symbols static libraries even if hidden
    add_definitions( -fvisibility=hidden -fvisibility-inlines-hidden )
endif()

#### SRC && CUSTOM DEFINITIONS
# Config file
#configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/config.h )

# SOURCES
add_subdirectory( src )

# TOOLS
if( BUILD_TOOLS )
    add_subdirectory( tools )
endif( BUILD_TOOLS )

# TESTS
if( BUILD_TESTS )
    # GTest
    set( GTEST_DIR ${3RDPARTY_DIR}/gtest-1.6.0 )

    enable_testing( true )
    # Configure GTest
    if( NOT DEFINED GTEST_DIR )
        message( FATAL_ERROR "Set GTEST_DIR to the path of your gtest source folder" )
    endif()

    if( CMAKE_USE_PTHREADS_INIT )  # The pthreads library is available and allowed.
        add_definitions(-DGTEST_HAS_PTHREAD=1 -DGTEST_USE_OWN_TR1_TUPLE=1 )
    else()
        add_definitions( -DGTEST_HAS_PTHREAD=0 -DGTEST_USE_OWN_TR1_TUPLE=1 )
    endif()

    set( GTEST_INCLUDE_DIRECTORIES
        ${GTEST_DIR}/include
        ${GTEST_DIR}
        ${GTEST_DIR}/src
    )
    include_directories( ${GTEST_INCLUDE_DIRECTORIES} )
    add_library( gtest STATIC ${GTEST_DIR}/src/gtest-all.cc ${GTEST_DIR}/src/gtest_main.cc )
    add_subdirectory( tests )

endif( BUILD_TESTS )

# DOC
if( BUILD_DOC )
    # add a target to generate API documentation with Doxygen
    find_package( Doxygen )
    if( DOXYGEN_FOUND )
        configure_file( ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in ${PROJECT_SOURCE_DIR}/doc/Doxyfile @ONLY )
        add_custom_target( doc
            ${DOXYGEN_EXECUTABLE} ${PROJECT_SOURCE_DIR}/doc/Doxyfile
            COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
        file( GLOB doc "${CMAKE_CURRENT_BINARY_DIR}/doc/*.*" )
        install( FILES ${doc} DESTINATION ${DOC_INSTALL_DIR} )
    else()
        message( SEND_ERROR "Doxygen not found" )
    endif( DOXYGEN_FOUND )
endif( BUILD_DOC )



