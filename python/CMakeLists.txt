
if( USE_DEFAULT_PYTHON )
    # Default python
    find_package( PythonLibs 2.7 REQUIRED )
    set( PYTHON_CMD python -c )
else()

    # Homebrew installed python
    set( PYTHON_CMD /usr/local/bin/python -c )
    set( PYTHON_INCLUDE_DIRS /usr/local/Frameworks/Python.framework/Headers )
    set( PYTHON_LIBRARIES /usr/local/Frameworks/Python.framework/Versions/Current/lib/libpython2.7.dylib )

endif()
include_directories( ${PYTHON_INCLUDE_DIRS} )
include(boost-python) # Special function to deal with boost python bindings

# Get Pyton site-packages folder
execute_process( COMMAND ${PYTHON_CMD} "from distutils.sysconfig import get_python_lib; print get_python_lib()"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE )

set( PYTHON_PROJECT_NAME ${PROJECT_NAME} )
set( PYTHON_INSTALL_PATH ${PYTHON_SITE_PACKAGES}/${PYTHON_PROJECT_NAME}/ )
set( PYTHON_LIB_INSTALL_DIR ${PYTHON_INSTALL_PATH} CACHE PATH "Installation prefix for exported python libraries" FORCE )
set( PYTHON_SRCS_INSTALL_DIR ${PYTHON_INSTALL_PATH} CACHE PATH "Installation prefix for exported python sources" FORCE )

add_subdirectory( wrapper )
add_subdirectory( graphilt )
