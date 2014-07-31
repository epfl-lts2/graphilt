function( boost_python_module NAME )
    set( DEP_LIBS
        ${Boost_PYTHON_LIBRARY}
        ${PYTHON_LIBRARIES}
        ${PROJECT_NAME}
    )
    add_library( ${NAME} SHARED
      ${ARGN}
    )
    set_target_properties( ${NAME}
      PROPERTIES
      OUTPUT_NAME ${NAME}
      COMPILE_FLAGS ""
      LINK_FLAGS -dynamic
      PREFIX ""
    )
    if( WIN32 )
      set_target_properties( ${NAME} PROPERTIES SUFFIX ".pyd" )
    elseif( APPLE OR ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
      # on mac osx, python cannot import libraries with .dylib extension
      set_target_properties( ${NAME} PROPERTIES SUFFIX ".so" )
    endif()
    target_link_libraries( ${NAME}
      ${DEP_LIBS}
    )

    # Installation of wrapped library and python files
    install( TARGETS ${NAME}
        LIBRARY DESTINATION ${PYTHON_LIB_INSTALL_DIR}
    )
endfunction()
