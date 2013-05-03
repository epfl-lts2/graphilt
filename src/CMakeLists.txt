# Sources add all sources here
set( SRCS
    util/log.cpp
    io/vector-io.h
    io/matrix-io.h
    io/parser/snap.h
    io/parsers.h
)

# Public headers
set( PUBLIC_HDRS
)

# All sources, add sources groupe here
set( ALL_SRCS
    ${SRCS}
)

# Set here libraries to link to the target
set( LIBS_TO_LINK
    ${Boost_LIBRARIES}
    ${LIBCPLUSPLUS}
    ${OPENCL_LIBRARIES}
)

#### TARGET
if( BUILD_SHARED_LIBS )
    add_library( ${PROJECT_NAME}
        SHARED
        ${ALL_SRCS}
    )

    set_target_properties( ${PROJECT_NAME} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_SOVERSION}
        INSTALL_NAME_DIR ${LIB_INSTALL_DIR}
        DEFINE_SYMBOL GRAPHILT_API_EXPORTS
    )

    target_link_libraries( ${PROJECT_NAME}
        ${LIBS_TO_LINK}
    )
endif( BUILD_SHARED_LIBS )

# STATIC
remove_definitions( -fvisibility=hidden -fvisibility-inlines-hidden )
add_library( ${PROJECT_NAME}_static
    STATIC
    ${ALL_SRCS}
)

target_link_libraries( ${PROJECT_NAME}_static
    ${LIBS_TO_LINK}
)

## INSTALL LIBS AND EXECUTABLES
install( TARGETS ${PROJECT_NAME}
    FRAMEWORK DESTINATION ${FRAMEWORK_INSTALL_DIR}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${BIN_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
)

## INSTALL PUBLIC HEADERS
set( DEST ${INCLUDE_INSTALL_DIR}/${PROJECT_NAME} )
install( FILES ${PUBLIC_HDRS} DESTINATION ${DEST} )



