# Sources add all sources here
set( SRCS
    main.cpp
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
)

#### TARGET
# Shared library if testing or building tools
add_library( lib${PROJECT_NAME}
    SHARED
    ${ALL_SRCS}
)

set_target_properties( lib${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_SOVERSION}
    INSTALL_NAME_DIR ${LIB_INSTALL_DIR}
    DEFINE_SYMBOL GRAPHILT_API_EXPORTS
)

target_link_libraries( lib${PROJECT_NAME}
    ${LIBS_TO_LINK}
)

# STATIC
remove_definitions( -fvisibility=hidden -fvisibility-inlines-hidden )
add_library( lib${PROJECT_NAME}_static
    STATIC
    ${ALL_SRCS}
)

target_link_libraries( lib${PROJECT_NAME}_static
    ${LIBS_TO_LINK}
)

# Executable
add_executable( ${PROJECT_NAME} 
    ${ALL_SRCS}
)
target_link_libraries( ${PROJECT_NAME}
    ${LIBS_TO_LINK}
)

