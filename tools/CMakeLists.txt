include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../src )

# For all tools
set( LIBS_TO_LINK
    ${Boost_LIBRARIES}
    ${LIBCPLUSPLUS}
    ${OPENCL_LIBRARIES}
    ${PROJECT_NAME}_static
)

# Device info
add_executable( PlateformInfo plateformInfo.cpp )
target_link_libraries( PlateformInfo
    ${LIBS_TO_LINK}
)

# Demo1
add_executable( Demo1 demo1.cpp )
target_link_libraries( Demo1
    ${LIBS_TO_LINK}
)

# Downsample
add_executable( Downsample downsample.cpp )
target_link_libraries( Downsample
    ${LIBS_TO_LINK}
)
