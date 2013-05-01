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
