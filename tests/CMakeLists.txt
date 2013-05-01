include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

macro(append_test test_name test_srcs)
    set(TESTS_SRCS ${TESTS_SRCS} ${test_srcs} )
    set(TESTS_NAME ${TESTS_NAME} ${test_name})
endmacro(append_test)

# APPEND TESTS HERE
append_test(ExampleTest Example.cpp)

list(LENGTH TESTS_SRCS len1)
math(EXPR len2 "${len1} - 1")

foreach(val RANGE ${len2})
    list(GET TESTS_SRCS ${val} test_src)
    list(GET TESTS_NAME ${val} test_name)
    add_executable( ${test_name} ${test_src} )
    target_link_libraries( ${test_name}
        gtest
        ${PROJECT_NAME}_static
        ${Boost_LIBRARIES}
        ${LIBCPLUSPLUS}
        ${OPENCL_LIBRARIES}
    )
    get_property( CONFIG_LOCATION TARGET ${test_name} PROPERTY LOCATION )
    add_test( NAME ${test_name} COMMAND ${test_name} )

endforeach()
