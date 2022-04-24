include(catch2)
FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/contrib)
include(Catch)

file(GLOB TEST_SOURCES CONFIGURE_DEPENDS "*.cpp")
add_executable(adamesh_tests ${TEST_SOURCES})
target_link_libraries(adamesh_tests PUBLIC
    wmtk::adamesh
    Catch2::Catch2
)
wmtk_copy_dll(adamesh_tests)

set_target_properties(adamesh_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")

catch_discover_tests(adamesh_tests)
