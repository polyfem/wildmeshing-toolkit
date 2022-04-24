
################################################################################
# FloatTetwild
################################################################################

include(gmp)
include(fenvelope)

file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS "src/*.h" "src/*.cpp" "src/*.hpp")

add_library(wmtk_tetwild ${LIB_SOURCES})
add_library(wmtk::tetwild ALIAS wmtk_tetwild)

target_include_directories(wmtk_tetwild PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/)

target_link_libraries(wmtk_tetwild PUBLIC
	Tracy::TracyClient
	wmtk::toolkit
	wmtk::data
	igl::predicates
	FastEnvelope
	gmp::gmp
)

add_executable(tetwild app/main.cpp)
target_compile_options(tetwild PRIVATE "-rdynamic")
target_link_libraries(tetwild PUBLIC wmtk::tetwild wmtk::sec_lib)
wmtk_copy_dll(tetwild)
