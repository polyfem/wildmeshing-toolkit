
################################################################################
# Triwild
################################################################################

include(gmp)
include(fenvelope)

add_library(wmtk_triwild
	Smooth.cpp
	TriWild.cpp
	TriWild.h
)
add_library(wmtk::triwild ALIAS wmtk_triwild)

target_include_directories(wmtk_triwild PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/)

target_link_libraries(wmtk_triwild PUBLIC
	wmtk::toolkit
	wmtk::data
	wmtk::sec_lib
	wmtk::remeshing_lib
	igl::predicates
	gmp::gmp
)

add_executable(triwild main.cpp)
target_compile_options(triwild PRIVATE "-rdynamic")
target_link_libraries(triwild PUBLIC wmtk::triwild)
wmtk_copy_dll(triwild)
