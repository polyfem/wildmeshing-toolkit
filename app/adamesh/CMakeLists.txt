
file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS
"src/*.cpp"
"src/*.hpp"
)

add_library(wmtk_adamesh "${LIB_SOURCES}")

add_library(wmtk::adamesh ALIAS wmtk_adamesh)

target_include_directories(wmtk_adamesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/)

target_link_libraries(wmtk_adamesh PUBLIC
	wmtk::toolkit
	wmtk::data
	igl::predicates
	FastEnvelope
	wmtk::tetwild
)


add_subdirectory(tests)