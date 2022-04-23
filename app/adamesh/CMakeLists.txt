
file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS
"*.cpp"
"*.hpp"
)

add_library(wmtk_adamesh "${LIB_SOURCES}")
# target_include_directories(wmtk_adamesh PUBLIC )

target_link_libraries(wmtk_adamesh PUBLIC
	wmtk::toolkit
	wmtk::data
	igl::predicates
	FastEnvelope
	wmtk::tetwild
)

add_subdirectory(tests)