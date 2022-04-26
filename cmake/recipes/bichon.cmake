if(TARGET bichon)
    return()
endif()

message(STATUS "Third-party: creating target 'Bichon'")

option(PRISM_LIB_ONLY "prism lib only" ON)
option(PRISM_TESTS "prism lib only" OFF)

include(FetchContent)
FetchContent_Declare(
    bichon
    GIT_REPOSITORY https://github.com/jiangzhongshi/bichon.git
    GIT_TAG polyshell
)
FetchContent_MakeAvailable(bichon)
