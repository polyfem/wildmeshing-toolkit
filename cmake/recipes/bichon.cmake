if(TARGET bichon)
    return()
endif()

message(STATUS "Third-party: creating target 'Bichon'")

include(FetchContent)
FetchContent_Declare(
    bichon
    GIT_REPOSITORY https://github.com/jiangzhongshi/bichon.git
    GIT_TAG polyshell
)
FetchContent_MakeAvailable(bichon)
