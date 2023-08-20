
# cnpy can work both with npy and npz files. At the moment we don't use it, but include it for future reference
FetchContent_Declare(
        cnpy
        GIT_REPOSITORY https://github.com/rogersce/cnpy.git
        GIT_TAG master # the project does not have a tag
)
FetchContent_MakeAvailable(cnpy)

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${cnpy_SOURCE_DIR}
)

add_subdirectory(tests)







