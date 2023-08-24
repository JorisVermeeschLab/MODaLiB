cmake_minimum_required(VERSION 3.25)
project(MODAlib VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF) # this ensures c++20 is used, otherwise gnu++20 will be used

###############################################################################
# Fetch  libraries:
###############################################################################

Include(FetchContent)

# Note: At the moment, these libraries are used in tests. They are not yet
#       available to setup.py for the main extension building

# Catch2 is used for writing unit tests.
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.4.0
)
FetchContent_MakeAvailable(Catch2)

# fmt is a powerful formatting library that makes formatting and writing text easy and fast.
FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 10.1.0
)
FetchContent_MakeAvailable(fmt)


###############################################################################
# Include subdirectories
###############################################################################

add_subdirectory(step1-bam-processing)
add_subdirectory(step2-fragment-level-analysis)
