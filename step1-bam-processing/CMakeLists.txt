
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF) # this ensures c++14 is used, otherwise gnu++14 will be used

include_directories(
        dependencies/htslib-1.13/
)

link_directories(
        dependencies/htslib-1.13
)

add_executable(BAMConvertor main.cpp)

target_link_libraries(BAMConvertor libhts.dylib)




