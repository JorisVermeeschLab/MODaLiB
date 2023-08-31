
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF) # this ensures c++14 is used, otherwise gnu++14 will be used

include_directories(
        C://Users/antoi/git/MODaLiB/step1-bam-processing/dependencies/htslib-1.18
)

link_directories(
        C://Users/antoi/git/MODaLiB/step1-bam-processing/dependencies/htslib-1.18
)

add_executable(BAMConvertor main.cpp)

target_link_libraries(BAMConvertor hts)

