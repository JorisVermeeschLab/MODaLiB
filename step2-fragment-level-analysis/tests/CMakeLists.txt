

# define directories. All directory paths include a trailing / at the end.
add_definitions(-DNUCLEOSOMES_DIR="${CMAKE_SOURCE_DIR}/nucleosomes/")
add_definitions(-DSTEP_1_OUT_DIR="${CMAKE_SOURCE_DIR}/step1-bam-processing/out/")
add_definitions(-DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")


add_executable(test_step2 test_step2.cpp)
target_link_libraries(test_step2 PRIVATE modalib_obj Catch2::Catch2WithMain fmt)