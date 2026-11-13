
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF) # this ensures c++14 is used, otherwise gnu++14 will be used

include_directories(
        dependencies/htslib-1.18
)

link_directories(
        dependencies/htslib-1.18
)

add_executable(BAMConvertor main.cpp
        params.h
        methylation/mm_ml.h
        fragment.h
        read_info.h
        reference.h
        serialization.h
        msi/vntr.h
        queue_map.h
        base.h
        pairing/single_read.h
        pairing/paired_reads.h
        strand.h
        end_motifs.h
        strand.cpp
        methylation/context.h
        methylation/context.cpp end_motifs.cpp read_info.cpp utils.h methylation/utils.h methylation/utils.cpp
        msi/vntr.cpp)

target_link_libraries(BAMConvertor hts)

