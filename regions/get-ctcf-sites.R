# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# BiocManager::install("AnnotationHub", update = FALSE) 
# BiocManager::install("GenomicRanges", update = FALSE)
# BiocManager::install("plyranges", update = FALSE)

suppressMessages(library(AnnotationHub))
ah <- AnnotationHub()

query_data <- subset(ah, preparerclass == "CTCF")

subset(query_data, species == "Homo sapiens" & genome == "hg38" & dataprovider == "JASPAR 2022")

CTCF_hg38_all <- query_data[["AH104727"]]
CTCF_hg38 <- query_data[["AH104729"]]

suppressMessages(library(plyranges))
CTCF_hg38_all <- CTCF_hg38_all %>% keepStandardChromosomes() %>% sort()
CTCF_hg38 <- CTCF_hg38 %>% keepStandardChromosomes() %>% sort()

write.table(CTCF_hg38_all %>% sort() %>% as.data.frame(), 
            file = "CTCF_hg38_all.bed",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(CTCF_hg38 %>% sort() %>% as.data.frame(), 
            file = "CTCF_hg38.bed",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
