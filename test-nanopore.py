import os


ROOT = os.path.dirname(os.path.abspath(__file__))
EXE_PATH = os.path.join(ROOT, 'step1-bam-processing', 'BAMConvertor')
REF_GENOME_PATH = os.path.join(
    ROOT, '../../../methylation_data/human_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa'
)

#filepath = '/lustre1/project/stg_00019/research/qfu/Project/Nanopore-cfDNA/data/bam/PAK67493_barcode01_sorted.bam'
#cmd = f'{EXE_PATH} {REF_GENOME_PATH} {filepath} -o nanopore_out --single_read --useMMML'
#os.system(cmd)

filepath = '/lustre1/project/stg_00019/research/qfu/Project/Nanopore-cfDNA/data/bam/PAK67493_barcode01_dorado_sorted.bam'
cmd = f'{EXE_PATH} {REF_GENOME_PATH} {filepath} -o nanopore_out2 --single_read --useMMML --mlCutoff 1'
os.system(cmd)

print('Finished.')
