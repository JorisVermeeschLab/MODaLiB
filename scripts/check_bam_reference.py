import pysam
import csv
import sys



not_found = set()
for metadata_filename in ["../config/pregnant_controls.csv"]:
    with open(metadata_filename, newline="") as f:
        reader = csv.DictReader(f, delimiter=",")
        rows = list(reader)
        filenames = [row["bam"] for row in rows]
        sample_ids = [row["sample_id"] for row in rows]
        references = [row["reference"] for row in rows]

    for sample_id, filename, ref in zip(sample_ids, filenames, references):
        filepath = "/lustre1/project/stg_00019/methylation_data/bam_files/inhouse_hg38/" + filename
        reference = ""
        print(filepath)
        with pysam.AlignmentFile(filepath, "rb") as bam:
            header = bam.header.to_dict()
            for pg in header.get("PG", []):
                if "CL" in pg:
                    if "bwameth" in pg["CL"]:
                        reference = pg["CL"]
        print(ref)
        print(reference)
        print(sample_id)
        if ref == "ref1":
            assert "/user/leuven/353/vsc35320/qfu-research/Project/methylation_data/resource/reference/hg38_lambda_puc19/GCA_000001405.15_GRCh38_no_alt_analysis_set_lambda_puc19_combined.fa" in reference or "/staging/leuven/stg_00124/methylation_data/resource/reference/hg38_lambda_puc19/GCA_000001405.15_GRCh38_no_alt_analysis_set_lambda_puc19_combined.fa" in reference or "/staging/leuven/stg_00019/methylation_data/hg38_lambda_puc19/GCA_000001405.15_GRCh38_no_alt_analysis_set_lambda_puc19_combined.fa" in reference
        else:
            assert "/staging/leuven/stg_00019/dhanya_research/human_genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" in reference
        #print(reference)

print(not_found)