import pysam
import csv
import argparse
import os

def extract_alleles(vcf_path, sample, output_file):
    vcf_catalog = {}
    if not os.path.exists(vcf_path):
        vcf_path = vcf_path + '.gz'
    assert os.path.exists(vcf_path)
    filename, filext = os.path.splitext(vcf_path)
    
    if ".vcf" in filext:
        pysam.tabix_index(vcf_path, preset="vcf", force=True)
        vcf_path = filename + ".vcf.gz"

    with pysam.VariantFile(vcf_path) as vcf_file:
        for record in vcf_file:
            desc_line = record.info
            var_id = desc_line.get("VARID")
            rep_id = desc_line.get("REPID")
            repeat = desc_line.get("RU")
            id_string = var_id if var_id else rep_id
            id_string = f"{id_string}:{repeat}"
            gt = record.samples[0].get("REPCN")
            if gt == "." or gt == "./.": 
                continue
            rep_counts = [int(x) for x in gt.split("/")]
            max_allele = max(rep_counts)
            min_allele = min(rep_counts) if len(rep_counts) > 1 else 'None'
            vcf_catalog[id_string] = (max_allele, min_allele)
    
    with open(output_file, 'w', newline='') as tsvfile:
        fieldnames = ['SAMPLE', 'ID_STRING', 'MAX_ALLELE', 'MIN_ALLELE']
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for id_string, (max_allele, min_allele) in vcf_catalog.items():
            writer.writerow({'SAMPLE': sample, 'ID_STRING': id_string, 
                           'MAX_ALLELE': max_allele, 'MIN_ALLELE': min_allele})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Extract allele counts from ExpansionHunter VCF output'
    )
    parser.add_argument('vcf_path', help='Path to VCF file')
    parser.add_argument('sample', help='Sample identifier')
    parser.add_argument('-o', '--output', 
                        help='Output TSV file (default: SAMPLE_allele_counts.tsv)')
    
    args = parser.parse_args()
    
    output_file = args.output if args.output else f'{args.sample}_allele_counts.tsv'
    extract_alleles(args.vcf_path, args.sample, output_file)