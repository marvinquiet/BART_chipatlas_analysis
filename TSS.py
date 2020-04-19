import os

chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
            'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
            'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX', 'chrY']

def get_TSS_from_refFlat(refFlat_file, output_file):
    # protein coding genes from HGNC
    protein_coding_genes = []
    with open('data/hg19_TSS/gene_with_protein_names.txt', 'r') as f:
        for line in f:
            protein_coding_genes.append(line.strip())

    with open(refFlat_file, 'r') as fin:
        with open(output_file, 'w') as fout:
            for line in fin:
                res = line.strip().split('\t')
                chr_id = res[0]
                gene_symbol = res[3]
                if gene_symbol not in protein_coding_genes or chr_id not in chr_list: # filter protein coding and chromosome
                    continue

                strand = res[4]
                if strand == '+':
                    TSS = int(res[1])
                else:
                    TSS = int(res[2])
    
                TSS_minus_5kb = 0 if TSS < 5000 else TSS-5000
                TSS_plus_5kb = TSS+5000
                fout.write('%s\t%d\t%d\t%s\t%s\n' % (chr_id, TSS_minus_5kb, TSS_plus_5kb, gene_symbol, strand))

def get_uniq_gene_bed(gene_sorted_file, output_file):
    fout = open(output_file, 'w')

    with open(gene_sorted_file, 'r') as fopen:
        last_gene = None
        for line in fopen:
            res = line.strip().split('\t')
            if last_gene is None:
                last_gene = res[3]

    fout.close()


if __name__ == '__main__':
    # refFlat_file = "data/hg19_TSS/hg19_refFlat.txt"
    TSS_file = "data/hg19_TSS/uniqueTSS.hg19.bed"
    TSS_5kb_file = "data/hg19_TSS/uniqueTSS.hg19.5kb.bed"
    get_TSS_from_refFlat(TSS_file, TSS_5kb_file)


