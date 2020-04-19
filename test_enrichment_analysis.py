import os
import scipy.stats as stats
import numpy as np
from operator import itemgetter
from datetime import datetime

import find_overlap_keep_info_NOT_sep_strand_lastColMarked

# === environmental settings
GENE_LIST_DIR = "data/fc1.5_genelist"
SRX_PEAK_DIR = "data/TFs_peaks_per_srx"
GENES_PEAK_FILE = "data/hg19_TSS/uniqueTSS.hg19.5kb.sorted.bed"

GENES_FILTERED_PEAKS_DIR = "data/filtered_fc1.5_gene_bed"
os.makedirs(GENES_FILTERED_PEAKS_DIR, exist_ok=True)

FISHER_RESULT_DIR = "data/result_20200417"
os.makedirs(FISHER_RESULT_DIR, exist_ok=True)

def get_input_genes(gene_list_file):
    genes = []
    with open(gene_list_file, 'r') as fopen:
        for line in fopen:
            genes.append(line.strip())
    return genes

def get_genes_peaks(gene_list_file):
    genes = get_input_genes(gene_list_file)
    gene_filtered_peaks_file = os.path.join(GENES_FILTERED_PEAKS_DIR, 
                                            os.path.basename(gene_list_file)+'_'
                                            +os.path.basename(GENES_PEAK_FILE))
    fout = open(gene_filtered_peaks_file, 'w')

    with open(GENES_PEAK_FILE, 'r') as fopen:
        for line in fopen:
            res = line.strip().split('\t')
            if res[3] in genes:
                fout.write(line)
    fout.close()
    return gene_filtered_peaks_file

def run_intersection(srx_peak, gene_filtered_peak, gene_bed_num, target_length):
    target_overlap = find_overlap_keep_info_NOT_sep_strand_lastColMarked.main(gene_filtered_peak, srx_peak) 
    overall_overlap = find_overlap_keep_info_NOT_sep_strand_lastColMarked.main(GENES_PEAK_FILE, srx_peak)

    nontarget_overlap = overall_overlap - target_overlap
    target_nonoverlap = target_length - target_overlap
    nontarget_nonoverlap = gene_bed_num - target_length - nontarget_overlap

    oddsratio, pvalue = stats.fisher_exact(
            [[target_overlap, target_nonoverlap], 
             [nontarget_overlap, nontarget_nonoverlap]])
    return np.log10(pvalue), target_overlap, nontarget_overlap

def run_analysis(gene_list_file):
    gene_filtered_peaks_file = get_genes_peaks(gene_list_file)
    gene_bed_num = int(os.popen("wc -l {}".format(GENES_PEAK_FILE)).read().strip().split()[0])
    target_length = int(os.popen("wc -l {}".format(gene_filtered_peaks_file)).read().strip().split()[0])
    
    analysis_res = []
    for srx_peak in os.listdir(SRX_PEAK_DIR):
        fisher_test_res = run_intersection(os.path.join(SRX_PEAK_DIR, srx_peak),
                                           gene_filtered_peaks_file,
                                           gene_bed_num,
                                           target_length)
        res_tuple = (srx_peak.replace('.bed', ''),
                fisher_test_res[0],
                fisher_test_res[1],
                target_length,
                fisher_test_res[2],
                gene_bed_num-target_length)
        # print(res_tuple)
        analysis_res.append(res_tuple)
    return analysis_res

def write_analysis(gene_list_file):
    if os.path.exists(os.path.join(FISHER_RESULT_DIR, gene_list_file.replace('.txt', '.res'))):
        return
    now = datetime.now()
    print("time %s, processing %s ..." % (now.strftime("%H:%M:%S"), gene_list_file))
    analysis_res = run_analysis(os.path.join(GENE_LIST_DIR, gene_list_file))
    now = datetime.now()
    print("time %s, finish processing %s ..." % (now.strftime("%H:%M:%S"), gene_list_file))
    sorted_analysis_res = sorted(analysis_res, key=itemgetter(1))

    gene_result_file = os.path.join(FISHER_RESULT_DIR, gene_list_file.replace('.txt', '.res'))
    with open(gene_result_file, 'w') as f:
        for res in sorted_analysis_res:
            f.write('%s\t%.2f\t%d/%d\t%d/%d\n' % (res[0], res[1], res[2], res[3], res[4], res[5]))

if __name__ == '__main__':
    # for gene_list_file in os.listdir(GENE_LIST_DIR):
    #     write_analysis(gene_list_file)
    from multiprocessing import Pool
    pool = Pool(10)
    pool.map(write_analysis, os.listdir(GENE_LIST_DIR))
