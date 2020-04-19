import os,sys,argparse,glob
import numpy as np
import pandas as pd
from GenomeData import *
import re,bisect
from operator import itemgetter

plus = re.compile("\+")
minus = re.compile("\-")

def union_islands(islandlist):
    #start/[0] position in islandlist must be sorted
    unionlist = []
    i = 1
    currentStart,currentEnd = islandlist[0][0],islandlist[0][1]
    while i<len(islandlist):
        compareStart,compareEnd = islandlist[i][0],islandlist[i][1]
        assert currentStart <= compareStart
        if compareStart > currentEnd:
            unionlist.append([currentStart,currentEnd])
            currentStart,currentEnd = compareStart,compareEnd 
            i+=1
        else:
            currentEnd = max(currentEnd,compareEnd)
            i+=1    
        #print(i)   
    unionlist.append([currentStart,currentEnd])
    return unionlist
    
def get_overlap_info_strID(islands,regions):
    # for each island in islands[chrom], check if 1-overlap 0-non-overlap with regions[chrom]
    # islands/regions represents for compare_regions/basic_region respectively
    #island_dict = {}
    #i=0
    overlapped = {}
    for chrom in islands:
        if chrom not in regions:
            pass
        else:
            # check the regions in regions[chrom] are non-overlap and sorted
            regionlist = []
            for region in regions[chrom]:
                regionlist.extend(region)
            for num in range(len(regionlist)-1):
                assert regionlist[num]<regionlist[num+1]
                            
            # for each island, check if overlapped with regions[chrom]
            for island in islands[chrom]:
                #ID = chrom+'_'+str(island[0])+'_'+str(island[1])+'_'+str(i)
                assert island[0]<=island[1]
                s = bisect.bisect_right(regionlist,island[0])
                e = bisect.bisect_left(regionlist,island[1])
                if s==e and s%2==0:
                    pass
                else:
                    if chrom not in overlapped:
                        overlapped[chrom] = []
                    overlapped[chrom].append(island)
                
    return overlapped

def get_overlap_info_numID(islands,regions):
    # for each island in islands[chrom], check if 1-overlap 0-non-overlap with regions[chrom]
    overlapped,nonoverlapped = {},{}
    for chrom in islands:
        if chrom not in regions:
            for island in islands[chrom]:
                if chrom not in nonoverlapped:
                    nonoverlapped[chrom] = []
                nonoverlapped[chrom].append(island)
        else:
            # check the regions in regions[chrom] are non-overlap and sorted
            regionStart,regionEnd = [],[]
            leftEnd=0
            for region in regions[chrom]:
                assert leftEnd <= region[0]<=region[1]            
                regionStart.append(region[0])
                regionEnd.append(region[1])
                leftEnd=region[1]
                            
            # for each island, check if overlapped with regions[chrom]
            for island in islands[chrom]:              
                assert island[0]<=island[1]
                e = bisect.bisect_right(regionEnd,island[0])
                s = bisect.bisect_left(regionStart,island[1])
                if s > e :
                    if chrom not in overlapped:
                        overlapped[chrom] = []
                    overlapped[chrom].append(island)
                else:
                    if chrom not in nonoverlapped:
                        nonoverlapped[chrom] = []
                    nonoverlapped[chrom].append(island)
                #i+=1
    return overlapped,nonoverlapped

def read_region_from_bed(infile,species,expand=0,write_out=None,end_col=2,region_sorted=True):
    if species in species_chroms.keys():
        chroms = species_chroms[species]
    else:
        print('Species not recognized!')
        exit(1)
    regions1,regions2 = {},{}
    with open(infile,'r') as inf:
        line = inf.readline()
        while line:
            #print(line)
            if not re.match("#",line):
                #print(line)
                sline = line.strip().split('\t')
                chrom = sline[0]
                start = int(sline[1])
                end = int(sline[end_col])
                #strand = sline[5]
                #if plus.match(strand):
                if 1:
                    if chrom in chroms:
                        if chrom not in regions1:
                            regions1[chrom]=[]
                        regions1[chrom].append([max(0,start-expand),end+expand,'\t'.join(sline[3:])])
                    else:
                        pass   
                elif minus.match(strand):
                    if chrom in chroms:
                        if chrom not in regions2:
                            regions2[chrom]=[]
                        regions2[chrom].append([max(0,start-expand),end+expand,'\t'.join(sline[3:])])
                    else:
                        pass                                        
            #print(regions1,regions2);exit(0)
            line = inf.readline()
            
    if region_sorted:
        for chrom in regions1:
            regions1[chrom] = sorted(regions1[chrom],key = itemgetter(0),reverse=False)
        for chrom in regions2:
            regions2[chrom] = sorted(regions2[chrom],key = itemgetter(0),reverse=False)
            
    if write_out:
        with open(write_out,'w') as outf:
            for chrom in regions1:
                for region in regions1[chrom]:
                    outf.write('{}\t{}\t{}\t{}\n'.format(chrom,region[0],region[1],region[2]))
            for chrom in regions2:
                for region in regions2[chrom]:
                    outf.write('{}\t{}\t{}\t{}\n'.format(chrom,region[0],region[1],region[2]))
                        
    return regions1,regions2


def write_regions_bed3(regions1,regions2,outfile):
    # write the regions(with chrom,start,end info) into bed3 format
    with open(outfile,'w') as outf:
        for chrom in regions1:
            for region in regions1[chrom]:
                outf.write('{}\t{}\t{}\t{}\t1\n'.format(chrom,region[0],region[1],region[2]))
        for chrom in regions2:
            for region in regions2[chrom]:
                outf.write('{}\t{}\t{}\t{}\t0\n'.format(chrom,region[0],region[1],region[2]))  


def count_overlap_in_region1(regions1):
    cnt = 0
    for chrom in regions1:
        for region in regions1[chrom]:
            cnt += 1
    return cnt

def process_peak_file(infile2):
    expand_region = 150
    with open(infile2) as inf:
        df = pd.read_csv(inf,sep='\t',header=None,index_col=3)
    
    df.columns = ['chr','start','end','n1','n2','score','n3','n4','dis']
    #df = df[df['score']>4]
    #df = df.sort_values(by=['score'],ascending=False)
    #df = df.iloc[:10000,:]
    #print(df)
    #exit()
    df['summit'] = df['start']+df['dis']
    df['left'] = df['summit']-expand_region
    df['left'] = df['left'].clip(0)
    df['right'] = df['summit']+expand_region
    df = df[['chr','left','right']]
    
    summit_file = 'summit.bed'
    df.to_csv(summit_file,index=False,sep='\t',header=None)
    return summit_file
    #exit()
    

def main(infile1,infile2,species='hg19',expand1=0,expand2=0):
    #infile2 = process_peak_file(infile2)
    compare_regions1,compare_regions2 = read_region_from_bed(infile1,species,expand1)
    basic_regions1,basic_regions2 = read_region_from_bed(infile2,species,expand2)
    for ele in basic_regions1:
        basic_regions1[ele] = union_islands(basic_regions1[ele])
    for ele in basic_regions2:
        basic_regions2[ele] = union_islands(basic_regions2[ele])

    overlapped1 = get_overlap_info_strID(compare_regions1,basic_regions1)        
    overlapped1_cnt = count_overlap_in_region1(overlapped1)
    return overlapped1_cnt
    # write_regions_bed3(overlapped1,nonoverlapped1,outfile1) 
    # process_overlapped_file(outfile1)       
        
        
if __name__ == '__main__':
    infile1 = "data/hg19_TSS/uniqueTSS.hg19.5kb.sorted.bed"
    infile2 = "data/TFs_peaks_per_srx/SRX317574.bed"
    print(main(infile1, infile2))
