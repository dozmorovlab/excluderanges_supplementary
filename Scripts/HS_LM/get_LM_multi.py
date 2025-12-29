'''
get_LM

Brydon Wall
'''

import os
import argparse
import subprocess
import urllib.request
import gzip
import shutil
import pandas as pd
import numpy as np
from datetime import datetime

def main(kmer, bridge, filter, cutoff):

    print(f'{datetime.now()}: Downloading required files')

    # Get fai index 
    url = 'https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai'
    fai_out_path = 'hg38.fai'
    if not os.path.exists(fai_out_path):
        urllib.request.urlretrieve(url, fai_out_path)

    # Get multi-track mappability
    url = f'https://bismap.hoffmanlab.org/raw/hg38/k{kmer}.umap.bedgraph.gz'
    gzip_mapp_path = f'k{kmer}.umap.bed.gz'
    if not os.path.exists(gzip_mapp_path):
        urllib.request.urlretrieve(url, gzip_mapp_path)

    # Unzip mappability
    print(f'{datetime.now()}: Unzipping mappability')
    with gzip.open(gzip_mapp_path, 'rb') as f_in:
        with open(f'k{kmer}.umap.bed', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    # Remove the gz file after extraction
    os.remove(gzip_mapp_path)

    print(f'{datetime.now()}: Calling regions')
    in_bed = pd.read_csv(
        f'k{kmer}.umap.bed',
        sep='\t',
        skiprows=1,
        names=['chr', 'start', 'end', 'score']
    )

    #cutoff = np.quantile(in_bed['score'], quantile)
    print(f'cutoff: {cutoff}')
    in_bed = in_bed[in_bed['score'] > cutoff]

    # Correct for read length
    #in_bed['end'] += (kmer - 1)
    
    in_bed = in_bed.drop(columns=['score'])
    in_bed.to_csv(f'k{kmer}.called.umap.bed', sep='\t', header=None, index=False)

    # Merge regions
    print(f'{datetime.now()}: Bridging regions')
    merge_command_1 = (
        f'bedtools merge -i k{kmer}.called.umap.bed -d {bridge} > k{kmer}.bridged_1.umap.bed'
    )
    print(merge_command_1)
    result = subprocess.run(merge_command_1, shell=True)
    if result.returncode != 0:
        raise RuntimeError('Error running bedtools merge')

    # Get complement
    print(f'{datetime.now()}: Getting complement')
    comp_command = (
        f'bedtools complement -i k{kmer}.bridged_1.umap.bed -g {fai_out_path} -L > k{kmer}.comp.umap.bed'
    )
    print(comp_command)
    result = subprocess.run(comp_command, shell=True)
    if result.returncode != 0:
        raise RuntimeError('Error running bedtools complement')
    
    # Merge regions
    print(f'{datetime.now()}: Bridging regions 2')
    merge_command_2 = (
        f'bedtools merge -i k{kmer}.comp.umap.bed -d {bridge} > k{kmer}.bridged_2.umap.bed'
    )
    print(merge_command_2)
    result = subprocess.run(merge_command_2, shell=True)
    if result.returncode != 0:
        raise RuntimeError('Error running bedtools merge')
    
    # Filter
    print(f'{datetime.now()}: Filtering Regions')
    in_bed = pd.read_csv(
        f'k{kmer}.bridged_2.umap.bed',
        sep='\t',
        skiprows=0,
        names=['chr', 'start', 'end']
    )
    in_length = in_bed['end'] - in_bed['start']
    in_bed = in_bed[in_length >= filter]
    in_bed.to_csv(f'k{kmer}_b{bridge}_f{filter}.multi.LM.bed', sep='\t', header=None, index=False)

    # Remove intermediate files
    #os.remove('temp.bed')
    #os.remove(f'k{kmer}.umap.bridged.bed')
    #os.remove(fai_out_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fetch and process UMAP single-read low mappability BED.')
    parser.add_argument(
        '-k',
        '--kmer',
        type=int,
        default=100,
        help='k-mer to fetch'
    )
    parser.add_argument(
        '-b',
        '--bridge',
        type=int,
        default=1_000,
        help='Bridge for merging regions'
    )
    parser.add_argument(
        '-f',
        '--filter',
        type=int,
        default=1_000,
        help='Filter regions with length smaller that this'
    )
    parser.add_argument(
        '-c',
        '--cutoff',
        type=float,
        default=0.990,
        help='Cutoff value'
    )
    args = parser.parse_args()
    
    main(args.kmer, args.bridge, args.filter, args.quantile)
