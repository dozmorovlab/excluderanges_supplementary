# Exclusion sets from realignment to the hg38 full plus sponge and the T2T-CHM13 assemblies

- `hg38_Blacklist_generated.bed` - excludable regions generated from 250 BAM files aligned to the hg38 autosomal assembly
- `hg38_full_sponge_Blacklist_generated.bed` - excludable regions from files aligned to the hg38 full plus sponge assembly
- `t2t_Blacklist_generated.bed` - excludable regions from files aligned to the T2T-CHM13 genome assembly
- `hg38_from_t2t_Blacklist_generated.bed` - hg38-liftOver `t2t_Blacklist_generated.bed` list using `chm13v2-grch38.chain`
- `hg38_from_t2t_Blacklist_unmapped.bed` - regions that cannot be liftOver'ed
- `01_Comparison.Rmd` - comparison of exclusion set characteristics

## Blacklist installation

git clone --recurse-submodules https://github.com/Boyle-Lab/Blacklist.git
cd Blacklist/bamtools/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$(cd ..; pwd)/install ..
make
make install
cd ../..

./lustre/home/mdozmorov/Work/Blacklist/Blacklist

## Mappability generation

https://github.com/hoffmangroup/newmap

Newmap installed as `conda create -n newmap bioconda::newmap -y`

**Generate** Mappability for hg38_full_sponge: /lustre/home/juicer/Blacklist/uint8_hg38_full_sponge (`01_newmap_index.sh`, `03_newmap_search.sh`). 
  - Input: Genome hg38 with contigs and sponge: /lustre/home/juicer/ExtData/UCSC/hg38_full_sponge/hg38_all.fa
  - Output: `uint8_hg38_full_sponge`

**Download** Mappability for T2T: https://bismap.hoffmanlab.org/raw/t2tv2/umap/uints/, (`download_t2t_umap_uint8.sh`)
  - Output: `uint8_t2t`

## Realignment

FASTQs: /lustre/home/juicer/Blacklist/fastqs

- `01_Additional_file_4_check.Rmd` - check if paired end experiments are paired. 
  - Input: `Additional_file_4.csv`, from `Human` sheet in `'/Users/mdozmorov/Documents/Work/GitHub/excluderanges.dev/manuscript/genome_biology_R1/Additional Files/Additional_file_4.xlsx'`
  - Output: `Additional_file_4_paired.csv` with paired-end data, paired FASTQs

- `/lustre/home/juicer/Blacklist/submit01_realignment_paired.sh` - realignment of paired FASTQs 
  - Input: `fastqs`
  - Output: `bams_hg38_full_sponge`, `bams_t2t`

## Blacklist generation

- `submit02_Blacklist_on_realigned.sh` - hg38 or t2t setting links corresponding mappability and bam files to the 'mappability' and 'input' folders in `/lustre/home/mdozmorov/Work/Blacklist` and runs array jobs 
  - Input: 'mappability' and 'input' folders
  - Output: files like `hg38_full_sponge_chr1.bed`, `t2t_chr1.bed`

### hg38

/Users/mdozmorov/Documents/Work/GitHub/excluderanges.dev/Brydon/biological_characterization/ex_data/250.bed
wc -l 250.bed
    1273 250.bed
### hg38_full_sponge realigned

scp mdozmorov@athena.hprc.vcu.edu:/lustre/home/mdozmorov/Work/Blacklist/hg38_full_sponge_chr*.bed .
<!--Combine and sort chromosomes-->
cat hg38_full_sponge_chr*.bed | sort -k1,1 -k2,2n - > hg38_full_sponge_Blacklist_generated.bed
<!--Region counts-->
wc -l hg38_full_sponge_Blacklist_generated.bed
    4029 hg38_full_sponge_Blacklist_generated.bed

### T2T realigned

scp mdozmorov@athena.hprc.vcu.edu:/lustre/home/mdozmorov/Work/Blacklist/t2t_chr*.bed results/
<!--Combine and sort chromosomes-->
cat t2t_chr*.bed | sort -k1,1 -k2,2n - > t2t_Blacklist_generated.bed
<!--Install liftOver-->
conda install bioconda::ucsc-liftover
<!--GRCh38/hg38 <- T2T-CHM13v2.0: chm13v2-grch38.chain, from https://github.com/marbl/CHM13-->
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain
<!--Replace "High Signal Region" to "High_Signal_Region", "Low Mappability" to "Low_Mappability"-->
<!--Convert t2t_Blacklist_generated.bed to hg38-->
liftOver t2t_Blacklist_generated.bed chm13v2-grch38.chain hg38_from_t2t_Blacklist_generated.bed hg38_from_t2t_Blacklist_unmapped.bed
<!--Region counts-->
wc -l t2t_Blacklist_generated.bed 
    3631 t2t_Blacklist_generated.bed
wc -l hg38_from_t2t_Blacklist_generated.bed 
    3192 hg38_from_t2t_Blacklist_generated.bed
wc -l hg38_from_t2t_Blacklist_unmapped.bed
     878 hg38_from_t2t_Blacklist_unmapped.bed
