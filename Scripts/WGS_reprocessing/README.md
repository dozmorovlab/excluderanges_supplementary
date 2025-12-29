# T2T variants within hg38 common regions

- `02_annotate_reannotate.sh` - commands to reannotate T2T-called SNPs focusing on regions common between T2T anf hg38 genome assemblies

- `hub_3671779_hgUniqueHg38.bed` - "CHM13 unique in comparison to GRCh38/hg38 and GRCh37/hg19", downloaded from https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=3510526275_j07KAhEyEZGgOM8wN9YlvhiFevai&db=hub_3671779_hs1&c=chr9&g=hub_3671779_hgUnique
scp hub_3671779_hgUniqueHg38.bed mdozmorov@athena.hprc.vcu.edu:/lustre/home/juicer/Blacklist/WGS_reprocessing_Mikhail/

<!--T2T variants-->
ln -s /lustre/home/juicer/Blacklist/WGS_reprocessing/T2T_CHM13/variant_calling/haplotypecaller/S1/S1.haplotypecaller.filtered.vcf.gz /lustre/home/juicer/Blacklist/WGS_reprocessing_Mikhail/T2T_CHM13/variant_calling/haplotypecaller/S1/S1.haplotypecaller.filtered.vcf.gz

ln -s /lustre/home/juicer/Blacklist/WGS_reprocessing/annotated /lustre/home/juicer/Blacklist/WGS_reprocessing_Mikhail/annotated

- `02_annotate_reannotate.sh` - reannotate T2T excluding CHM13 unique regions.
  - Input: `T2T_CHM13.annotated.vcf.gz`, `hub_3671779_hgUniqueHg38.bed.gz`
  - Output: `reannotated/ids/T2T_CHM13/autosomes.indels.ids.txt`, `autosomes.snps.ids.txt`

<!--Compare ID counts-->
wc -l /lustre/home/juicer/Blacklist/WGS_reprocessing/annotated/ids/T2T_CHM13/autosomes.indels.ids.txt
wc -l /lustre/home/juicer/Blacklist/WGS_reprocessing/annotated/ids/T2T_CHM13/autosomes.snps.ids.txt
wc -l /lustre/home/juicer/Blacklist/WGS_reprocessing_Mikhail/reannotated/ids/T2T_CHM13/autosomes.indels.ids.txt
wc -l /lustre/home/juicer/Blacklist/WGS_reprocessing_Mikhail/reannotated/ids/T2T_CHM13/autosomes.snps.ids.txt
701416 /lustre/home/juicer/Blacklist/WGS_reprocessing/annotated/ids/T2T_CHM13/autosomes.indels.ids.txt
3325440 /lustre/home/juicer/Blacklist/WGS_reprocessing/annotated/ids/T2T_CHM13/autosomes.snps.ids.txt
675580 /lustre/home/juicer/Blacklist/WGS_reprocessing_Mikhail/reannotated/ids/T2T_CHM13/autosomes.indels.ids.txt
3132361 /lustre/home/juicer/Blacklist/WGS_reprocessing_Mikhail/reannotated/ids/T2T_CHM13/autosomes.snps.ids.txt

<!--GIAB IDs-->
/lustre/home/juicer/Blacklist/WGS_reprocessing/annotated/ids/HG001_GRCh38_1_22_v4.2.1_benchmark/autosomes.indels.ids.txt
/lustre/home/juicer/Blacklist/WGS_reprocessing/annotated/ids/HG001_GRCh38_1_22_v4.2.1_benchmark/autosomes.snps.ids.txt

<!--Make figures-->
scp mdozmorov@athena.hprc.vcu.edu:/lustre/home/juicer/Blacklist/WGS_reprocessing/03_figures.R .
<!--Adjust paths and rename-->
mv 03_figures.R 03_figures_reannotated.R
scp 03_figures_reannotated.R mdozmorov@athena.hprc.vcu.edu:/lustre/home/juicer/Blacklist/WGS_reprocessing_Mikhail/

<!--New counts-->
Rscript -e "source('03_figures_reannotated.R')"
Output in `/lustre/home/juicer/Blacklist/WGS_reprocessing_Mikhail/plots`

scp -r mdozmorov@athena.hprc.vcu.edu:/lustre/home/juicer/Blacklist/WGS_reprocessing_Mikhail/plots .

<!--Original counts-->
scp -r mdozmorov@athena.hprc.vcu.edu:/lustre/home/juicer/Blacklist/WGS_reprocessing/plots plots_original



