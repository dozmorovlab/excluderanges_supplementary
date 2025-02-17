# HS_LM

A directory of scripts used to generate the High Signal, Low Mappability, and combined HS_LM exclusion sets.

[36_no_dup_no_res.txt](36_no_dup_no_res.txt), [101_accessions.txt](101_accessions.txt): ENCODE accession numbers for 36bp BAM files used for ChIP-seq controls. Accessions were cross-referenced with their constituent FASTQ files to remove redundancy / duplicate FASTQs (no_dup) and ensure non-restricted base calls (no_res).

[BedCombine_example.sh](BedCombine_example.sh): An example script using [BedCombine.py](BedCombine.py).

[BedCombine.py](BedCombine.py): A python script that takes a High Signal (HS) and Low Mappability (LM) BED file and combines them with annotations indicating HS, LM, or both.

[get_LM_multi.py](get_LM_multi.py): A python script that creates the Low Mappability BED, downloading and using the references from the [Hoffman Lab](https://bismap.hoffmanlab.org/)

[MACS3_call.sh](MACS3_call.sh): A script to call High Signal regions given control signal data.

[realign_36.sh](realign_36.sh), [realign_101.sh](realign_101.sh), [realign_36_template.sh](realign_36_template.sh), [realign_101_template.sh](realign_101_template.sh): Scripts and template scripts to realign BAM files in parellel on a slurm HPRC cluster.
