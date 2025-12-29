k="100"
aligner="bowtie2"
d="1000"
q="99.90"
sample="${aligner}_${k}"
python BedCombine.py \
    -s "./results/HS/chr1/${sample}_1_${q}_chr1.bed" \
    -m "./results/LM/chr1/k${k}_LM_chr1.bed" \
    -o "./results/Both/chr1/${sample}_${d}_${q}_chr1_combined.bed" \
    -d ${d}