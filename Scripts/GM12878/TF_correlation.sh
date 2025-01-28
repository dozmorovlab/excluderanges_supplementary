
#!/bin/bash
#usage: plotCorrelation [-h] --corData FILE --corMethod {spearman,pearson} --whatToPlot {heatmap,scatterplot}
#                       [--plotFile FILE] [--skipZeros] [--labels sample1 sample2 [sample1 sample2 ...]]
#                       [--plotTitle PLOTTITLE] [--plotFileFormat FILETYPE] [--removeOutliers] [--version]
#                       [--outFileCorMatrix FILE] [--plotHeight PLOTHEIGHT] [--plotWidth PLOTWIDTH] [--zMin ZMIN]
#                       [--zMax ZMAX] [--colorMap] [--plotNumbers] [--xRange XRANGE XRANGE]
#                       [--yRange YRANGE YRANGE] [--log1p]

conda activate deeptools

out_dir="./plots"
mkdir -p "${out_dir}"

# Specify the directory to loop through
directory="./"

# Loop through all .npz files in the directory
for in_file in "$directory"/*.npz; do

    echo "Processing: ${in_file}"

    out_prefix=$(basename "$in_file" .npz)

    plotCorrelation \
        --corData "${in_file}" \
        --corMethod pearson \
        --whatToPlot heatmap \
        --plotFile "${out_dir}/${out_prefix}.pearson.heatmap.svg" \
        --outFileCorMatrix "${out_dir}/${out_prefix}.pearson.matrix.tsv"

    plotCorrelation \
        --corData "${in_file}" \
        --corMethod spearman \
        --whatToPlot heatmap \
        --plotFile "${out_dir}/${out_prefix}.spearman.heatmap.svg" \
        --outFileCorMatrix "${out_dir}/${out_prefix}.spearman.matrix.tsv"

    # No outliers --------------------------------------------------

    plotCorrelation \
        --corData "${in_file}" \
        --corMethod pearson \
        --whatToPlot heatmap \
        --plotFile "${out_dir}/${out_prefix}.no_outliers.pearson.heatmap.svg" \
        --outFileCorMatrix "${out_dir}/${out_prefix}.no_outliers.pearson.matrix.tsv" \
        --removeOutliers

    plotCorrelation \
        --corData "${in_file}" \
        --corMethod spearman \
        --whatToPlot heatmap \
        --plotFile "${out_dir}/${out_prefix}.no_outliers.spearman.heatmap.svg" \
        --outFileCorMatrix "${out_dir}/${out_prefix}.no_outliers.spearman.matrix.tsv" \
        --removeOutliers

done
