#seurat4-scrublet-filter.sh
#must be run in conda environment containg R and Seurat v4
#Example command call: bash /diskmnt/Projects/Users/austins2/tools/seurat4-scrublet-no-filter-cellranger-arc2.sh /diskmnt/Projects/PDX_scRNA_analysis/matching_snRNAseq/PDAC/cellranger-arc2 /diskmnt/Projects/PDX_scRNA_analysis/matching_snRNAseq/PDAC/scrublet
#file structure for when command is run:
#incoming-path/
            # /cellranger #provide the fullpath here
            #            /Sample-1
            #                     /outs
            #            /sample-2
            #                     /outs
            #            /3rd-sample
            #                       /outs
            # /scrublet #this directory can actually be anywhere you want it to be so long as it exists. It also can be named something else. provide the full path here second
            #            /Sample-1
            #                     /Sample-1_scrublet_output_table.csv
            #            /sample-2
            #                     /sample-2_scrublet_output_table.csv
            #            /3rd-sample
            #                       /3rd-sample_scrublet_output_table.csv
#
folder=$(pwd)
cellrangerpath=$1 
scrublet=$2
for dir in $cellrangerpath/*; do
    dir=${dir%*/}
    dir="${dir##*/}"
    echo $dir
    if [ -d "$cellrangerpath"/"$dir" ]
    then
        if [ -d "$cellrangerpath"/"$dir"/outs ]
        then
            echo "Cellranger for $dir exists"
            if [ -d "$scrublet"/"$dir" ]
            then
                echo "Scrublet for $dir exists"
                Rscript /diskmnt/Projects/Users/austins2/tools/seurat4_no_removal_scrublet_sctv2_cellranger-arc2.R -i $cellrangerpath/$dir/outs -s $dir -o $dir --scrublet $scrublet/$dir
            else
                echo "Scrublet for $dir does not exist at the provided location"
            fi
        else
            echo "Cellranger outs for $dir does not exist at the provided location"
        fi
    fi
done


