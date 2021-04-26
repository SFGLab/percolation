pipe.sh file is a script processing bed and bigwig files for different features. pipe.py script run under it.

pipe.py  is a script to calculate mean of the peak signal, sum of the peak signal, minimal peak, maximal peak from files being the result of intersection of loop files and feature bed files.

To run pipe.sh a ./4columns folder is needed which contain feature peak values in four columns bed format <chr, start, end, peak_value> and a ./bigwig folder with files describing features in bigwig format (signal for whole genome not only peaks). The loop file (<loop_name>.bed)in bed format (chrom, start, end) is also required in the folder ./loop.

The output feature table file will be written to ./output/<loop_name>/ folder.

Required tools:
bedtools (https://bedtools.readthedocs.io/)
bigWigToBedGraph (UCSC tool, http://hgdownload.soe.ucsc.edu/admin/exe/) 
