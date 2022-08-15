# HiCSampler

This script is designed to quickly randomly sample reads in Hi-C data, to create maps with lower sequencing depth. This can be useful to make two maps have approximately the same sequencing depth, or to evaluate the complexity of the libraries and determine how adding more reads will affect the data. HiCSampler works on binned data by assigning a probability of keeping each read contributing to the total binned signal depending on the desired subsampling ratio. 

HiCSampler takes Hi-C files in [juicer](https://github.com/aidenlab/juicer) format (.hic), using [straw](https://github.com/aidenlab/straw) to quickly dump the reads. It outputs a .hic file of the subsampled reads, and therefore also requires [juicer_tools](https://github.com/aidenlab/juicertools). 

## Dependencies:
<br>python3<br>numpy<br>argparse<br>hicstraw<br>pandas<br>glob<br>juicer tools<br>

## Usage: 

python subsampling.py [-h] -hic FILE -res LIST -ratio FLOAT -o OUTPUT -juicer JUICER_TOOLS [-sizes chrSizesFile] [-cpu INT]

-h, --help show this help message and exit

-hic, --hicfile FILE Input .hic file

-res LIST Resolution(s) to process the Hi-C file, can be comma separated list. Will perform subsampling at the fines scale, but create the coarser bins when creating the resultant .hic file.

-ratio FLOAT A decimal of value <1 specifying the fraction of reads you wish to keep

-o OUTPUT A directory name for the output files

-juicer JUICER_TOOLS Path to juicer tools

-sizes chrSizesFile Specify the path to a chromosome size file in 2 column format with the chromosome name and size on each row.

-cpu INT An integer specify the number of Gigs of RAM to use. Default is 4GB.

## Example:

python subsampling.py -i InputFile.hic -j juicer_tools_1.14.08.jar -res 1000,5000,10000,25000,50000,100000,250000,500000,1000000 -ratio 0.5 -o Subsample50percent -sizes hg19.sizes

The above command will create a Hi-C map with approximately 50 percent of the reads chosen radomly. It will do this calculation at 1 kb, and write into the specified output directory in juicers short score format. It will create the .hic format using all resolutions listed. 

