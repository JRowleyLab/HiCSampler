#!usr/bin/env python
import numpy as np
import os
import straw
import pandas as pd
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-hic", "--hicfile", dest="hic", help="input .hic file\n\n", metavar="HiC File", required=True)
parser.add_argument("-res", dest="res", help="input the resolutions to process the HiC file, can be comma spearated list\n\n", metavar="LIST", required=True)
parser.add_argument("-ratio", dest="ratio", help="input the ratio for subsampling\n\n", metavar="FLOAT", default="1.0",  required=True)
parser.add_argument("-o", dest="output", help="path to output files\n\n", metavar="OUTPUT PATH", required=True)
parser.add_argument("-juicer", dest="juicer", help="path to juicer tools\n\n", metavar="JUICER_TOOLS PATH", required=True)
parser.add_argument("-sizes", dest="sizes", help="chromosome size file with two columns corresponding to chromosome and size respectively\n\n", metavar="SIZE FILE")
parser.add_argument("-cpu", dest="cpu", help="number of Gigs of RAM to use; default is 4GB\n\n", metavar="INT", default=4)

args = parser.parse_args()

if not os.path.exists(args.output):
    os.mkdir(args.output)
    
chrom = [str(i) for i in range(1,22)]
chrom += ["X", "Y"]

inputres = args.res
myres = inputres.split(",")
myres = list(map(int, myres))
highres = min(myres)

for curchrom in chrom:
    dumped = straw.straw("NONE", args.hic, curchrom, curchrom, "BP", highres)
    rowidxs = np.array(dumped[0][:])
    colidxs = np.array(dumped[1][:])
    scores = np.array(dumped[2][:])
    length = len(scores)
    str1 = list(); frag1=list(); frag2=list(); chr_ls = list()
    
    for idx in range(length):
        rand_reads = np.random.random_sample(int(scores[idx]))
        read_count = len(np.where(rand_reads <= float(args.ratio))[0])
        scores[idx] = read_count
        str1.append(0)
        frag1.append(1)
        frag2.append(2)
        chr_ls.append(curchrom)
   
    data = {"rowidxs":rowidxs, "colidxs":colidxs, "scores":scores}
    df = pd.DataFrame(data)
    df.insert(0, "str1", str1, True)
    df.insert(1, "chr1", chr_ls, True)
    df.insert(3, "frag1", frag1, True)
    df.insert(4, "str2", str1, True)
    df.insert(5, "chr2", chr_ls, True)
    df.insert(7, "frag2", frag2, True)
    output_path = os.path.join(args.output, curchrom+".ssfv")
    df.to_csv(output_path, sep=" ", header=False, index=False)
    
read_files = glob.glob(os.path.join(args.output,"*.ssfv"))
output_path = os.path.join(args.output, "subsampled")
with open(output_path, "wb") as outfile:
    for f in read_files:
        with open(f, "rb") as infile:
            outfile.write(infile.read())
hicoutpath = output_path + ".hic"
del_files = os.path.join(args.output, "*.ssfv")
my_cmd = "rm " + del_files
os.system(my_cmd)
cpu = "-Xmx"+args.cpu+"g"
mycmd =  'java ' + cpu + ' -jar ' + args.juicer + ' pre ' + '-r ' + args.res + ' ' + output_path + ' ' + hicoutpath + ' ' + args.sizes
os.system(mycmd)
