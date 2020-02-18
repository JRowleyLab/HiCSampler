import numpy as np
import argparse
import gzip
import os
from straw import straw
from itertools import repeat
import time

def kdiag(mat_dim, offset):
    '''
    The "kdiag" function is used to extract the kth diagonal co-ordinates of a square matrix 
    having the dimensions of "mat_dim" and extracting the "offset" diagonal.
    '''
    row_idx, col_idx = np.diag_indices(mat_dim)
    if offset < 0:
        return row_idx[-offset:],col_idx[:offset]
    elif offset > 0:
        return row_idx[:-offset],col_idx[offset:]
    else:
        return row_idx,col_idx
'''
"randpsn" is a function used to create a poisson distribution, centered around the varaible "read".
'''
randpsn = lambda read:np.random.poisson(int(round(read)),1)[0]

def progressbar(total, j):
    if j == 0:
        print("#0", sep="", end="")
    percent = round(((j+1)/total)*100)
    if progressbar.check != percent:    
        if 0 <= percent <= 10:
            print(f"\b#{percent}", end="", sep="")
        elif 10 < percent <= 100:
            print(f"\b\b#{percent}", end="", sep="")
        progressbar.check += 1

class HiCPSN:
    '''
    "HiCPSN" is used to create the HiC map in NumPy matrix form.
    The size of the matrix is decided by the "chrsize" file inputted.
    The matrix is intialised by the "straw_result".
    Each cell in the matrix is sclaed by the given "ratio".
    '''
    def __init__(self, straw_result, chrom, ratio, res, chrsize):
        self.rowidxs = np.array(straw_result[0][:])
        self.colidxs = np.array(straw_result[1][:])
        self.scores = np.array(straw_result[2][:])
        self.nrows = int(chrsize/res)+1
        self.ncols = int(chrsize/res)+1
        if ratio > 1:
            self.rowidxs = self.rowidxs/res
            self.colidxs = self.colidxs/res
            print("The size of the matrix is: ", self.nrows, self.ncols)
            self.observed = np.zeros((self.nrows, self.ncols))
            print("Intializing HiC Matrix")
            start = time.time()
            progressbar.check = 0
            length = len(self.scores)
            for ele in range(length):
                progressbar(length, ele)
                self.observed[int(self.rowidxs[ele])][int(self.colidxs[ele])] = self.scores[ele]
            stop = time.time()
            print("\nHiC Matrix Intialized for chromosome: ", chrom)
            print(f"Time taken (min): {(stop-start)/60}")
            print("Scaling the HiC matrix with the value of ratio: ", ratio)
            start = time.time()
            self.observed = self.observed*ratio
            stop = time.time()
            print(f"Time taken for scaling the matrix with ratio {ratio} is (mins): {(stop-start)/60}")
            del ele; del straw_result;

class HiCSampler:
    '''
    "HiCSampler" is used for randomly assiging the reads from the poisson distribution in each cell in the matrix.
    As the matrix will be symmetric, the calculations are done only on the upper triangular half of the diagonal and 
    the results are copied to the lower triangular half. The diagonals with different offests are extracted using the 
    "kdiag" function. Random assignment of reads from the poisson distribution to these diagonal elements is performed 
    using the "np_randnorm" function. To remove the effects of the 1 added to the matrix, the matrix is subtracted by the
    value of the "ratio".
    '''
    def __init__(self, straw_result, chrom, chrsize, ratio=1.0, res=50000):
        self.hicpsn = HiCPSN(straw_result, chrom, ratio, res, chrsize)
        self.chrom = chrom
        self.res = res
        self.ratio = ratio
        if ratio > 1:
            print("Up-sampling reads based on Poisson distribution")
            start = time.time()
            progressbar.check = 0
            length = self.hicpsn.ncols
            for col in range(length):
                progressbar(length, col)
                diag_ele = np.diag(self.hicpsn.observed, col)
                expected = np.average(diag_ele)
                diag_ele = diag_ele + int(round(expected*ratio))
                tmp = np.array(list(map(randpsn, diag_ele)))
                tmp = tmp - int(round(expected*ratio))
                diag_ele = diag_ele - int(round(expected*ratio))
                tmp = np.where(tmp < 0, diag_ele, tmp)
                self.hicpsn.observed[kdiag(self.hicpsn.ncols, col)] = tmp
                #self.hicpsn.observed[kdiag(self.hicpsn.ncols, -col)] = tmp
            stop = time.time()
            print("Random Poisson distribution assigned for chromosome: ", self.chrom)
            print(f"Time taken (min): {(stop-start)/60}")
            del tmp; del diag_ele;
            
        elif 0 < ratio <=1 :
            print("Down-sampling reads randomly")
            start = time.time()
            progressbar.check = 0
            length = len(self.hicpsn.scores)
            for idx in range(length):
                progressbar(length, idx)
                rand_reads = np.random.random_sample(int(self.hicpsn.scores[idx]))
                read_count = len(np.where(rand_reads <= ratio)[0])
                self.hicpsn.scores[idx] = read_count
            stop = time.time()
            print("\nRandom reads assigned for chromosome: ", self.chrom)
            print(f"Time taken (min): {(stop-start)/60}")
            del rand_reads
    
    def readShortFormattedLine(self):
        '''
        The columns just need to have
        <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <score>
        You have the chromosome, the positions, and score. The strands can just be 0, and the fragments just need to be distinct from each other i.e. frag1=1 frag2=2.
        '''
        frag1=1; frag2=2;
        chr1=self.chrom; chr2=self.chrom
        str1=0; str2=0;
        if self.ratio > 1:
            for col in range(self.hicpsn.ncols):
                diag_cords = kdiag(self.hicpsn.ncols, col)
                for ele in range(len(diag_cords[0])):
                    pos1 = diag_cords[0][ele]*self.res; pos2 = diag_cords[1][ele]*self.res; score = self.hicpsn.observed[diag_cords[0][ele]][diag_cords[1][ele]]
                    if score < 1:
                        continue
                    yield '{0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(str1, chr1, pos1, frag1, str2, chr2, pos2, frag2, score)
        
        elif 0 < self.ratio <= 1:
            for idx in range(len(self.hicpsn.scores)):
                pos1 = self.hicpsn.rowidxs[idx]; pos2 = self.hicpsn.colidxs[idx]; score = self.hicpsn.scores[idx]
                if score < 1:
                        continue
                yield '{0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(str1, chr1, pos1, frag1, str2, chr2, pos2, frag2, score)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--hicfile", dest="hic_file", help="Input Hi-C file (-i) or short score format files (-d)", metavar="FILE")
    group.add_argument("-d", "--dir", dest="sht_scr_dir", help="Directory contatining processed files", metavar="Directory")
    parser.add_argument("-o", "--output", dest="output", help="Creates output directory of short with score format Hi-C file for each chromosome", metavar="Output directory name", required=True)
    parser.add_argument('--ratio', dest='ratio', default=1.0, help='ratio of sub sample (ratio>1.0) will upsample/boost the input Hi-C file')
    parser.add_argument('-s', '--size', help="Input chromosome size file if i/p is HiC file (-i)", metavar="FILE")
    parser.add_argument('--res', dest='res', default = 50000, metavar='int', help="Resolution to process the Hi-C file")
    parser.add_argument('-g', dest='gzip', action='store_true', help="Output is gzipped")
    args = parser.parse_args()
    
    if args.hic_file == None and args.size:
        parser.error("Enter chromosome size file -s, if input is HiC file -i")
        
    print("Entered Inputs....HiC file: {}, Short score format directory: {}, Output directory: {}, ratio: {}, Chromosome size file: {}, Resolution: {}, gzipped: {}".format(args.hic_file, args.sht_scr_dir, args.output, args.ratio, args.size, args.res, args.gzip))
    
    '''
    To parse through the chromosomes in .size files.
    '''        
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    
    if args.hic_file:
        sizels = list()
        with open (args.size, "r") as f:
            for line in f:
                sizels.append(line.split())
        
        for chrm in sizels:
            print("Inputs to straw ",'NONE', args.hic_file, chrm[0], chrm[0], 'BP', args.res)
            start = time.time()
            try:
                result = straw('NONE', args.hic_file, chrm[0], chrm[0], 'BP', int(args.res))
                print("Processing Chromosome: ", chrm[0])
            except:
                print(f"Error in processing chromosome {chrm[0]}")
            stop = time.time()
            print(f"Time taken for dumping reads (min): {(stop-start)/60}")
            try:
                del mixer
            except:
                print("No HiCSampler object intialized")
    
            mixer = HiCSampler(result, chrm[0], ratio=float(args.ratio), res=int(args.res), chrsize=int(chrm[1]))
            del result
            file = args.output+"/chr"+chrm[0]
            print("Writing short score format file")
            start = time.time()
            if args.gzip:
                with gzip.open(file, 'wt') as f:
                    for i in mixer.readShortFormattedLine():
                        f.write(i)
                        f.write('\n')
            else: 
                with open(file, 'w+') as f:
                    for i in mixer.readShortFormattedLine():
                        f.write(i)
                        f.write('\n')
            stop = time.time()
            print(f"Time taken for writing short score format file (mins): {(stop-start)/60}")
            print("Processed Chromosome: ", chrm[0])
    
    elif args.sht_scr_dir:
        for chrfile in os.listdir(args.sht_scr_dir):
            result0 = list(); result1 = list(); result2 = list(); result = list()
            print("Dumping reads")
            start = time.time()
            with open (os.path.join(args.sht_scr_dir, chrfile), "r") as f:
                for line in f:
                    line_values = line.split()
                    result0.append(int(line_values[2]))
                    result1.append(int(line_values[6]))
                    result2.append(int(float(line_values[8])))
                    chrom = line_values[1]
            chrsize = max(max(result0), max(result1))
            result.append(result0); result.append(result1); result.append(result2)
            stop = time.time()
            print("Processing Chromosome: ", chrom)
            print(f"Time taken for dumping reads (min): {(stop-start)/60}")
            try:
                del mixer
            except:
                print("No HiCSampler object intialized")
            mixer = HiCSampler(result, chrom, ratio=float(args.ratio), res=int(args.res), chrsize=chrsize)
            del result
            file = args.output+"/chr"+str(chrom)
            print("Writing short score format file")
            start = time.time()
            if args.gzip:
                with gzip.open(file, 'wt') as f:
                    for i in mixer.readShortFormattedLine():
                        f.write(i)
                        f.write('\n')
            else: 
                with open(file, 'w+') as f:
                    for i in mixer.readShortFormattedLine():
                        f.write(i)
                        f.write('\n')
            stop = time.time()
            print(f"Time taken for writing short score format file (mins): {(stop-start)/60}")
            print("Processed Chromosome: ", chrom)

