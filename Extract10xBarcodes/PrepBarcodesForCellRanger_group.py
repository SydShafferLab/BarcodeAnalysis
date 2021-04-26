# Run this script by entering the 10xbc environment and then runing: python /path/to/PrepBarcodesForCellRanger.py

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from fuzzywuzzy import fuzz
import matplotlib as plt
from matplotlib import pyplot
import numpy as np
import glob
import os
import subprocess
from time import sleep
import statistics
import shutil

##### Users should modulate these parametes to run this script

#path to all barcode analyssi scripts
scripts="/home/gharm/BarcodeAnalysis"

#Folder that contains all folders containing FASTQ files generated from sequencing the barcodes
Fastqfolder="/project/shafferslab/Guillaume/10X_exp1_reanalysis/20190808_10X1_BC_r1_r2_seq1_fastq/outs/fastq_path/10x"

#folder you want outputs go go into (dont make this folder, this scipt will make it)
Outfolder = "/project/shafferslab/Guillaume/10X_exp1_reanalysis/BarcodeProcessing/Barcode_output"

#length to keep from sequenced barcode (this is actual bacode, does not include strt seq below)
bclen = 70

#common sequence right before starcode starts
strtseq = "GGACGAGCTGTACAAGTAGG"

#allowed number of mismatches between barcodes to be called the same (starcode input)
sc_mm = 8

#Define which samples should be run together in starcode (any sample that could contain cells from the same lineage should be grouped together)
#NOTE: names used here should match the names that will be assinged to samples in line 77-81 in this script (should be everythign befor lane (_L0))
GSAMP = [["R1Enriched1_S1","R1Enriched2_S2"],["R1Mix1_S3","R1Mix2_S4"],["R2Enriched1_S5","R2Enriched2_S6"],["R2Mix1_S7","R2Mix2_S8"]]

#### END OF USER DEFINED PARAMETERS ##### NOTE: there are many more parameters that can be changed but this required you to go into the code below and understand all the different funtions

#define funtion to determine file has something in the first line (consider files with nothing in the first line as empty)
def empty(fname):
    with open(fname) as f:
       return f.readline() == ""

#define any new paths
sc_in = Outfolder + "/starcode_inputs/"
sc_out = Outfolder + "/starcode_outputs/"
mod_R2 = Outfolder + "/Modified_fastq/"
CellR = Outfolder + "/CellRanger_inputs/"
CellRfq = CellR + "FASTQ/"

# ADd starcode to PATH
os.environ["PATH"] += os.pathsep + scripts + '/starcode/'

# Make any necessary files
os.mkdir(Outfolder)
os.mkdir(sc_in)
os.mkdir(sc_out)
os.mkdir(mod_R2)
os.mkdir(CellR)
os.mkdir(CellRfq)

#unzip all files created by 10x for barcode runs
gunzipCommand = ['gunzip', '-r', Fastqfolder]
subprocess.call(gunzipCommand)
print("unziped")

### GET ALL BARCODES TO START WITH SAME SEQUENCE AND TRIM THE BARCODES

#get all Read2 fastq file paths
all_R2 = glob.glob(Fastqfolder + "/**/*_R2*.fastq", recursive = True)

#define samples
samples = []
for paths in all_R2:
    samples.append(paths.split("/")[-1].split("_L0")[0])
samples = list(set(samples))
samples.sort()

#loop through all the Read 2 fastq files (which contain barcode sequences)
for sample in samples:
    s_fastq = []
    for path in all_R2:
        if "/"+ sample in path:
            s_fastq.append(path)
    #get the sequences in all the fastqs
    #conatains all sequences for a sample
    seqs = []
    for fsmp in s_fastq:
        for record in SeqIO.parse(fsmp, "fastq"):
                seqs.append(str(record.seq))


    #determine which seqeunces are actual barcodes by determining which start with the constant sequence before the barcode

    # variables: the constant sequence before the barcode, percent match to call it a barcode (here it is 70)

    # since we stagger priming sight, here we determine what the correct stagger is for a given fastq file
    strt_ind = []
    for lines in seqs:
        c_ind = len(lines.split(strtseq)[0])
        if c_ind < len(lines):
            strt_ind.append(c_ind)

    bc_strt = statistics.mode(strt_ind)

    # modify sequences so that those that have a start sequence (where there can be erros) get replaces with a pefect start sequence, and get rid of stagger so that all sequence start with the perfect start sequences
    modseq1 = []
    for i in seqs:
        strt = i[bc_strt:len(strtseq)]
        pctmatch = (fuzz.ratio(strtseq,strt))

        if pctmatch >= 70:
            trim = i[bc_strt + len(strtseq) :]
            modseq1.append(strtseq + trim)
        else :
            modseq1.append(i)


    #trim barcodes
    # variable: how long do you want your barcode
    sc_input = []
    for i in modseq1:
        sc_input.append(i[0:bclen+len(strtseq)])

    #write files with these edited barcodes ( these are used as the input into starcode)
    f = open(sc_in + "/" "sc_input" + "_" + sample +".txt","w")
    f.write('\n'.join(sc_input))
    f.close()
    sleep(20)



# Get all of the files in the sc_input folder and concatenate them
for grp in GSAMP:
    count = 0
    with open(sc_in + '/sc_in_'+ "_".join(grp) + '_comb.txt', 'w') as outfile:

        for smp in grp:
            smpf = (sc_in + "sc_input_" + smp +".txt")
            with open(smpf) as infile:
                for line in infile:
                    outfile.write(line)
                count = count + 1
                if count < len(grp):
                    outfile.write("\n") # only add line if it is not the last sample

## RUN STARCODE
#get directories in files written in the step above
all_sc_in = glob.glob(sc_in + "sc_in_*_comb.txt")

#for all the starcode input files run starcode
for files in all_sc_in:
    outpath = sc_out + "sc_output_" + files.split("/")[-1].split("sc_in_")[-1]
    starcodeCommand = ['bsub', '-M','204800','-e' ,"Error_" + files.split("_")[-1], '-o', 'Output_' + files.split("_")[-1] ,'starcode', '-d', str(sc_mm), '-t', '20', '-i', files, '-o', outpath, '--seq-id']
    subprocess.call(starcodeCommand)
    sleep(1)

#Using a while loop to wait for starcode files to be produced. If starcodes arent appearing in a reasonable rate this errors out assuming an issue, but if there are large files you may need to change how long this code waits until erroring out

#wait 10 min to make sure that all the empty starcdoes files are generated
sleep(600)

#define paths of all sc outputs
all_sc_out = glob.glob(sc_out + "sc_output*.txt")

#wait for starcode to finish
while sum([empty(file) for file in all_sc_out]) > 0:
    sleep(600)


### REPLACE ORGINAL SEQUENCES WITH MODIFIED SEQUENCES IN FASTQ FILES

for grp in GSAMP:
    #define sc output file for this group of samples
    for path in all_sc_out:
        if "/sc_output_" + "_".join(grp) in path:
            scfile = path

    #read in starcode file and make a list of its lines
    starcode_file = open(scfile, "r")
    sc_lines = starcode_file.readlines()


    #Identify the fastq files that need to be concatenated
    R2paths = []
    for smp in grp:
        for path in all_R2:
            if "/"+ smp in path:
                R2paths.append(path)
    #read in fastq files of given group of samples (concateniating all fastq files from same sample into one list)
    cat_fastq = []
    l_fastq = []
    for p in R2paths:
        fastq_file = open(p, "r")
        l_fastq.append(p + "|" + str(len(cat_fastq)))
        for line in fastq_file:
            cat_fastq.append(line)



    for i in sc_lines:
        #get the barcode sequence in the line
        split_sc = i.split("\t")
        bcseq = split_sc[0]
        #get the list of lines that need to be replaces with bcseq
        rind = split_sc[2].split(",")
        rind = np.array(list(map(int, rind)))

        #convert index from sc_output file indexes to fastq file indexes
        rind=rind*4-3
        for j in rind:
            cat_fastq[j] = bcseq + "\n"
            #edit quality score to be the same length as the sequence
            cat_fastq[j+2] = cat_fastq[j+2][0:len(bcseq)] + "\n"


    # split list back out into individual fastq files
    fq_list = []
    for z in range(1,4*len(grp)):
        fq_list.append(cat_fastq[int(l_fastq[z-1].split("|")[-1]) : int(l_fastq[z].split("|")[-1])])
    fq_list.append(cat_fastq[int(l_fastq[z].split("|")[-1]) : ])



    #write edited fastq files
    #ind to determine which file to write
    ind = -1
    for el in l_fastq:
        ind = ind + 1
        n_fastq = open(mod_R2 + el.split("|")[0].split("/")[-1], "w")
        n_fastq.writelines(fq_list[ind])
        n_fastq.close()

### USING GENERATED FILES CREATE OUTPUTS TO RUN CELLRANGE FEATURE BARCODE ANALYSIS

#make Feature Reference CSV

#list of all final barcodes
fbc = []
for fname in all_sc_out:
    with open(fname) as infile:
        for line in infile:
            if strtseq in line:
                fbc.append(line[len(strtseq):].split("\t")[0])
#keep only unique barcode sequences
fbc = list(set(fbc))

#cnt for nameing barcodes
cnt = 0
with open(CellR + "FeatureReference.csv" , 'w+') as outfile:
    outfile.write("id,name,read,pattern,sequence,feature_type\n")
    for bc in fbc:
        cnt = cnt + 1
        full_line = "L" + str(cnt) + "," + "Lin" + str(cnt) + "," + "R2," + "5P" + strtseq + "(BC)," + bc + "," + "Custom\n"
        outfile.write(full_line)

#create fastq directory with the unchaged Read 1 files and the modified read2 files
all_modR2 = glob.glob(mod_R2 + "*_R2*.fastq", recursive = True)
#get all the fastq files that were not R2 files
all_nR2 = []
all_nR2.extend(glob.glob(Fastqfolder + "/**/*_R1*.fastq", recursive = True))
all_nR2.extend(glob.glob(Fastqfolder + "/**/*_I1*.fastq", recursive = True))
all_nR2.extend(glob.glob(Fastqfolder + "/**/*_I2*.fastq", recursive = True))

for sample in samples:
    sf = CellRfq + sample
    os.mkdir(sf)
    for R2file in all_modR2:
        if "/"+ sample in R2file:
            shutil.copy(R2file, sf)
    for R1file in all_nR2:
        if "/"+ sample in R1file:
            shutil.copy(R1file, sf)
