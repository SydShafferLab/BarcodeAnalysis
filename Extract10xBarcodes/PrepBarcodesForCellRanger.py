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
Fastqfolder="/project/shafferslab/Guillaume/10X_exp1_reanalysis/BarcodeProcessing/Barcode_Seq_Raw/BaseSpace_FASTQ"
#folder you want outputs go go into (dont make this folder, this scipt will make it)
Outfolder = "/project/shafferslab/Guillaume/10X_exp1_reanalysis/BarcodeProcessing/Barcode_output"

#length to keep from sequenced barcode (this is actual bacode, does not include strt seq below)
bclen = 60

#common sequence right before starcode starts
strtseq = "AGGACGAGCTGTACAAGTAGG"

#allowed number of mismatches between barcodes to be called the same (starcode input)
sc_mm = 8



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

#add starcode to PATH
os.environ["PATH"] += os.pathsep + scripts + '/starcode/'

#make any necessary files
os.mkdir(Outfolder)
os.mkdir(sc_in)
os.mkdir(sc_out)
os.mkdir(mod_R2)
os.mkdir(CellR)
os.mkdir(CellRfq)


### GET ALL BARCODES TO START WITH SAME SEQUENCE AND TRIM THE BARCODES

#get all Read2 fastq file paths
all_R2 = glob.glob(Fastqfolder + "/**/*_R2*.fastq", recursive = True)

#define samples
samples = []
for paths in all_R2:
    samples.append(paths.split("/")[-1].split("_")[0])
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


### RUN STARCODE
#get directories in files written in the step above
all_sc_in = glob.glob(sc_in + "sc_input*.txt")

#for all the starcode input files run starcode
for files in all_sc_in:
    outpath = sc_out + "sc_output_" + files.split("/")[-1].split("sc_input_")[-1]
    starcodeCommand = ['bsub', '-M','32000','-e' ,"Error_" + files.split("_")[-1], '-o', 'Output_' + files.split("_")[-1] ,'starcode', '-d', str(sc_mm), '-t', '20', '-i', files, '-o', outpath, '--seq-id']
    subprocess.call(starcodeCommand)
    sleep(.1)


# Using a while loop to wait for starcode files to be produced. If starcodes arent appearing in a reasonable rate this errors out assuming an issue, but if there are large files you may need to change how long this code waits until erroring out
all_sc_out = glob.glob(sc_out + "sc_output*.txt")

cnt = 0
while sum([empty(file) for file in all_sc_out]) > 0:
    #old number of empty files
    o_nef = sum([empty(file) for file in all_sc_out])
    sleep(300)
    #new number of empty files
    n_nef = sum([empty(file) for file in all_sc_out])

    #determine in the number of emty files is changing
    df = o_nef - n_nef

    if df <= 0 :
        cnt = cnt + 1
        if cnt >= 12 and df <= 0:
            print("ERROR: stracode files not being produced at expected rate, the input may be wrong or if inputing large files into starcod you may need to modify script to increase time allowed to produce these fiels")
            exit()

#### REPLACE ORGINAL SEQUENCES WITH MODIFIED SEQUENCES IN FASTQ FILES
all_sc_out = glob.glob(sc_out + "sc_output*.txt")


for sample in samples:
    #define sc output file for that sample
    for path in all_sc_out:
        if "/sc_output_"+ sample in path:
            scfile = path

    #read in starcode file and make a list of its lines
    starcode_file = open(scfile, "r")
    sc_lines = starcode_file.readlines()

    #read in fastq files of given sample (concateniating all fastq files from same sample into one list)
    cat_fastq = []
    l_fastq = []
    for path in all_R2:
        if "/"+ sample in path:
            fastq_file = open(path, "r")
            l_fastq.append(path + "|" + str(len(cat_fastq)))
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


    F1 = cat_fastq[int(l_fastq[0].split("|")[-1]) : int(l_fastq[1].split("|")[-1])-1]
    F2 = cat_fastq[int(l_fastq[1].split("|")[-1]) : int(l_fastq[2].split("|")[-1])-1]
    F3 = cat_fastq[int(l_fastq[2].split("|")[-1]) : int(l_fastq[3].split("|")[-1])-1]
    F4 = cat_fastq[int(l_fastq[3].split("|")[-1]) : ]

    fq_list = [F1,F2,F3,F4]

    #write edited fastq files
    #ind to determine which file to write
    ind = -1
    for el in l_fastq:
        print(el)
        ind = ind + 1
        n_fastq = open(mod_R2 + el.split("|")[0].split("/")[-1], "w")
        n_fastq.writelines(fq_list[ind])
        n_fastq.close()


### USING GENERATED FILES CREATE OUTPUTS TO RUN CELLRANGE FEATURE BARCODE ANALYSIS

#make Feature Reference CSV
all_scout = glob.glob(sc_out + "sc_output*.txt", recursive = True)
#cnt for nameing barcodes
cnt = 0
with open(CellR + "FeatureReference.csv" , 'w+') as outfile:
    for fname in all_scout:
        with open(fname) as infile:
            for line in infile:
                if strtseq in line:
                    cnt = cnt + 1
                    full_line = "L" + str(cnt) + "," + "Lin" + str(cnt) + "," + "R2," + "5P" + strtseq + "(BC)," + line.split("\t")[0] + "," + "Custom\n"
                    outfile.write(full_line)

#create fastq directory with the unchaged Read 1 files and the modified read2 files
all_modR2 = glob.glob(mod_R2 + "*_R2*.fastq", recursive = True)
all_R1 = glob.glob(Fastqfolder + "/**/*_R1*.fastq", recursive = True)

for sample in samples:
    sf = CellRfq + sample
    os.mkdir(sf)
    for R2file in all_modR2:
        if "/"+ sample in R2file:
            shutil.copy(R2file, sf)
    for R1file in all_R1:
        if "/"+ sample in R1file:
            shutil.copy(R1file, sf)
