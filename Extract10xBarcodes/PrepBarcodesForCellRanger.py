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
Fastqfolder="/project/shafferslab/Guillaume/10X_exp1/Barcode_Seq_Raw/BaseSpace_FASTQ"
#folder you want outputs go go into (dont make this folder, this scipt will make it)
Outfolder = "/project/shafferslab/Guillaume/10X_exp1/TestOut"

#length to keep from sequenced barcode (this is actual bacode, does not include strt seq below)
bclen = 60

#common sequence right before starcode starts
strtseq = "AGGACGAGCTGTACAAGTAGG"

#allowed number of mismatches between barcodes to be called the same (starcode input)
sc_mm = 8



#### END OF USER DEFINED PARAMETERS ##### NOTE: there are many more parameters that can be changed but this required you to go into the code below and understand all the different funtions


#define any new paths
sc_in = Outfolder + "/starcode_inputs/"
sc_out = Outfolder + "/starcode_outputs/"
mod_R2 = Outfolder + "/Modified_fastq/"
CellR = Outfolder + "/CellRanger_inputs/"
CellRfq = CellR + "/FASTQ/"

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

#loop through all the Read 2 fastq files (which contain barcode sequences)
for files in all_R2:
    seqs = []
    for record in SeqIO.parse(files, "fastq"):
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
    # may need to change this to actually extract the sample name
    f = open(sc_in + "/" "sc_input" + "_" + files.split("/")[-1].split(".")[0] +".txt","w")
    f.write('\n'.join(sc_input))


### RUN STARCODE
#get directories in files written in the step above
all_sc_in = glob.glob(sc_in + "sc_input*.txt")

#for all the starcode input files run starcode
for files in all_sc_in:
    outpath = sc_out + "sc_output_" +files.split("/")[-1].split("sc_input_")[-1]
    starcodeCommand = ['bsub', '-M','32000','-e' ,"Error_" + files.split("_")[-1], '-o', 'Output_' + files.split("_")[-1] ,'starcode', '-d', str(sc_mm), '-t', '20', '-i', files, '-o', outpath, '--seq-id']
    subprocess.call(starcodeCommand)

# Using a while loop to wait for starcode files to be produced. If starcodes arent appearing in a reasonable rate this errors out assuming an issue, but if there are large files you may need to change how long this code waits until erroring out
all_sc_out = glob.glob(sc_out + "sc_output*.txt")

cnt = 0
while (len(all_sc_out) < len(all_sc_in)):

    onf = len(all_sc_out)
    sleep(60)
    all_sc_out = glob.glob(sc_out + "sc_output*.txt")
    df = len(all_sc_out)-onf

    if df <= 0 :
        cnt = cnt + 1
        if cnt >= 10 and df <= 0:
            print("ERROR: stracode files not being produced at expected rate, the input may be wrong or if inputing large files into starcod you may need to modify script to increase time allowed to produce these fiels")
            exit()

#### REPLACE ORGINAL SEQUENCES WITH MODIFIED SEQUENCES IN FASTQ FILES
all_sc_out = glob.glob(sc_out + "sc_output*.txt")

for scfiles in all_sc_out:

    #read in starcode file and make a list of its lines
    starcode_file = open(scfiles, "r")
    sc_lines = starcode_file.readlines()

    #read in fastq file
    c_fastq = scfiles.split("/")[-1].split("sc_output_")[-1].split(".")[0]

    for R2 in all_R2:
        if c_fastq in R2:
            c_fastq_path = R2

    fastq_file= open(c_fastq_path, "r")
    fq_lines = fastq_file.readlines()

    for i in sc_lines:
        split_sc = i.split("\t")
        #get the barcode sequence in the line
        bcseq = split_sc[0]
        #get the list of lines that need to be replaces with bcseq
        rind = split_sc[2].split(",")
        rind = np.array(list(map(int, rind)))

        #convert index from sc_output file indexes to fastq file indexes
        rind=rind*4-3
        for j in rind:
            fq_lines[j] = bcseq + "\n"
            #edit quality score to be the same length as the sequence
            fq_lines[j+2] = fq_lines[j+2][0:len(bcseq)] + "\n"

    n_fastq = open(mod_R2 + scfiles.split("/")[-1].split("sc_output_")[-1].split(".")[0] + ".fastq", "w")
    n_fastq.writelines(fq_lines)
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

for R2file in all_modR2:
        shutil.copy(R2file, CellRfq)

for R1file in all_R1:
        shutil.copy(R1file, CellRfq)
