# Run this script by entering the 10xbc environment and then running: python /path/to/StepOnePrepBarcodesCellRanger.py

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from rapidfuzz import fuzz
import matplotlib as plt
from matplotlib import pyplot
import numpy as np
import glob
import os
import subprocess
from time import sleep
import statistics
import shutil
import re

from Functions.cDNA_10x_starcode_prep import cDNA_10x_starcode_prep
from Functions.gDNA_starcode_prep import gDNA_starcode_prep

from Functions.cDNA_10x_after_starcode import cDNA_10x_after_starcode
from Functions.gDNA_after_starcode import gDNA_after_starcode

print(" ")
print("Running...")
print(" ")


# Find paths_and_variables.json file
path_to_script = os.path.abspath(os.getcwd())

path_to_script = '/'.join(path_to_script.split('/')[:-1])
path = path = os.path.expanduser(path_to_script + "/paths_and_variables.json") 


# read paths_and_variables.json file
with open(path, 'r') as myfile:
    data=myfile.read()

result_dict =  jstyleson.loads(data) # Raise Exception

scripts=result_dict['scripts']         #path to all barcode analyssi scripts
Fastqfolder10x=result_dict['Fastqfolder10x'] #Folder that contains all folders containing FASTQ files generated from sequencing the barcodes
FastqfoldergDNA=result_dict['FastqfoldergDNA']#Folder that contains all folders containing gDNA FASTQ files generated from sequencing the barcodes.
Outfolder= result_dict['Outfolder']    #folder you want outputs go go into (dont make this folder, this scipt will make it)
strtseq= result_dict['strtseq']        #common sequence right before starcode starts
barcodeSource = result_dict['barcodeSource']    #determine whether the data has barcodes from 10x ("10x"), gDNA ("gDNA"), or both "both"
GSAMP= result_dict['GSAMP']            #Define which samples should be run together in starcode 
bclen = result_dict['bclen']           #length to keep from sequenced barcode 
strtseq =  result_dict['strtseq']      #common sequence right before starcode starts
strtseq_revcomp =  result_dict['strtseq_revcomp'] #rev_comp common sequence right before starcode starts
startseqMatch =  result_dict['startseqMatch']# The percentage match you for startseq to be called as correct in a barcode
sc_mm =  result_dict['sc_mm']          #allowed number of mismatches between barcodes to be called the same (starcode input)


#define funtion to determine file has something in the first line (consider files with nothing in the first line as empty)
def empty(fname):
    with open(fname) as f:
       return f.readline() == ""

#define funtion to determine if folder exist
def does_folder_exist(path_to_folder):
    if not os.path.exists(path_to_folder):
        os.mkdir(path_to_folder)
    else:
        raise Exception("folder {} already exists".format(path_to_folder))



#--------------------------------------------------------------------------
# Define all paths and create folders

#define any new paths
sc_in  = Outfolder + "/starcode_inputs/"
sc_out = Outfolder + "/starcode_outputs/"
mod_R2 = Outfolder + "/Modified_fastq/"
CellR  = Outfolder + "/CellRanger_inputs/"
CellRfq = CellR + "FASTQ/"

if barcodeSource == 'both' or barcodeSource == '10x':
    filt_haveStart_10x = Outfolder + '/fastq_with_startseq_10x'
    filt_WSN_10x = Outfolder + "/filtered_fastq_WSN_10x_Final_BC/"

if barcodeSource == 'both' or barcodeSource == 'gDNA':
    filt_haveStart_gDNA = Outfolder + '/fastq_with_startseq_gDNA'
    filt_WSN_gDNA = Outfolder + "/filtered_fastq_WSN_gDNA_Final_BC/"


# Add starcode to PATH
os.environ["PATH"] += os.pathsep + scripts + '/starcode/'


# Make any necessary files
path_to_folders = [sc_in,sc_out,mod_R2,filt_WSN_gDNA,CellR,CellRfq]

if barcodeSource == 'both' or barcodeSource == '10x':
    path_to_folders.append(filt_haveStart_10x,filt_haveStart_gDNA)

if barcodeSource == 'both' or barcodeSource == 'gDNA':
    path_to_folders.append(filt_haveStart_gDNA,filt_WSN_gDNA)

# checking whether folder/directory exists
for path_to_folder in path_to_folders:
    does_folder_exist(path_to_folder)


if barcodeSource == 'both' or barcodeSource == '10x':
    #unzip all files created by 10x for barcode runs
    gunzipCommand = ['gunzip', '-r', Fastqfolder10x]
    subprocess.call(gunzipCommand)

if barcodeSource == 'both' or barcodeSource == 'gDNA':
    #unzip all files created by 10x for barcode runs
    gunzipCommand = ['gunzip', '-r', FastqfoldergDNA]
    subprocess.call(gunzipCommand)








#--------------------------------------------------------------------------

## Set up sequences for starcode

print("     Setting up sequences for starcode")
print(" ")

if barcodeSource == 'both' or barcodeSource == '10x':

    cDNA_10x_starcode_prep(Fastqfolder10x,GSAMP,bclen,filt_haveStart_10x)


if barcodeSource == 'both' or barcodeSource == 'gDNA':
 
    gDNA_starcode_prep(FastqfoldergDNA,GSAMP,bclen)








#--------------------------------------------------------------------------

## Get all of the files in the sc_input folder and concatenate them
counter = 1
for grp in GSAMP:
    count = 0
    with open(sc_in + '/sc_in_group'+ str(counter) + '_comb.txt', 'w') as outfile:

        for smp in grp:
            smpf = (sc_in + "sc_input_" + smp +".txt")
            with open(smpf) as infile:
                for line in infile:
                    outfile.write(line)
                count = count + 1
                if count < len(grp):
                    outfile.write("\n") # only add line if it is not the last sample
    counter = counter+1









#--------------------------------------------------------------------------

## RUN STARCODE
print("     Running starcode")
print(" ")

#get directories in files written in the step above
all_sc_in = glob.glob(sc_in + "sc_in_*_comb.txt")

#for all the starcode input files run starcode
for files in all_sc_in:
    outpath = sc_out + "sc_output_" + files.split("/")[-1].split("sc_in_")[-1]
    starcodeCommand = ['starcode', '-d', str(sc_mm), '-t', '20', '-i', files, '-o', outpath, '--seq-id']
    subprocess.call(starcodeCommand)
    sleep(1)

#define paths of all sc outputs
all_sc_out = glob.glob(sc_out + "sc_output*.txt")

#wait for starcode to finish
while sum([empty(file) for file in all_sc_out]) > 0:
    sleep(10)









#--------------------------------------------------------------------------
## get all Read1 or Read2 fastq file paths after starcode for gDNA or cDNA


if barcodeSource == 'both' or barcodeSource == '10x':

    #get all Read2 fastq file paths after starcode for 10x
    all_R2_10x_start_unfilt = glob.glob(filt_haveStart_10x + "/**/*_R2*.fastq", recursive = True)
    all_R2_10x_start_unfilt.sort()

    # Remove any Read2 fastq files that you dont care about
    all_R2_10x_start_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_R2_10x_start_unfilt):
                all_R2_10x_start_temp.append([s for s in all_R2_10x_start_unfilt if str(smp) in s])

    all_R2_10x_start = [item for sublist in all_R2_10x_start_temp for item in sublist]

    #define samples
    samples_R2_10x_start = []
    for paths in all_R2_10x_start:
        samples_R2_10x_start.append(paths.split("/")[-1].split("_L0")[0])
    samples_R2_10x_start = list(set(samples_R2_10x_start))
    samples_R2_10x_start.sort()


if barcodeSource == 'both' or barcodeSource == 'gDNA':  
    #get all Read1 fastq file paths after starcode for gDNA
    all_R1_gDNA_start_unfilt = glob.glob(filt_haveStart_gDNA + "/**/*_R1*.fastq", recursive = True)

    # Remove any Read2 fastq files that you dont care about
    all_R1_gDNA_start_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_R1_gDNA_start_unfilt):
                all_R1_gDNA_start_temp.append([s for s in all_R1_gDNA_start_unfilt if str(smp) in s])

    all_R1_gDNA_start = [item for sublist in all_R1_gDNA_start_temp for item in sublist]

    #define samples
    samples_R1_gDNA_start = []
    for paths in all_R1_gDNA_start:
        samples_R1_gDNA_start.append(paths.split("/")[-1].split("_L0")[0])
    samples_R1_gDNA_start = list(set(samples_R1_gDNA_start))
    samples_R1_gDNA_start.sort()










#--------------------------------------------------------------------------

### REPLACE ORGINAL SEQUENCES WITH MODIFIED SEQUENCES IN FASTQ FILES
print("     Modifying fastq files")
print(" ")

counter = 1
for grp in GSAMP:
    count = 0
    #define sc output file for this group of samples
    for path in all_sc_out:
        if "/sc_output_group" + str(counter) +"_comb.txt" in path:
            scfile = path
            print(scfile)
    counter = counter + 1
    
    #read in starcode file and make a list of its lines
    starcode_file = open(scfile, "r")
    sc_lines = starcode_file.readlines()


    #Identify the fastq files that need to be concatenated
    fastq_paths = []
    for smp in grp:

        if barcodeSource == 'both' or barcodeSource == '10x':
            for path in all_R2_10x_start:
                if "/"+ smp in path:
                    fastq_paths.append(path)
        if barcodeSource == 'both' or barcodeSource == 'gDNA':     
            for path in all_R1_gDNA_start:
                if "/"+ smp in path:
                    fastq_paths.append(path)
                
    
    #read in fastq files of given group of samples (concatenating all fastq files from same sample into one list)
    cat_fastq = []
    l_fastq = []
    for p in fastq_paths:
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











#--------------------------------------------------------------------------

### USING GENERATED FILES CREATE OUTPUTS TO RUN CELLRANGER FEATURE BARCODE ANALYSIS



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

#cnt for naming barcodes
cnt = 0
with open(CellR + "FeatureReference.csv" , 'w+') as outfile:
    outfile.write("id,name,read,pattern,sequence,feature_type\n")
    for bc in fbc:
        cnt = cnt + 1
        full_line = "L" + str(cnt) + "," + "Lin" + str(cnt) + "," + "R2," + "5P" + strtseq + "(BC)," + bc + "," + "Custom\n"
        outfile.write(full_line)


# If cDNA samples are used, filter
if barcodeSource == 'both' or barcodeSource == '10x':
    cDNA_10x_after_starcode(mod_R2,GSAMP,filt_WSN_10x)

# If gDNA samples are used, filter
if barcodeSource == 'both' or barcodeSource == 'cDNA':
    gDNA_after_starcode(mod_R2,GSAMP,filt_WSN_gDNA)













#--------------------------------------------------------------------------
print("     Filtering barcodes based on WSN")
print(" ")


# Filter the feature reference file to remove all of the barcodes that do not file WSN

#read in feature reference file and make a list of its lines
fr_file = CellR + "FeatureReference.csv"

featureReference_file = open(fr_file, "r")
featureReference_lines = featureReference_file.readlines()

iters = 0
with open(CellR + "FeatureReference_filtered.csv" , 'w+') as outfile:
    outfile.write(featureReference_lines[0]) # Adds the titles back in
    for line in featureReference_lines[1:]:
        if 'AAAA' in line.split(',')[4]:
            sleep(0)
        elif 'TTTT' in line.split(',')[4]:
            sleep(0)
        elif 'GGGG' in line.split(',')[4]:
            sleep(0)
        elif 'CCCC' in line.split(',')[4]:
            sleep(0)
        elif 'NN' in line.split(',')[4]:
            sleep(0)
        else:
            outfile.write(line)
        
        iters = iters+1


if barcodeSource == 'both' or barcodeSource == '10x':
    #create fastq directory with all of the modified files for starcode
    all_modR2 = glob.glob(filt_WSN_10x + "*_R2*.fastq", recursive = True)

    #get all the modified fastq files that were not R2 files
    all_nR2 = []
    all_nR2.extend(glob.glob(filt_WSN_10x + "/**/*_R1*.fastq", recursive = True))
    all_nR2.extend(glob.glob(filt_WSN_10x + "/**/*_I1*.fastq", recursive = True))
    all_nR2.extend(glob.glob(filt_WSN_10x + "/**/*_I2*.fastq", recursive = True))

    for sample in samples_R2_10x:
        sf = CellRfq + sample
        os.mkdir(sf)
        for R2file in all_modR2:
            if "/"+ sample in R2file:
                shutil.copy(R2file, sf)
        for R1file in all_nR2:
            if "/"+ sample in R1file:
                shutil.copy(R1file, sf)







#--------------------------------------------------------------------------

# Properly format the gDNA output files
# Currently is output in --seq-id format which means the read (by numerical order) that goes into 
#sequence is listed after

count = 1
for grp in GSAMP:
    #define sc output file for this group of samples
    for path in all_sc_out:
        if "/sc_output_group" + str(count) +"_comb.txt" in path:
            scfile = path

    count = count+1
    #read in starcode file and make a list of its lines
    starcode_file = open(scfile, "r")
    sc_lines = starcode_file.readlines()
    
    ind_starts = []
    for i in l_fastq:
        ind_starts.append(i.split('|')[1])
    ind_starts = [(int(i)/4)+1 for i in ind_starts]
    ind_starts.append(len(cat_fastq)/4+1) # Have to add the one since starcode outputs are 1-based instead of 0-based
    
    counter = -1
    for smp in grp:
        counter = counter+1
        ind_to_test = counter*4

        with open(sc_out + "sc_output_counts_" + smp +".txt", 'w') as outfile:
            for i in sc_lines:
                #get the barcode sequence in the line
                split_sc = i.split("\t")
                bcseq = split_sc[0]
                #get the list of lines that need to be replaces with bcseq
                rind = split_sc[2].split(",")
                rind = np.array(list(map(int, rind)))
                num_in_bounds_mask = np.logical_and(rind >= int(ind_starts[ind_to_test]), rind < int(ind_starts[ind_to_test+4]))
                num_in_bounds = rind[num_in_bounds_mask]
                outfile.write("{}\t{}\n".format(bcseq, len(num_in_bounds)))





#--------------------------------------------------------------------------
## Remove old location for Cellranger input fastq files
sleep(10)
os.remove(filt_WSN_10x)


print("     Step 3 is done :D ")
print(" ")

