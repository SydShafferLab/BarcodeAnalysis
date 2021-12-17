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
barcodeSource = result_dict['barcodeSource']    #determine whether the data has barcodes from 10x and gDNA "both" or only from 10x "10x"
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


#define any new paths
sc_in = Outfolder + "/starcode_inputs/"
sc_out = Outfolder + "/starcode_outputs/"
mod_R2 = Outfolder + "/Modified_fastq/"
filt_haveStart_10x = Outfolder + '/fastq_with_startseq_10x'
if barcodeSource == 'both':
    filt_haveStart_gDNA = Outfolder + '/fastq_with_startseq_gDNA'
filt_WSN_10x = Outfolder + "/filtered_fastq_WSN_10x_Final_BC/"
if barcodeSource == 'both':
    filt_WSN_gDNA = Outfolder + "/filtered_fastq_WSN_gDNA_Final_BC/"
CellR = Outfolder + "/CellRanger_inputs/"
CellRfq = CellR + "FASTQ/"

# Add starcode to PATH
os.environ["PATH"] += os.pathsep + scripts + '/starcode/'



# Make any necessary files
path_to_folders = [sc_in,sc_out,mod_R2,filt_haveStart_10x,filt_haveStart_gDNA,filt_WSN_gDNA,CellR,CellRfq]

if barcodeSource == 'both':
    path_to_folders.append(filt_haveStart_gDNA,filt_WSN_gDNA)

# checking whether folder/directory exists
for path_to_folder in path_to_folders:
    does_folder_exist(path_to_folder)

#unzip all files created by 10x for barcode runs
gunzipCommand = ['gunzip', '-r', Fastqfolder10x]
subprocess.call(gunzipCommand)
print("unzipped")



#get all Read2 fastq file paths
all_R2_10x_unfilt = glob.glob(Fastqfolder10x + "/**/*_R2*.fastq", recursive = True)
all_R2_10x_unfilt.sort()

# Remove any Read2 fastq files that you dont care about
all_R2_10x_temp = [] 
for grp in GSAMP:
    for smp in grp:
        if str(smp) in str(all_R2_10x_unfilt):
            all_R2_10x_temp.append([s for s in all_R2_10x_unfilt if str(smp) in s])

all_R2_10x = [item for sublist in all_R2_10x_temp for item in sublist]

#define samples
samples_R2_10x = []
for paths in all_R2_10x:
    print(paths)
    samples_R2_10x.append(paths.split("/")[-1].split("_L0")[0])
samples_R2_10x = list(set(samples_R2_10x))
samples_R2_10x.sort()



#loop through all the Read 2 fastq files from 10x (which contain barcode sequences) and write the textfile for starcode inputs
readsKeep_10x = []
counter = 0
file_endpoints_10x = []
for sample in samples_R2_10x:
    file_endpoints_10x.append([])
    file_endpoints_10x[counter].insert(0,0)
    print(sample)
    s_fastq = []
    for path in all_R2_10x:
        if "/"+ sample in path:
            s_fastq.append(path)
    # get JUST the sequences in all the fastqs contains all sequences for a sample
    seqs = []
    for fsmp in s_fastq:
        for record in SeqIO.parse(fsmp, "fastq"):
            seqs.append(str(record.seq))
        file_endpoints_10x[counter].append(len(seqs))

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
    modseq1 = [] # For writing just the good sequences
    modseq2 = [] # For writing good and bad sequnces
    readsKeep_10x.append([])
    for i in seqs:
        strt = i[bc_strt:(len(strtseq)+bc_strt)]
        pctmatch = (fuzz.ratio(strtseq,strt))

        if pctmatch >= startseqMatch:
            trim = i[bc_strt + len(strtseq) :]
            modseq1.append(strtseq + trim)
            modseq2.append(strtseq + trim)
            readsKeep_10x[counter].append(True)
        else :
            readsKeep_10x[counter].append(False)
            modseq2.append(i)

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
    
    
    # get all of the fastq lines from the fastqs
    lines = [] 
    for fsmp in s_fastq:
        with open(fsmp) as f:
            for line in f:
                lines.append(line)
                
    # Write the edited sequence and trimmed file  to a fastq folder 
    lane_counter = 0
    for fsamp in s_fastq: 
        print(fsamp)
        #write files with these edited barcodes ( these are used as the input into starcode)
        with open(filt_haveStart_10x + "/" + fsamp.split('/')[-1],"w") as f:
            print('Beginning count: ' + str(file_endpoints_10x[counter][lane_counter]))
            print('End count: '+ str(file_endpoints_10x[counter][lane_counter+1]))
            for j in range(file_endpoints_10x[counter][lane_counter],file_endpoints_10x[counter][lane_counter+1]):
                if readsKeep_10x[counter][j] == 1:
                    f.write(lines[(4*j)])
                    f.write(modseq2[j][0:bclen+len(strtseq)]+'\n') # pulls the trimmed sequence
                    f.write(lines[(4*j)+2])
                    f.write(lines[(4*j)+3][bc_strt:(bclen+len(strtseq)+bc_strt)]+'\n')
        
        lane_counter = lane_counter + 1

    counter = counter+1
  
#get all 10x Read1 fastq file paths
all_R1_10x_unfilt = glob.glob(Fastqfolder10x + "/**/*_R1*.fastq", recursive = True)
all_R1_10x_unfilt.sort()

# Remove any Read1 fastq files that you dont care about
all_R1_10x_temp = [] 
for grp in GSAMP:
    for smp in grp:
        if str(smp) in str(all_R1_10x_unfilt):
            all_R1_10x_temp.append([s for s in all_R1_10x_unfilt if str(smp) in s])

all_R1_10x = [item for sublist in all_R1_10x_temp for item in sublist]

#define samples Read1
samples_R1_10x = []
for paths in all_R1_10x:
    print(paths)
    samples_R1_10x.append(paths.split("/")[-1].split("_L0")[0])
samples_R1_10x = list(set(samples_R1_10x))
samples_R1_10x.sort()

# Loop through the corresponding Read1 files and remove the bad reads 
counter = 0
for sample in samples_R1_10x:
    s_fastq = []
    for path in all_R1_10x:
        if "/"+ sample in path:
            s_fastq.append(path)

    # get all of the fastq lines from the fastqs
    lines = [] 
    for fsmp in s_fastq:
        with open(fsmp) as f:
            for line in f:
                lines.append(line)
                
    # Write the edited sequence and trimmed file  to a fastq folder 
    lane_counter = 0
    for fsamp in s_fastq: 
        print(fsamp)
        #write files with these edited barcodes ( these are used as the input into starcode)
        with open(filt_haveStart_10x + "/" + fsamp.split('/')[-1],"w") as f:
            print('Beginning count: ' + str(file_endpoints_10x[counter][lane_counter]))
            print('End count: '+ str(file_endpoints_10x[counter][lane_counter+1]))
            for j in range(file_endpoints_10x[counter][lane_counter],file_endpoints_10x[counter][lane_counter+1]):
                if readsKeep_10x[counter][j] == 1:
                    f.write(lines[(4*j)])
                    f.write(lines[(4*j)+1]) # pulls the trimmed sequence
                    f.write(lines[(4*j)+2])
                    f.write(lines[(4*j)+3])
        
        lane_counter = lane_counter + 1

    
    counter = counter+1

#get all 10x Index1 fastq file paths
all_I1_10x_unfilt = glob.glob(Fastqfolder10x + "/**/*_I1*.fastq", recursive = True)
all_I1_10x_unfilt.sort()
# Remove any Read1 fastq files that you dont care about
all_I1_10x_temp = [] 
for grp in GSAMP:
    for smp in grp:
        if str(smp) in str(all_I1_10x_unfilt):
            all_I1_10x_temp.append([s for s in all_I1_10x_unfilt if str(smp) in s])

all_I1_10x = [item for sublist in all_I1_10x_temp for item in sublist]

#define samples Read1
samples_I1_10x = []
for paths in all_I1_10x:
    print(paths)
    samples_I1_10x.append(paths.split("/")[-1].split("_L0")[0])
samples_I1_10x = list(set(samples_I1_10x))
samples_I1_10x.sort()

# Loop through the corresponding index1 files and remove the bad reads 
counter = 0
for sample in samples_I1_10x:
    s_fastq = []
    for path in all_I1_10x:
        if "/"+ sample in path:
            s_fastq.append(path)

    # get all of the fastq lines from the fastqs
    lines = [] 
    for fsmp in s_fastq:
        with open(fsmp) as f:
            for line in f:
                lines.append(line)
                
    # Write the edited sequence and trimmed file  to a fastq folder 
    lane_counter = 0
    for fsamp in s_fastq: 
        print(fsamp)
        #write files with these edited barcodes ( these are used as the input into starcode)
        with open(filt_haveStart_10x + "/" + fsamp.split('/')[-1],"w") as f:
            print('Beginning count: ' + str(file_endpoints_10x[counter][lane_counter]))
            print('End count: '+ str(file_endpoints_10x[counter][lane_counter+1]))
            for j in range(file_endpoints_10x[counter][lane_counter],file_endpoints_10x[counter][lane_counter+1]):
                if readsKeep_10x[counter][j] == 1:
                    f.write(lines[(4*j)])
                    f.write(lines[(4*j)+1]) # pulls the trimmed sequence
                    f.write(lines[(4*j)+2])
                    f.write(lines[(4*j)+3])
        
        lane_counter = lane_counter + 1

    
    counter = counter+1
  
#get all 10x Index2 fastq file paths
all_I2_10x_unfilt = glob.glob(Fastqfolder10x + "/**/*_I2*.fastq", recursive = True)
all_I2_10x_unfilt.sort()

# Remove any Read1 fastq files that you dont care about
all_I2_10x_temp = [] 
for grp in GSAMP:
    for smp in grp:
        if str(smp) in str(all_I2_10x_unfilt):
            all_I2_10x_temp.append([s for s in all_I2_10x_unfilt if str(smp) in s])

all_I2_10x = [item for sublist in all_I2_10x_temp for item in sublist]

#define samples Read1
samples_I2_10x = []
for paths in all_I2_10x:
    print(paths)
    samples_I2_10x.append(paths.split("/")[-1].split("_L0")[0])
samples_I2_10x = list(set(samples_I2_10x))
samples_I2_10x.sort()

# Loop through the corresponding index2 files and remove the bad reads 
counter = 0
for sample in samples_I2_10x:
    s_fastq = []
    for path in all_I2_10x:
        if "/"+ sample in path:
            s_fastq.append(path)

    # get all of the fastq lines from the fastqs
    lines = [] 
    for fsmp in s_fastq:
        with open(fsmp) as f:
            for line in f:
                lines.append(line)
                
    # Write the edited sequence and trimmed file  to a fastq folder 
    lane_counter = 0
    for fsamp in s_fastq: 
        print(fsamp)
        #write files with these edited barcodes ( these are used as the input into starcode)
        with open(filt_haveStart_10x + "/" + fsamp.split('/')[-1],"w") as f:
            print('Beginning count: ' + str(file_endpoints_10x[counter][lane_counter]))
            print('End count: '+ str(file_endpoints_10x[counter][lane_counter+1]))
            for j in range(file_endpoints_10x[counter][lane_counter],file_endpoints_10x[counter][lane_counter+1]):
                if readsKeep_10x[counter][j] == 1:
                    f.write(lines[(4*j)])
                    f.write(lines[(4*j)+1]) # pulls the trimmed sequence
                    f.write(lines[(4*j)+2])
                    f.write(lines[(4*j)+3])
        
        lane_counter = lane_counter + 1

    
    counter = counter+1

if barcodeSource == 'both':
    #get all gDNA Read1 fastq file paths
    all_R1_gDNA_unfilt = glob.glob(FastqfoldergDNA + "/**/*_R1*.fastq", recursive = True)
    all_R1_gDNA_unfilt.sort()
    
    # Remove any Read1 fastq files that you dont care about
    all_R1_gDNA_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_R1_gDNA_unfilt):
                all_R1_gDNA_temp.append([s for s in all_R1_gDNA_unfilt if str(smp) in s])

    all_R1_gDNA = [item for sublist in all_R1_gDNA_temp for item in sublist]

    #define samples Read1
    samples_R1_gDNA = []
    for paths in all_R1_gDNA:
        print(paths)
        samples_R1_gDNA.append(paths.split("/")[-1].split("_L0")[0])
    samples_R1_gDNA = list(set(samples_R1_gDNA))
    samples_R1_gDNA.sort()
    
    
    #loop through all the Read 1 fastq files from gDNA (which contain barcode sequences) and write the textfile for starcode inputs
    readsKeep_gDNA = []
    counter = 0
    file_endpoints_gDNA = []
    for sample in samples_R1_gDNA:
        file_endpoints_gDNA.append([])
        file_endpoints_gDNA[counter].insert(0,0)
        print(sample)
        s_fastq = []
        for path in all_R1_gDNA:
            if "/"+ sample in path:
                s_fastq.append(path)
        # get JUST the sequences in all the fastqs contains all sequences for a sample
        seqs = []
        for fsmp in s_fastq:
            for record in SeqIO.parse(fsmp, "fastq"):
                seqs.append(str(record.seq.reverse_complement()))
            file_endpoints_gDNA[counter].append(len(seqs))

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
        modseq1 = [] # For writing just the good sequences
        modseq2 = [] # For writing good and bad sequnces
        readsKeep_gDNA.append([])
        for i in seqs:
            strt = i[bc_strt:(len(strtseq)+bc_strt)]
            pctmatch = (fuzz.ratio(strtseq,strt))

            if pctmatch >= startseqMatch:
                trim = i[bc_strt + len(strtseq) :]
                modseq1.append(strtseq + trim)
                modseq2.append(strtseq + trim)
                readsKeep_gDNA[counter].append(True)
            else :
                readsKeep_gDNA[counter].append(False)
                modseq2.append(i)

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
    
    
        # get all of the fastq lines from the fastqs
        lines = [] 
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)
                
        # Write the edited sequence and trimmed file  to a fastq folder 
        lane_counter = 0
        for fsamp in s_fastq: 
            print(fsamp)
            #write files with these edited barcodes ( these are used as the input into starcode)
            with open(filt_haveStart_gDNA + "/" + fsamp.split('/')[-1],"w") as f:
                print('Beginning count: ' + str(file_endpoints_gDNA[counter][lane_counter]))
                print('End count: '+ str(file_endpoints_gDNA[counter][lane_counter+1]))
                for j in range(file_endpoints_gDNA[counter][lane_counter],file_endpoints_gDNA[counter][lane_counter+1]):
                    if readsKeep_gDNA[counter][j] == 1:
                        f.write(lines[(4*j)])
                        f.write(modseq2[j][0:bclen+len(strtseq)]+'\n') # pulls the trimmed sequence
                        f.write(lines[(4*j)+2])
                        f.write(lines[(4*j)+3][::-1][bc_strt+1:(bclen+len(strtseq)+bc_strt)+1]+'\n') # Also reverses the phred score sequence
        
            lane_counter = lane_counter + 1

    
        counter = counter+1
        
    
    
    #get all gDNA Index1 fastq file paths
    all_I1_gDNA_unfilt = glob.glob(FastqfoldergDNA + "/**/*_I1*.fastq", recursive = True)
    all_I1_gDNA_unfilt.sort()

    # Remove any Read1 fastq files that you dont care about
    all_I1_gDNA_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_I1_gDNA_unfilt):
                all_I1_gDNA_temp.append([s for s in all_I1_gDNA_unfilt if str(smp) in s])

    all_I1_gDNA = [item for sublist in all_I1_gDNA_temp for item in sublist]

    #define samples Read1
    samples_I1_gDNA = []
    for paths in all_I1_gDNA:
        print(paths)
        samples_I1_gDNA.append(paths.split("/")[-1].split("_L0")[0])
    samples_I1_gDNA = list(set(samples_I1_gDNA))
    samples_I1_gDNA.sort()
    
    # Loop through the corresponding index1 files and remove the bad reads 
    counter = 0
    for sample in samples_I1_gDNA:
        s_fastq = []
        for path in all_I1_gDNA:
            if "/"+ sample in path:
                s_fastq.append(path)

        # get all of the fastq lines from the fastqs
        lines = [] 
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)
                
        # Write the edited sequence and trimmed file  to a fastq folder 
        lane_counter = 0
        for fsamp in s_fastq: 
            print(fsamp)
            #write files with these edited barcodes ( these are used as the input into starcode)
            with open(filt_haveStart_gDNA + "/" + fsamp.split('/')[-1],"w") as f:
                print('Beginning count: ' + str(file_endpoints_gDNA[counter][lane_counter]))
                print('End count: '+ str(file_endpoints_gDNA[counter][lane_counter+1]))
                for j in range(file_endpoints_gDNA[counter][lane_counter],file_endpoints_gDNA[counter][lane_counter+1]):
                    if readsKeep_gDNA[counter][j] == 1:
                        f.write(lines[(4*j)])
                        f.write(lines[(4*j)+1]) # pulls the trimmed sequence
                        f.write(lines[(4*j)+2])
                        f.write(lines[(4*j)+3])
        
            lane_counter = lane_counter + 1

    
        counter = counter+1
        
    #get all gDNA Index2 fastq file paths
    all_I2_gDNA_unfilt = glob.glob(FastqfoldergDNA + "/**/*_I2*.fastq", recursive = True)
    all_I2_gDNA_unfilt.sort()

    # Remove any Read1 fastq files that you dont care about
    all_I2_gDNA_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_I2_gDNA_unfilt):
                all_I2_gDNA_temp.append([s for s in all_I2_gDNA_unfilt if str(smp) in s])

    all_I2_gDNA = [item for sublist in all_I2_gDNA_temp for item in sublist]

    #define samples Read1
    samples_I2_gDNA = []
    for paths in all_I2_gDNA:
        print(paths)
        samples_I2_gDNA.append(paths.split("/")[-1].split("_L0")[0])
    samples_I2_gDNA = list(set(samples_I2_gDNA))
    samples_I2_gDNA.sort()
    
    # Loop through the corresponding index1 files and remove the bad reads 
    counter = 0
    for sample in samples_I2_gDNA:
        s_fastq = []
        for path in all_I2_gDNA:
            if "/"+ sample in path:
                s_fastq.append(path)

        # get all of the fastq lines from the fastqs
        lines = [] 
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)
                
        # Write the edited sequence and trimmed file  to a fastq folder 
        lane_counter = 0
        for fsamp in s_fastq: 
            print(fsamp)
            #write files with these edited barcodes ( these are used as the input into starcode)
            with open(filt_haveStart_gDNA + "/" + fsamp.split('/')[-1],"w") as f:
                print('Beginning count: ' + str(file_endpoints_gDNA[counter][lane_counter]))
                print('End count: '+ str(file_endpoints_gDNA[counter][lane_counter+1]))
                for j in range(file_endpoints_gDNA[counter][lane_counter],file_endpoints_gDNA[counter][lane_counter+1]):
                    if readsKeep_gDNA[counter][j] == 1:
                        f.write(lines[(4*j)])
                        f.write(lines[(4*j)+1]) # pulls the trimmed sequence
                        f.write(lines[(4*j)+2])
                        f.write(lines[(4*j)+3]) 
        
            lane_counter = lane_counter + 1

    
        counter = counter+1
  
# Get all of the files in the sc_input folder and concatenate them
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

## RUN STARCODE
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
    print(paths)
    samples_R2_10x_start.append(paths.split("/")[-1].split("_L0")[0])
samples_R2_10x_start = list(set(samples_R2_10x_start))
samples_R2_10x_start.sort()

print(samples_R2_10x_start)

if barcodeSource == 'both':  
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
        print(paths)
        samples_R1_gDNA_start.append(paths.split("/")[-1].split("_L0")[0])
    samples_R1_gDNA_start = list(set(samples_R1_gDNA_start))
    samples_R1_gDNA_start.sort()

    print(samples_R1_gDNA_start)

### REPLACE ORGINAL SEQUENCES WITH MODIFIED SEQUENCES IN FASTQ FILES
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
        print(smp)
        for path in all_R2_10x_start:
            if "/"+ smp in path:
                fastq_paths.append(path)
        if barcodeSource == 'both':     
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

#get all Read2 fastq file paths after starcode
all_R2_10x_mod_unfilt = glob.glob(mod_R2 + "/**/*_R2*.fastq", recursive = True)
all_R2_10x_mod_unfilt.sort()

# Remove any Read2 fastq files that you dont care about
all_R2_10x_mod_temp = [] 
for grp in GSAMP:
    for smp in grp:
        if str(smp) in str(all_R2_10x_mod_unfilt):
            all_R2_10x_mod_temp.append([s for s in all_R2_10x_mod_unfilt if str(smp) in s])

all_R2_10x_mod = [item for sublist in all_R2_10x_mod_temp for item in sublist]

#define samples
samples_R2_10x_mod = []
for paths in all_R2_10x_mod:
    print(paths)
    samples_R2_10x_mod.append(paths.split("/")[-1].split("_L0")[0])
samples_R2_10x_mod = list(set(samples_R2_10x_mod))
samples_R2_10x_mod.sort()

print(samples_R2_10x_mod)

# for 10x samples, build a boolean of sequences that follow WSN nucleotide motif
file_endpoints_10x_postsc = []
bool_WSN_10x = []
counter = 0
for sample in samples_R2_10x_mod:
    print(sample)
    s_fastq = []
    for path in all_R2_10x_mod:
        if "/"+ sample in path:
            s_fastq.append(path)
            
    print(len(s_fastq))
    
    # Build boolean that sequences with runs of 4 or containing N are not used
    bool_WSN_10x.append([])
    file_endpoints_10x_postsc.append([])
    seqs = []
    for fsmp in s_fastq:
        for record in SeqIO.parse(fsmp, "fastq"):
            if 'AAAA' in record.seq:
                bool_WSN_10x[counter].append(0)
            elif 'TTTT' in record.seq:
                bool_WSN_10x[counter].append(0)
            elif 'GGGG' in record.seq:
                bool_WSN_10x[counter].append(0)
            elif 'CCCC' in record.seq:
                bool_WSN_10x[counter].append(0)
            elif 'NN' in record.seq:
                bool_WSN_10x[counter].append(0)
            else:
                bool_WSN_10x[counter].append(1)
            seqs.append(record.seq)
        print(len(seqs))
        file_endpoints_10x_postsc[counter].append(len(seqs))
        
    file_endpoints_10x_postsc[counter].insert(0,0)
    counter = counter + 1

# Use the booleans to write new files only using only good read2s in 10x samples
counter = 0
for sample in samples_R2_10x_mod:
    print(sample)
    s_fastq = []
    for path in all_R2_10x_mod:
        if "/"+ sample in path:
            s_fastq.append(path)
    
    # get all of the fastq lines from the modified fastqs
    lines = [] # advance counter as I read lines and only append those where the boolean is 1
    for fsmp in s_fastq:
        with open(fsmp) as f:
            for line in f:
                lines.append(line)

    lane_counter = 0
    for fsamp in s_fastq: 
        print(fsamp)
        #write files with these edited barcodes ( these are used as the input into starcode)
        f = open(filt_WSN_10x + "/" + fsamp.split('/')[-1],"w")
        print('Beginning count: ' + str(file_endpoints_10x_postsc[counter][lane_counter]))
        print('End count: '+ str(file_endpoints_10x_postsc[counter][lane_counter+1]))
        for j in range(file_endpoints_10x_postsc[counter][lane_counter],file_endpoints_10x_postsc[counter][lane_counter+1]):
            #add print statements
            if bool_WSN_10x[counter][j] == 1:
                f.write(lines[(4*j)])
                f.write(lines[(4*j)+1])
                f.write(lines[(4*j)+2])
                f.write(lines[(4*j)+3])
        f.close()
        
        lane_counter = lane_counter + 1
        
        
        
    counter = counter + 1

#get all Read1 fastq file paths after starcode
all_R1_10x_mod_unfilt = glob.glob(filt_haveStart_10x + "/**/*_R1*.fastq", recursive = True)
all_R1_10x_mod_unfilt.sort()

# Remove any Read1 fastq files that you dont care about
all_R1_10x_mod_temp = [] 
for grp in GSAMP:
    for smp in grp:
        if str(smp) in str(all_R1_10x_mod_unfilt):
            all_R1_10x_mod_temp.append([s for s in all_R1_10x_mod_unfilt if str(smp) in s])

all_R1_10x_mod = [item for sublist in all_R1_10x_mod_temp for item in sublist]

#define samples
samples_R1_10x_mod = []
for paths in all_R1_10x_mod:
    print(paths)
    samples_R1_10x_mod.append(paths.split("/")[-1].split("_L0")[0])
samples_R1_10x_mod = list(set(samples_R1_10x_mod))
samples_R1_10x_mod.sort()

print(samples_R1_10x_mod)

# Use the booleans to write new files only using only good read2s in 10x samples
counter = 0
for sample in samples_R1_10x_mod:
    print(sample)
    s_fastq = []
    for path in all_R1_10x_mod:
        if "/"+ sample in path:
            s_fastq.append(path)
    
    # get all of the fastq lines from the modified fastqs
    lines = [] # advance counter as I read lines and only append those where the boolean is 1
    for fsmp in s_fastq:
        with open(fsmp) as f:
            for line in f:
                lines.append(line)

    lane_counter = 0
    for fsamp in s_fastq: 
        print(fsamp)
        #write files with these edited barcodes ( these are used as the input into starcode)
        f = open(filt_WSN_10x + "/" + fsamp.split('/')[-1],"w")
        print('Beginning count: ' + str(file_endpoints_10x_postsc[counter][lane_counter]))
        print('End count: '+ str(file_endpoints_10x_postsc[counter][lane_counter+1]))
        for j in range(file_endpoints_10x_postsc[counter][lane_counter],file_endpoints_10x_postsc[counter][lane_counter+1]):
            #add print statements
            if bool_WSN_10x[counter][j] == 1:
                f.write(lines[(4*j)])
                f.write(lines[(4*j)+1])
                f.write(lines[(4*j)+2])
                f.write(lines[(4*j)+3])
        f.close()
        
        lane_counter = lane_counter + 1
        
        
        
    counter = counter + 1

#get all Read1 fastq file paths after starcode
all_I1_10x_mod_unfilt = glob.glob(filt_haveStart_10x + "/**/*_I1*.fastq", recursive = True)
all_I1_10x_mod_unfilt.sort()

# Remove any Read1 fastq files that you dont care about
all_I1_10x_mod_temp = [] 
for grp in GSAMP:
    for smp in grp:
        if str(smp) in str(all_I1_10x_mod_unfilt):
            all_I1_10x_mod_temp.append([s for s in all_I1_10x_mod_unfilt if str(smp) in s])

all_I1_10x_mod = [item for sublist in all_I1_10x_mod_temp for item in sublist]

#define samples
samples_I1_10x_mod = []
for paths in all_I1_10x_mod:
    print(paths)
    samples_I1_10x_mod.append(paths.split("/")[-1].split("_L0")[0])
samples_I1_10x_mod = list(set(samples_I1_10x_mod))
samples_I1_10x_mod.sort()

print(samples_I1_10x_mod)

# Use the booleans to write new files only using only good read2s in 10x samples
counter = 0
for sample in samples_I1_10x_mod:
    print(sample)
    s_fastq = []
    for path in all_I1_10x_mod:
        if "/"+ sample in path:
            s_fastq.append(path)
    
    # get all of the fastq lines from the modified fastqs
    lines = [] # advance counter as I read lines and only append those where the boolean is 1
    for fsmp in s_fastq:
        with open(fsmp) as f:
            for line in f:
                lines.append(line)

    lane_counter = 0
    for fsamp in s_fastq: 
        print(fsamp)
        #write files with these edited barcodes ( these are used as the input into starcode)
        f = open(filt_WSN_10x + "/" + fsamp.split('/')[-1],"w")
        print('Beginning count: ' + str(file_endpoints_10x_postsc[counter][lane_counter]))
        print('End count: '+ str(file_endpoints_10x_postsc[counter][lane_counter+1]))
        for j in range(file_endpoints_10x_postsc[counter][lane_counter],file_endpoints_10x_postsc[counter][lane_counter+1]):
            #add print statements
            if bool_WSN_10x[counter][j] == 1:
                f.write(lines[(4*j)])
                f.write(lines[(4*j)+1])
                f.write(lines[(4*j)+2])
                f.write(lines[(4*j)+3])
        f.close()
        
        lane_counter = lane_counter + 1
        
        
        
    counter = counter + 1

#get all Read1 fastq file paths after starcode
all_I2_10x_mod_unfilt = glob.glob(filt_haveStart_10x + "/**/*_I2*.fastq", recursive = True)
all_I2_10x_mod_unfilt.sort()

# Remove any Read1 fastq files that you dont care about
all_I2_10x_mod_temp = [] 
for grp in GSAMP:
    for smp in grp:
        if str(smp) in str(all_I2_10x_mod_unfilt):
            all_I2_10x_mod_temp.append([s for s in all_I2_10x_mod_unfilt if str(smp) in s])

all_I2_10x_mod = [item for sublist in all_I2_10x_mod_temp for item in sublist]

#define samples
samples_I2_10x_mod = []
for paths in all_I2_10x_mod:
    print(paths)
    samples_I2_10x_mod.append(paths.split("/")[-1].split("_L0")[0])
samples_I2_10x_mod = list(set(samples_I2_10x_mod))
samples_I2_10x_mod.sort()

print(samples_I2_10x_mod)

# Use the booleans to write new files only using only good read2s in 10x samples
counter = 0
for sample in samples_I2_10x_mod:
    print(sample)
    s_fastq = []
    for path in all_I2_10x_mod:
        if "/"+ sample in path:
            s_fastq.append(path)
    
    # get all of the fastq lines from the modified fastqs
    lines = [] # advance counter as I read lines and only append those where the boolean is 1
    for fsmp in s_fastq:
        with open(fsmp) as f:
            for line in f:
                lines.append(line)

    lane_counter = 0
    for fsamp in s_fastq: 
        print(fsamp)
        #write files with these edited barcodes ( these are used as the input into starcode)
        f = open(filt_WSN_10x + "/" + fsamp.split('/')[-1],"w")
        print('Beginning count: ' + str(file_endpoints_10x_postsc[counter][lane_counter]))
        print('End count: '+ str(file_endpoints_10x_postsc[counter][lane_counter+1]))
        for j in range(file_endpoints_10x_postsc[counter][lane_counter],file_endpoints_10x_postsc[counter][lane_counter+1]):
            #add print statements
            if bool_WSN_10x[counter][j] == 1:
                f.write(lines[(4*j)])
                f.write(lines[(4*j)+1])
                f.write(lines[(4*j)+2])
                f.write(lines[(4*j)+3])
        f.close()
        
        lane_counter = lane_counter + 1
        
        
        
    counter = counter + 1

# If gDNA samples are used, filter
if barcodeSource == 'both':
    #get all Read1 fastq file paths after starcode
    all_R1_gDNA_mod_unfilt = glob.glob(mod_R2 + "/**/*_R1*.fastq", recursive = True)
    all_R1_gDNA_mod_unfilt.sort()

    # Remove any Read1 fastq files that you dont care about
    all_R1_gDNA_mod_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_R1_gDNA_mod_unfilt):
                all_R1_gDNA_mod_temp.append([s for s in all_R1_gDNA_mod_unfilt if str(smp) in s])

    all_R1_gDNA_mod = [item for sublist in all_R1_gDNA_mod_temp for item in sublist]

    #define samples
    samples_R1_gDNA_mod = []
    for paths in all_R1_gDNA_mod:
        print(paths)
        samples_R1_gDNA_mod.append(paths.split("/")[-1].split("_L0")[0])
    samples_R1_gDNA_mod = list(set(samples_R1_gDNA_mod))
    samples_R1_gDNA_mod.sort()

    print(samples_R1_gDNA_mod)
    
    # for gDNA samples, build a boolean of sequences that follow WSN nucleotide motif
    file_endpoints_gDNA_postsc = []
    bool_WSN_gDNA = []
    counter = 0
    for sample in samples_R1_gDNA_mod:
        print(sample)
        s_fastq = []
        for path in all_R1_gDNA_mod:
            if "/"+ sample in path:
                s_fastq.append(path)
            
        print(len(s_fastq))
    
        # Build boolean that sequences with runs of 4 or containing N are not used
        bool_WSN_gDNA.append([])
        file_endpoints_gDNA_postsc.append([])
        seqs = []
        for fsmp in s_fastq:
            for record in SeqIO.parse(fsmp, "fastq"):
                if 'AAAA' in record.seq:
                    bool_WSN_gDNA[counter].append(0)
                elif 'TTTT' in record.seq:
                    bool_WSN_gDNA[counter].append(0)
                elif 'GGGG' in record.seq:
                    bool_WSN_gDNA[counter].append(0)
                elif 'CCCC' in record.seq:
                    bool_WSN_gDNA[counter].append(0)
                elif 'NN' in record.seq:
                    bool_WSN_gDNA[counter].append(0)
                else:
                    bool_WSN_gDNA[counter].append(1)
                seqs.append(record.seq)
            print(len(seqs))
            file_endpoints_gDNA_postsc[counter].append(len(seqs))
        
        file_endpoints_gDNA_postsc[counter].insert(0,0)
        counter = counter+1
        
    # Use the booleans to write new files only using only good reads in gDNA samples
    counter = 0
    for sample in samples_R1_gDNA_mod:
        print(sample)
        s_fastq = []
        for path in all_R1_gDNA_mod:
            if "/"+ sample in path:
                s_fastq.append(path)
    
        # get all of the fastq lines from the modified fastqs
        lines = [] # advance counter as I read lines and only append those where the boolean is 1
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)

        lane_counter = 0
        for fsamp in s_fastq: 
            print(fsamp)
            #write files with these edited barcodes ( these are used as the input into starcode)
            f = open(filt_WSN_gDNA + "/" + fsamp.split('/')[-1],"w")
            print('Beginning count: ' + str(file_endpoints_gDNA_postsc[counter][lane_counter]))
            print('End count: '+ str(file_endpoints_gDNA_postsc[counter][lane_counter+1]))
            for j in range(file_endpoints_gDNA_postsc[counter][lane_counter],file_endpoints_gDNA_postsc[counter][lane_counter+1]):
                #add print statements
                if bool_WSN_gDNA[counter][j] == 1:
                    f.write(lines[(4*j)].strip())
                    f.write(str(Seq(lines[(4*j)+1]).reverse_complement())+"\n")
                    f.write(lines[(4*j)+2].strip())
                    f.write(lines[(4*j)+3][::-1]+"\n")
            f.close()
        
            lane_counter = lane_counter + 1
        counter = counter + 1
            
            
    #get all Index1 fastq file paths after starcode
    all_I1_gDNA_mod_unfilt = glob.glob(filt_haveStart_gDNA + "/**/*_I1*.fastq", recursive = True)
    all_I1_gDNA_mod_unfilt.sort()

    # Remove any Index1 fastq files that you dont care about
    all_I1_gDNA_mod_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_I1_gDNA_mod_unfilt):
                all_I1_gDNA_mod_temp.append([s for s in all_I1_gDNA_mod_unfilt if str(smp) in s])

    all_I1_gDNA_mod = [item for sublist in all_I1_gDNA_mod_temp for item in sublist]

    #define samples
    samples_I1_gDNA_mod = []
    for paths in all_I1_gDNA_mod:
        print(paths)
        samples_I1_gDNA_mod.append(paths.split("/")[-1].split("_L0")[0])
    samples_I1_gDNA_mod = list(set(samples_I1_gDNA_mod))
    samples_I1_gDNA_mod.sort()

    print(samples_I1_gDNA_mod)
    
    # Use the booleans to write new files only using only good indices in gDNA samples
    counter = 0
    for sample in samples_I1_gDNA_mod:
        print(sample)
        s_fastq = []
        for path in all_I1_gDNA_mod:
            if "/"+ sample in path:
                s_fastq.append(path)
    
        # get all of the fastq lines from the modified fastqs
        lines = [] # advance counter as I read lines and only append those where the boolean is 1
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)

        lane_counter = 0
        for fsamp in s_fastq: 
            print(fsamp)
            #write files with these edited barcodes ( these are used as the input into starcode)
            f = open(filt_WSN_gDNA + "/" + fsamp.split('/')[-1],"w")
            print('Beginning count: ' + str(file_endpoints_gDNA_postsc[counter][lane_counter]))
            print('End count: '+ str(file_endpoints_gDNA_postsc[counter][lane_counter+1]))
            for j in range(file_endpoints_gDNA_postsc[counter][lane_counter],file_endpoints_gDNA_postsc[counter][lane_counter+1]):
                #add print statements
                if bool_WSN_gDNA[counter][j] == 1:
                    f.write(lines[(4*j)])
                    f.write(lines[(4*j)+1])
                    f.write(lines[(4*j)+2])
                    f.write(lines[(4*j)+3])
            f.close()
        
            lane_counter = lane_counter + 1  
        counter = counter+1
    
    
    #get all Index1 fastq file paths after starcode
    all_I2_gDNA_mod_unfilt = glob.glob(filt_haveStart_gDNA + "/**/*_I2*.fastq", recursive = True)
    all_I2_gDNA_mod_unfilt.sort()

    # Remove any Index1 fastq files that you dont care about
    all_I2_gDNA_mod_temp = [] 
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_I2_gDNA_mod_unfilt):
                all_I2_gDNA_mod_temp.append([s for s in all_I2_gDNA_mod_unfilt if str(smp) in s])

    all_I2_gDNA_mod = [item for sublist in all_I2_gDNA_mod_temp for item in sublist]

    #define samples
    samples_I2_gDNA_mod = []
    for paths in all_I2_gDNA_mod:
        print(paths)
        samples_I2_gDNA_mod.append(paths.split("/")[-1].split("_L0")[0])
    samples_I2_gDNA_mod = list(set(samples_I2_gDNA_mod))
    samples_I2_gDNA_mod.sort()

    print(samples_I2_gDNA_mod)
    
    # Use the booleans to write new files only using only good indices in gDNA samples
    counter = 0
    for sample in samples_I2_gDNA_mod:
        print(sample)
        s_fastq = []
        for path in all_I2_gDNA_mod:
            if "/"+ sample in path:
                s_fastq.append(path)
    
        # get all of the fastq lines from the modified fastqs
        lines = [] # advance counter as I read lines and only append those where the boolean is 1
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)

        lane_counter = 0
        for fsamp in s_fastq: 
            print(fsamp)
            #write files with these edited barcodes ( these are used as the input into starcode)
            f = open(filt_WSN_gDNA + "/" + fsamp.split('/')[-1],"w")
            print('Beginning count: ' + str(file_endpoints_gDNA_postsc[counter][lane_counter]))
            print('End count: '+ str(file_endpoints_gDNA_postsc[counter][lane_counter+1]))
            for j in range(file_endpoints_gDNA_postsc[counter][lane_counter],file_endpoints_gDNA_postsc[counter][lane_counter+1]):
                #add print statements
                if bool_WSN_gDNA[counter][j] == 1:
                    f.write(lines[(4*j)])
                    f.write(lines[(4*j)+1])
                    f.write(lines[(4*j)+2])
                    f.write(lines[(4*j)+3])
            f.close()
        
            lane_counter = lane_counter + 1  
            
        counter = counter+1
        
# Filter the feature reference file to remove all of the barcodes that do not file WSN

#read in feature reference file and make a list of its lines
fr_file = CellR + "FeatureReference.csv"
print(fr_file)
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
    print(iters)

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

# Properly format the gDNA output files
# Currently is output in --seq-id format which means the read (by numerical order) that goes into 
#sequence is listed after

count = 1
for grp in GSAMP:
    #define sc output file for this group of samples
    for path in all_sc_out:
        if "/sc_output_group" + str(count) +"_comb.txt" in path:
            scfile = path
    print(scfile)
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
        print(smp)
        print(ind_to_test)
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