# Script will take in gDNA outputs from sc_outputs and analyze them


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

import jstyleson
import glob
from matplotlib_venn import venn3
import json
import warnings

import pandas as pd

#---------- User ---------------
path_to_feature_ref = "/Volumes/GoogleDrive/My Drive/Box_20220118/Projects/Proliferative_Senesent/Data/Split_barcode/starcode_outputs/FeatureReference_filtered.csv"

path_to_sc_output_folder = "/Volumes/GoogleDrive/My Drive/Box_20220118/Projects/Proliferative_Senesent/Data/Split_barcode/starcode_outputs"

output_folder = "/Users/raul/Desktop/"

strtseq = "GCTGTACAAGTAGGAT"

GSAMP = [["RAR2","RAR2_1","RAR2_2"],["RAR3","RAR3_1","RAR3_2"]]


#------------------------------



print("Running ")
print(" ")


#get all gDNA Read1 fastq file paths
all_sc_output = glob.glob(path_to_sc_output_folder + "/sc_output_counts_*.txt", recursive = True)

counter = 0
for grp in GSAMP:

    # Remove any files not in grp
    subset_sc_output = []
    for smp in grp:

    	hold_path = []
    	for path in all_sc_output:
    		if smp+"_S" in path:                #is this true
    			subset_sc_output.append(path)


    #------ Lets make the custom FeatureRef ------
    df_ref = pd.read_csv(path_to_feature_ref)
    ind_dict_ref = dict((k,i) for i,k in enumerate(df_ref['sequence']))

    for path in subset_sc_output:

    	df_output = pd.read_csv(path,header = None, sep='\t', engine='python')

    	df_output[0] = df_output [0].str[len(strtseq):]

    	ind_dict_output = dict((k,i) for i,k in enumerate(df_output[0]))
    	inter = set(df_ref['sequence']).intersection(df_output[0])
    	indices_ref = [ ind_dict_ref[x] for x in inter ]
    	indices_out = [ ind_dict_output[x] for x in inter ]

    	hold_zero = [0]*len(df_ref['sequence'])
    	for x in inter:
    		hold_zero[ind_dict_ref[x]] = df_output[1][ind_dict_output[x] ] 

    	df_ref[path.split('/')[-1]] = hold_zero

    counter += 1
    filepath = output_folder + "/" + '/FeatureReference_filtered_group'+ str(counter) + '_counts.csv'
    df_ref.to_csv(filepath)  



print("     Done :D ")
print(" ")

#---------plot-----------

# plot_f = plt.figure()
# set1 = set(seq_base_run_no)
# set2 = set(fuzzy_no[cnt_i])
# set3 = set(qscore_no)

# ax = plt.gca() 

# total = len(set1.union(set2,set3))

# venn3([set1, set2,set3], ('Base_repeat', 'Match_to_StartSeq','Qscore'),subset_label_formatter=lambda x: f"{(x/total):1.0%}");

# plt.suptitle(sample + "_venn" +  "   pct_startseq_match=" + str(fuzz_i), y=0.99);
# plt.title("Percent of total reads: Qscore=" + str("{:3.0f}".format(100*len(qscore_no)/len(seqs))) + " Base_repeat=" + str("{:3.0f}".format(100*len(seq_base_run_no)/len(seqs))) + " Match_to_StartSeq=" + str("{:3.0f}".format(100*len(fuzzy_no[cnt_i])/len(seqs))) );
# pdf.savefig(plot_f);
# plt.close()
