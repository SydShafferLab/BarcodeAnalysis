#/bin/bash

#to run this simply log on to the consign PMACS server and run this line (replacing PATH with the path to this file): bsub -e Error.e -o Output.o sh /PATH/ExtractBarcodePipeline.sh

bsub -Is bash
source $HOME/my_python/2.7.5/bin/activate
PATH=$PATH:/home/gharm/starcode/

########## Fill out these parameters NOTE: this only specifies some of the parameters you can chage look at TimeMachine Wiki to see further specification on bens python commands

#If you got FASTQ  files from BaseSpace write "bs" if you got them by running bcl2fastq write "bf"
fastq_origin="bf"

# path to fastq files
fasq_path="/project/shafferslab/Guillaume/Todd_Barcode/RawData/FirstExperiment/seq2/FASTQ/"

# path to project folder
project_folder="/project/shafferslab/Guillaume/Todd_Barcode"

#path to Ben's TimeMachine scripts
TM_sctipts="/BarcodeAnalysis/PythonFunctions/"

#path to staggerfile
##the stagger file is a csv file containing 2 columns, the first column contins the sample name and the second contains the number of stagger bases in the i5 primer used for that sample
staggerFile="/project/shafferslab/Guillaume/Todd_Barcode/repo/Replicate1seq1/staggerfile.csv"

#if you want the program to only look at 3' end enter "before" if you want it to look for it at both the 3'and 5' end enter "both"
check_vector="before"

#input how many bases of barcode is sequenced (for default enter 100)
length="25"

#specify Levenshtein distance (number of mismatches accepatable to consider barcodes to be the same)(the default is "8")
ld="4"

#Decide whether final barcode count is just based on number of time a barcode is seen (input: "Reads") or if you want to take into account the UMI (input: "UMI")
countType="Reads"

######## End of user modulated parameters
echo 1
#organize fastq files properly
if [[ $fastq_origin = "bs" ]]; then
  #Ben script step #0
  mkdir $project_folder/repo/raw/
  cd $project_folder/raw/
  python $TM_sctipts/stepZeroReorganizeBasespaceFiles/reorganizeBaseSpaceOutput.py $fasq_path
else
  if [[ $fastq_origin = "bf" ]]; then
    #get number of fastq files
    num_files=$(ls -l $fasq_path | wc -l)
    for (( i=1 ; ((i-$num_files-1)) ; i++ ))
    do
      mkdir $project_folder/repo/raw/${i}/
      find $fasq_path -name "*Guillaume_Sample_${i}_*" -exec cp '{}' $project_folder/repo/raw/${i}/ \;
    done;
  else
    echo "ERROR: source of FASTQ file improperly defined"
  fi
fi

# at this point you should have all fastq files in a folder called "raw" in your project folder. In "raw folder each sample should have a folder with the sample name and the FASTQ files corresponding to that sample in the project_folder

##Extract the Barcodes from the FASQ sample (Ben script step #1)

python $TM_sctipts/stepOneExtractBarcodes/allExtractBarcodes.py $project_folder/repo/ --staggerFile $staggerFilev -r -l $length --check_vector $check_vector


##Identify barcodes that are similar (differences likely due to PCR errors or sequencing) and count them as 1 barcode (Ben script step #2)

python $TM_sctipts/stepTwoRunStarcode/allSubmitStarcode.py $project_folder/repo/Replicate1seq1/ -d $ld --check_vector $check_vector -c $countType
