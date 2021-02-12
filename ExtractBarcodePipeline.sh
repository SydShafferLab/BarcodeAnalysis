#/bin/bash

#to run this simply log on to the consign PMACS server and run this line (replacing PATH with the path to this file): bsub -e Error.e -o Output.o sh /PATH/ExtractBarcodePipeline.sh

bsub -Is bash

########## Fill out these parameters NOTE: this only specifies some of the parameters you can chage look at TimeMachine Wiki to see further specification on bens python commands

#path to "BarcodeAnalysis" file downloaded from github
scripts="/home/gharm/BarcodeAnalysis"

# path to project folder (this should be a directory you make containing your Fastq files and your stagger file. Outputs will also go in this file)
project_folder="/project/shafferslab/Guillaume/ToddBC2"

# path to fastq files
fastq_path="/project/shafferslab/Guillaume/ToddBC2/FASTQ"

#path to staggerfile
##the stagger file is a csv file containing 2 columns, the first column contains the sample name and the second contains the number of stagger bases in the i5 primer used for that sample
staggerFile="/project/shafferslab/Guillaume/ToddBC2/staggerfile.csv"

#If you got FASTQ  files from BaseSpace write "bs" if you got them by running bcl2fastq write "bf"
fastq_origin="bs"

#if you want the program to only look at 3' end of the barcode enter "before" if you want it to look for it at both the 3'and 5' end of the barcode enter "both"
check_vector="both"

#input how many bases of barcode is sequenced (for default enter 100)
length="75"

#specify Levenshtein distance (number of mismatches accepatable to consider barcodes to be the same)(the default is "8")
ld="8"

#Decide whether final barcode count is based on number of time a barcode is seen (input: "Reads") or if you want to take into account the UMI (input: "UMI")
countType="Reads"

######## End of user modulated parameters

##build out paths based on user inputs
#source $scripts/bcenv/bin/activate
source $scripts/bcenv/bin/activate
PATH=$PATH:${scripts}/starcode/
TM_scripts=$scripts/PythonFunctions
#put this path in environment
export TM_scripts=$TM_scripts

##Start Runing code


echo 1
#organize fastq files properly
if [[ $fastq_origin = "bs" ]]; then
  #Ben script step #0
  mkdir $project_folder/repo/
  mkdir $project_folder/repo/raw/
  cp -r $fastq_path/. $project_folder/repo/raw/
  rm $project_folder/repo/raw/*.json
  chmod -R 775 $project_folder/repo/
  rm -r $project_folder/repo/raw/boot
  rm -r $project_folder/repo/raw/dev
  rm -r $project_folder/repo/raw/proc
  python $TM_scripts/stepZeroReorganizeBasespaceFiles/reorganizeBaseSpaceOutput.py $project_folder/repo/raw/
else
  if [[ $fastq_origin = "bf" ]]; then
    #get number of fastq files
    num_files=$(ls -l $fastq_path | wc -l)
    for (( i=1 ; ((i-$num_files-1)) ; i++ ))
    do
      mkdir $project_folder/repo/raw/${i}/
      find $fastq_path -name "*Guillaume_Sample_${i}_*" -exec cp '{}' $project_folder/repo/raw/${i}/ \;
    done;
  else
    echo "ERROR: source of FASTQ file improperly defined"
  fi
fi

# at this point you should have all fastq files in a folder called "raw" in your project folder. In "raw folder each sample should have a folder with the sample name and the FASTQ files corresponding to that sample in the project_folder

##Extract the Barcodes from the fastq sample (Ben script step #1)

python $TM_scripts/stepOneExtractBarcodes/allExtractBarcodes.py $project_folder/repo/ --staggerFile $staggerFile -r -l $length --check_vector $check_vector


##Identify barcodes that are similar (differences likely due to PCR errors or sequencing) and count them as 1 barcode (Ben script step #2)

python $TM_scripts/stepTwoRunStarcode/allSubmitStarcode.py $project_folder/repo/ -d $ld --check_vector $check_vector -c $countType
