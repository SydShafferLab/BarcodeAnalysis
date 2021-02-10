# BarcodeAnalysis

NOTE: This repository is just mean to link existing packages (that I did not write) to make identifying barcodes as easy as possible. All this code int the "PythonFuntions" folder is written by Ben Emert and the original code as well as it's documentation can be found here:https://github.com/arjunrajlaboratory/timemachine. I have also included the package that Ben's code relies on (starcode) so that there is only one download step. Starcode was written by Guillaume Filion et. al and the original code as well as its documentation can be found here:https://github.com/gui11aume/starcode.



## Setting up the environment (you should only need to do this once)
   All this code is run on the pmacs cluster as starcode is very memory intensive.
 1. Get if you dont already have it I recoment you download cyberduck (https://cyberduck.io) to easily interface with the cluster
 2. Once cyberduck is installed, open it, hit open connection, in the top drop down menu select "SFTP(SSH File Transfer Protocol)", for server type in "mercury.pmacs.upenn.edu", enter your PMACS username and password, and hit connect.
 3. Download all the files in this github repository by hitting the green "Code" button on the top right of this page
 4. Drag and drop the downloaded file "BarcodeAnalysis-main" into cyberduck (this should cause it to download on your cluster home directory)


## Preparing files
  To run barcode analysis you need to put your FASTQ files generated from sequencing your barcodes as well as a file we call a "staggerfile" on the cluster
  1. Here is the explination of the stagger file taken directly from the wiki written by Ben found here (https://github.com/arjunrajlaboratory/timemachine)
  
"If using our original Rewind vector and primers for library prep, then your sequencing templates will all contain the same sequence of DNA used for PCR priming. Unfortunately, on Illumina instruments, if all sequencing templates contain the same nucleotides at the same positions, this can create errors during sequencing or fail the run entirely. To mitigate this issue, each primer contains a different number of extra bases to stagger the primer region during sequencing (see timeMachineVectorPrimers.csv)
Based on which primers were used, we recommend specifying the number of "stagger" bases expected for each sample in a stagger file. Otherwise, the default is to assume 0 stagger bases which may lead to reads containing true barcodes being filtered out and not counted."

  2. to make this file simply have the first colum as sample names, and the second column as the stagger number for that sample and save this file as a .csv. The file "staggerfile.csv" included in this repository is an example of this file.
  3. For your fastq files, you should be able to download them from base space or get them as an output from bcl2fastq
  4. Using cyberduck create a folder for this specific project on the cluster
  5. to upload the stagger file you can simply drag and drop the file into the folder
  6. to upload the fastq you can try and drag and drop, but this is usually too slow. A better way for large files is to write this command in terminal `rsync -rav /path/containing/fast1/files/ <username>@mercury.pmacs.upenn.edu:/path/to/project/folder` replacing with correct paths and your username. it will then prompt you for your PMACS password and allow you to start transfering
  7. In the same project folder add a copy of the "ExtractBarcodePipeline.sh" file
  8. open the copied "ExtractBarcodePipeline.sh" file by right clicking on it and selecting "edit with..." and select a program with which to edit the file (text edit is fine, but I prefer Atom:https://atom.io).
  9. In this file edit the paths and parameters indicated under the line "########## Fill out these parameters ". for most applications you should not need to touch anyting under the line "######## End of user modulated parameters"
  10. once you have adjusted all the parameters, save the file
  
  ## Running barcode detection
 1. Open a terminal window and sigh into your pmacs by running the line `ssh <username>@consign.pmacs.upenn.edu` (replacing <username> with your PMACS username) hit enter, and then enter your PMACS password and hit enter
 2. Now activate a node on the cluster by running the command `bsub -Is bash`
 3. Next activate the line `bsub -e Error.e -o Output.o sh /PATH/ExtractBarcodePipeline.sh` replacing "PATH" with the correct directory to your edited ExtractBarcodePipeline.sh file wich should be in your project folder



