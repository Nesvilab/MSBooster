# MSBooster
MSBooster allows users to add deep learning-based features to .pin files before Percolator PSM rescoring

## System requirements
### All software dependencies and OS
• FragPipe, either version 17 or 18 (Download at https://github.com/Nesvilab/FragPipe/releases)
### Versions the software has been tested on
•	Windows

	• v10/11
	• Java(TM) SE Runtime Environment (build 16.0.1+9-24)
• Linux, Architecture: amd64

	• Java 11.0.16.1, OpenJDK 64-Bit Server VM, Red Hat, Inc.
      
## Installation guide
### Instructions
•	If using MSBooster in FragPipe, follow the instructions from the download link above. In this example, we use FragPipe-jre-18.0.zip. Extract it from the zip file and place the FragPipe-fre-18.0 folder in the “software” folder

•	If running standalone MSBooster, the instructions above to download FragPipe must still be followed. While standalone MSBooster does not need FragPipe, it requires the DIA-NN prediction tool in this tutorial which can be found at fragpipe/tools/diann/1.8.1/win/DiaNN.exe
 
### Typical install time
•	FragPipe: Less than 2 minutes
  
## Demo
**NOTE**: The 2 mzML files and MSBooster standalone jar file are at https://www.dropbox.com/sh/zr69smrov3qhgm2/AACgiJPMjp2erSMgyBqyBBmQa?dl=0. Please download and place in the “software” folder

### Instructions to run on demo data (Windows)
#### FragPipe GUI
•	If running in FragPipe, MSBooster is already available and loaded

•	In the Workflow tab, select the Default workflow and click “Load workflow”
 
•	Load the two mzML files under the Workflow tab using “Add files”
 
•	Pin files from the MSFragger search are already available in this demo folder, so uncheck “Run MSFragger” in the MSFragger tab
 
•	Although not used by MSBooster, FragPipe requires a fasta file. It is provided in the “software” folder (2022-03-18-decoys-reviewed-contam-UP000005640.fas)
 
•	In the Validation tab, uncheck all boxes except for “Run Validation Tools” and “Run MSBooster”. To uncheck other tools such as ProteinProphet, “Run Validation Tools” may need to be unchecked and checked again to allow interaction with other checkboxes. “Predict RT” and “Predict spectra” are also checked by default, but either one can be turned off. “Use correlated features” is unchecked by default, but is demonstrated in the manuscript and you can also test it
 
•	In the Run tab, set “Output dir” to the location of your demo folder (contains mzml, pin, pepxml files; fragger.params). Then click “RUN”
 
•	On a normal desktop (Intel(R) Core(TM) i7-8700 CPU @ 3.20GHz   3.19 GHz, 32.0 GB), MSBooster in FragPipe took 0.3 minutes to annotate the 2 pin files provided

#### MSBooster command line
•	Make sure the MSBooster-1.1.6-jar-with-dependencies.jar file and two mzML files are in the “software” folder, or at least the same folder as the two pin files (download link at beginning of Demo section). The FragPipe-jre-18.0 folder should be in the same folder so that DiaNN.exe is available for predictions

•	Open the command line by typing “cmd” in the Windows search bar

•	Navigate to the folder of this demo with the “cd” command

•	Optional: open “msbooster_params_template.txt”. There are various parameters to change, such as file locations
 
•	Enter “java -jar MSBooster-1.1.6-jar-with-dependencies.jar --paramsList msbooster_params_template.txt”

•	In the parameter file, parameters can be customized

	•	numThreads ishow many threads you want to use on the computer to run MSBooster (default: available threads minus 1)
	•	renamePin = 1 ensures the original pin files are not overwritten; instead, new ones with the “_edited” suffix are added. Setting this to 0 will rewrite the old pin files
	•	mzmlDirectory is where the mzml files are stored. If running MSBooster from the demo folder, “.” designates they are in the same folder
	•	pinPepXMLDirectory is where the pin/pepxml files are stored. This could be set as a folder name, or individual pin file name separated by spaces. While on the first run of MSBooster in the demo this could be set as “.” like mzmlDirectory is, it will produced “_edited” pin files. Upon running MSBooster a second time, MSBooster will search for mzml files that include “_edited” in the name and will fail because they are nonexistent. Therefore, it is safer to have the pin files listed out.
	•	If information on other parameters is needed, run “java -jar MSBooster-1.1.6-jar-with-dependencies.jar --help”
	
## Expected output
### MSBooster
•	MSBooster-annotated pin files with the suffixed “_edited.pin”

•	spectraRT.tsv and spectraRT_full.tsv: input files for DIA-NN prediction model

•	spectraRT.predicted.bin: predictions from DIA-NN to be used by MSBooster for feature calculation. Not in human readable format
	
### FragPipe-specific output
•	Log*.txt: whatever is printed by FragPipe when run 

•	Lcms-files*.fp-manifest: lists which mzml files were used

•	msbooster_params.txt: input parameters for MSBooster

### Expected run time for “demo”
•	Less than 30 seconds

## Instructions for use
### How to run the software on your data
•	If running through FragPipe, MSBooster is automatically given the correct parameters to work within the pipeline. Please refer to https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html 

•	If running on the command line, change the parameter file (in the demo) to use new file locations or customizable parameters. They are listed and described in the –help section. They can also be individually written on the command line, but it may be easier to organize by writing them in the parameter file and passing it with “—paramsList”.
