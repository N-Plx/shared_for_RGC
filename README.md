# shared_for_RGC
Useful things to share with rgc

# FCup reading
You can run it using clas12root. It can be run on ifarm with :

>source /group/clas12/packages/setup.csh
>
>module load clas12
>
>clas12root FCup_reading.cc("path_to_files",minimum run, maximum run)

- I am using the gmn train files to read the FCup information but please change line 9 as you need if you don't want to read them from here (the information is the same in all trains and recon). So as it is, the first argument is a path to where the gmn files are located. An example for it would be "/cache/clas12/rg-c/production/summer22/pass1/10.5gev/C/dst/train/gmn/".

- Second and last arguments are the run range you want to look at.

- Output is a text file with the FCup current in nC for each helicity state and each run.

# PbPt

This is the code to compute the elastic asymmetries for PbPt extraction of the polarized target.

### Step 1: Building the elastic events from RGC data.

> clas12root -l analysis_elastic(int run_number, string target_type)

Input: the run number you want to compute and the target type which can be “ND3”, “NH3”, “C”,… (it is only used to read in the proper folder for RGC data).

Outputs a ROOT file with two ttrees:
+ ”1electron” contains the elastic ep → e’p’ events and the corresponding exclusivity variables.
+ ”Scaler info” contains the accumulated Fcup charges for that run

This ROOT file is the input for step 2. 

### Step 2: Computing PbPt.

> root -l PbPt(string analysis_folder,string filename_C, string filename_signal, string target)

Input: 
+ string analysis_folder: folder containing the analysis files from step1.
+ string filename_C: text file containing a list of C runs you want to process (needed for dilution factor). Runs should be on one line, separated by commas.
+ string filename_signal: text file containing a list of signal runs you want to process. Runs should be on one line, separated by commas.
+ string target is the target type (it is only used to distinguish between files that contain the exclusivity cuts and are labeled depending on the target type)

Output: 
+ Plot of the asymmetries in Q2 bins: theoretical, measured for signal, measured for carbon and a ratio plot of the measured and theoretical asymmetry as a visual check.
+ Plot of the dilution factor in Q2 bins.
+ Text files containing in the first line the PbPt result and its error. Lines 2 to 4 are sanity checks for the beam charge asymmetry and carbon asymmetry. Following lines are all the asymmetries in Q2 bins.

Notes:
- The “Particle” class is used for handling particle kinematics/4-vectors.

Documentation for details of the analysis: [CLAS wiki](https://clasweb.jlab.org/wiki/index.php/Elastic_Analysis_for_PbPt_Extraction)

# Utils Folder
- utils for dealing with files
- cosmetics for plots
