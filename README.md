# shared_for_RGC
some small useful things to share with rgc

# FCup reading
You can run it using clas12root.

>source /group/clas12/packages/setup.csh
>
>module load clas12
>
>clas12root FCup_reading.cc("path_to_files",minimum run, maximum run)

- I am using the gmn train files to read the FCup information but please change line 9 as you need if you don't want to read them from here (the information is the same in all trains and recon). So as it is, the first argument is a path to where the gmn files are located. An example for it would be "/volatile/clas12/rg-c/production/dst/8.7.0_TBT/dst/train/gmn/".

- Second and last arguments are the run range you want to look at.

- Output is a text file with the FCup current in nC for each helicity state and each run.
