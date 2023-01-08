For the PESC problem, under a linux operating system, please run the script file 'compile_pesc.sh' to generate the executable pragram 'PBTSPESC'. 

To solve the PESC problem with N unit spheres by using the pragram 'PBTSPESC', please submit the scribt file 'N.sh' as follows: 
-----------------------------------------
chmod 777 PBTSPESC
chmod 777 N.sh 
sbatch N.sh
-----------------------------------------

The program 'PBTSPESC' contains three input parameters which are listed in each script file 'N.sh', where N is the number of unit spheres, 
NumberOfRuns is the number of times of running the program,
and TimeLimit is the time limit (in seconds) for each run of program:
————————————————————————————————————————
./PBTSPESC   N   NumberOfRuns  TimeLimit 
————————————————————————————————————————