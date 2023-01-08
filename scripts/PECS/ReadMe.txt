For the PECS problem, under a linux operating system, please run the script file 'compile_pecs.sh' to generate the executable pragram 'PBTSPECS'. 

To solve the PECS problem with N circles by using the pragram 'PBTSPECS', please submit the scribt file 'N.sh' as follows: 
-----------------------------------------
chmod 777 PBTSPECS
chmod 777 N.sh 
sbatch N.sh
-----------------------------------------

The program PBTSPECS contains three input parameters which are listed in each script file 'N.sh', where N is the number of unit circles, 
NumberOfRuns is the number of times of running the program,
and TimeLimit is the time limit (in seconds) for each run of program:
————————————————————————————————————————
./PBTSPECS   N   NumberOfRuns  TimeLimit 
————————————————————————————————————————