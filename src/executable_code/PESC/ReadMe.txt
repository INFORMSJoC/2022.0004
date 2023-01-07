1) Input parameters 
The program PBTSPESC contains following input parameters, where N is the number of unit spheres, 
NumberOfRuns is the number of times of running the algorithm,
and TimeLimit is the time limit (in seconds) for each run of algorithm:
————————————————————————————————————————
./PBTSPESC   N   NumberOfRuns  TimeLimit 
————————————————————————————————————————


2) Submission of jobs
The parameters of program can be located in a script file named 'N.sh', where N is the number of unit spheres. 
Under a linux operating system, the job 'N.sh' can be submitted as follows:
-----------------------------------------
chmod 777 PBTSPESC 
chmod 777 N.sh 
sbatch N.sh
-----------------------------------------
As an example, we provide a script file '72.sh'.

3) Computational results: 
The computational results are summarized in a text file named "CompuptationalResesults.txt", 
and for each instance the following information is given:
——————————————————————————————————————————————————————————————————————————
N     L_{best}      L_{avg}      L_{worst}      SR       Sigma    Time_avg
——————————————————————————————————————————————————————————————————————————
where N is the number of unit spheres, L_{best}, L_{avg} and L_{worst} represent respectively the best, average, and worst results over 'NumberOfRuns' runs . 
"SR" is the number of times that L_{best} is obtained, "Sigma" is the standard deviation of objective values obtained（L）, 
and "Time_avg" represents the average computational time (in seconds) for each run of algorithm.  


4) Solution file：
The best solution found is stored in a text file named "N.txt", where N is the number of unit spheres.
The format of solution file is as follows: 
The first line contains the number of unit spheres (N) and the size of cube container (L). 
As for the remaining parts, each line contains the X-coordinate, Y-coordinate and Z-coordinate of the center of one unit sphere. 
