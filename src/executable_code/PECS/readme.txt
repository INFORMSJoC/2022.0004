1) Input parameters 
The program PBTSPECS contains following input parameters, where N is the number of unit circles, 
NumberOfRuns is the number of times of running the algorithm,
and TimeLimit is the time limit (in seconds) for each run of algorithm:
————————————————————————————————————————
./PBTSPECS   N   NumberOfRuns  TimeLimit 
————————————————————————————————————————


2) Submission of jobs
The parameters of program can be located in a script file named 'N.sh', where N is the number of unit circles. 
Under a linux operating system, the job 'N.sh' can be submitted as follows:
-----------------------------------------
chmod 777 PBTSPECS
chmod 777 N.sh 
sbatch N.sh
-----------------------------------------


3) Computational results: 
The computational results are summarized in a text file named "CompuptationalResesults.txt", 
and for each instance the following information is given:
——————————————————————————————————————————————————————————————————————————
N     L_{best}      L_{avg}      L_{worst}      SR       Sigma    Time_avg
——————————————————————————————————————————————————————————————————————————
where N is the number of unit circles, L_{best}, L_{avg} and L_{worst} represent respectively the best, average, and worst results over 'NumberOfRuns' runs. 
"SR" is the number of times that L_{best} is obtained, "Sigma" is the standard deviation of objective values obtained（L）, 
and "Time_avg" represents the average computational time (in seconds) for each run of algorithm.   


4) Solution file：
The best solution found is stored in a text file named "N.txt", where N is the number of unit circles.
The format of solution file is as follows: 
The first line contains the number of unit circles (N) and the size of square container (L). 
As for the remaining parts, each line contains the X-coordinate and Y-coordinate of the center of one unit circle. 
