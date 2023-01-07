#!/bin/bash
#This is a script file to compile .cpp files

SRC=main.cpp
SRC1=CG_descent.cpp
SRC2=ObjGrad.cpp
SRC3=Variables.cpp
SRC4=Shake.cpp
SRC7=Threshold.cpp
EXE=PBTSPECS

srun g++ ${SRC} ${SRC1} ${SRC2} ${SRC3} ${SRC4} ${SRC7} -O3 -o ${EXE}


