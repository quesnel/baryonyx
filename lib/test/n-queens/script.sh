#!/bin/bash

> solutions.txt

for ((i = 0; i < 30; i++))
do
	cp n-queens-problem-$i.lp pb.lp
	cplex < n-queens-problem.commands
	mv cplex.log cplex-$i.log
	cat cplex-$i.log | grep ' Objective = ' | awk '{print $NF}' >> solutions.txt
done
