#! /bin/bash

grad=$1

line_init=`grep -n "The current gradient" $1 | awk -F ':' '{print $1}'` 
line_init=$(( ${line_init} +2 ))
line_final=`grep -n "The atomic numbers and current" $1 | awk -F ':' '{print $1}'` 
line_final=$(( ${line_final} -2 ))
line_end=$(( ${line_final} +1 ))


sed -n "${line_init},${line_final}p;${line_end}q" $1
