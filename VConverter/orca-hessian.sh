#! /bin/bash
############################################# 
###                                       ###    
###   VConverter script for Orca v5       ###
###	Frequency .hess output > VCM FORMAT   ### 
###					                      ###
###     Julien Eng		 	              ###
#############################################

file=$1

size=`grep -A 1 '$hessian' $1 | tail -n 1` 
block=$(( ${size}/5 ))		    #Number of blocks of 6columns the hessian is split in
rest=$(( ${size} - ${block}*5))	#How many columns on the last block (if any)

### D E B U G 
#echo "size: ${size}"
#echo "block: ${block}"
#echo "rest: ${rest}"
###

for (( lines=1 ; lines <=${size} ; lines++ )) ; do
    for (( k=1 ; k<=${block} ; k++ )) ; do
                line2grep=$(( 2+${lines}+(${k}-1)*(${size}+1) ))
        for col in {1..5} ; do
            colb=$(( ${col} + 1 ))
            counter=$(( (${k} -1)*5 + ${col} ))
            hess[${counter}]=`grep -A ${line2grep} '$hessian' $file | tail -n 1 | awk -v col="${colb}" -F ' ' '{print $col}'`
        done
    done
    if [[ ${rest} -ne 0 ]] ; then
    line2grep=$(( 2+${lines}+ ${block}*(${size}+1) ))
        for (( col=1 ; col<=${rest} ; col++ )) ; do
            colb=$(( ${col} + 1 ))
            counter=$(( ${block}*5 + ${col} ))
            hess[${counter}]=`grep -A ${line2grep} '$hessian' $file | tail -n 1 | awk -v col="${colb}" -F ' ' '{print $col}'`
        done
    fi

    for (( out=1 ; out<=${size} ; out++ )) ; do
       printf "%12.8f" ${hess[${out}]}
    done
        printf "\n"
done
