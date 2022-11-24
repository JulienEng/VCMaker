#! /bin/bash
############################################# 
###                                       ###    
###   VConverter script for Q-Chem vXXX   ###
###	FORCE output > VCM FORMAT         ### 
###					  ###
###     Julien Eng		 	  ###
#############################################



file=$1


## Determination of the number of atoms of the molecule.
counter=0
check=0
while [[ "${check}" == "${counter}" ]] ; do
	n_at=${check}
	counter=$(( ${counter} +1 ))
	check=`grep -A $(( ${counter} + 2 )) 'Standard Nuclear Orientation' $file | head -n $(( ${counter} +3 ))  | tail -n 1 | awk -F ' ' '{print $1}'`
done
########################################################

## Determination of the type of Gradient
if grep 'Full Analytical Gradient' $file > /dev/null 2>&1 ; then
	analytical=.TRUE.
	grepstr="Full Analytical Gradient"
elif grep 'Gradient of the state energy' $file > /dev/null 2>&1 ; then
	numerical=.TRUE.
	grepstr="Gradient of the state energy"
fi
########################################################

ncol=`grep -A 1 "${grepstr}" $file | tail -n 1 | awk -F ' ' '{print $NF}'`

size=$(( ${n_at} )) 			#Size of the Gradient (3 x size)
block=$(( (${n_at}+4)/${ncol} ))		#Number of blocks of 6columns the hessian is split in
floor=$(( (${n_at})/${ncol} ))		#FLOOOOOOR
rest=$(( ${n_at} - ${floor}*${ncol}))	#How many columns on the last block (if any)
l2grep=$(( ${block}*(4) ))	#Number of lines to grep.

###echo "####DEBUG####"
###echo "#Col: $ncol"
###echo "Size: $size"
###echo "Block: $block"
###echo "Floor: $floor"
###echo "Rest: $rest"
###echo "L2Grep: $l2grep"

	for (( j=1; j<=${floor}; j++ )) ; do
		for (( k=1 ; k<=${ncol} ; k++ )) ; do
			for (( i=1; i<=3; i++ )) ; do
		col=$(( (${j}-1)*${ncol} + ${k} ))
		array=`grep -A ${l2grep} "${grepstr}" $file | tail -n ${l2grep} | sed -n $(( (${j} -1)*(4)+ ${i}+1 ))p | awk -v col=$(( $k +1)) '{print $col}'`
		printf "%12.8f \n" ${array}
			done
		done
	done
	if [ ${rest} -gt 0 ] ; then
		for (( k=1 ; k<=${rest} ; k++ )) ; do
			for (( i=1; i<=3; i++ )) ; do
		col=$(( (${floor})*${ncol} + ${k} ))
		array=`grep -A ${l2grep} "${grepstr}" $file | tail -n ${l2grep} | sed -n $(( ${floor}*(4)+ ${i}+1 ))p | awk -v col=$(( $k +1)) '{print $col}'`
		printf "%12.8f \n" ${array}
			done
		done
	fi
