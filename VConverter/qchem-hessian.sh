# /bin/bash
############################################# 
###                                       ###    
###   VConverter script for Q-Chem vXXX   ###
###	Frequency output > VCM FORMAT     ### 
###       Requires VIBMAT_PRINT 6         ###
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
str2grep="Hessian of the SCF Energy"
chk=`grep "${str2grep}" $file`
if [[ ! ${chk} ]] ; then
str2grep="Final Hessian"
chk=`grep "${str2grep}" $file`
   if [[ ! ${chk} ]] ; then
   echo "No Hessian found in ${file}. Leaving now.."
   exit
   fi
fi

size=$(( 3*${n_at} )) 			#Size of the Hessian (size x size)
block=$(( (3*${n_at}+5)/6 ))		#Number of blocks of 6columns the hessian is split in
floor=$(( (3*${n_at})/6 ))		#FLOOOOOOR
rest=$(( ${n_at}*3 - ${floor}*6))	#How many columns on the last block (if any)
l2grep=$(( ${block}*(3*${n_at}+1) ))	#Number of lines to grep.

for (( i=1; i<=${size}; i++ )) ; do
	for (( j=1; j<=${floor}; j++ )) ; do
		for k in {1..6} ; do
		col=$(( (${j}-1)*6 + ${k} ))
		array=`grep -A ${l2grep} "${str2grep}" $file | tail -n ${l2grep} | sed -n $(( (${j} -1)*(3*${n_at}+1)+ ${i}+1 ))p | awk -v col=$(( $k +1)) '{print $col}'`
		printf "%12.8f" ${array}
		done
	done
	if [ ${rest} -gt 0 ] ; then
		for (( k=1 ; k<=${rest} ; k++ )) ; do
		col=$(( (${floor})*6 + ${k} ))
		array=`grep -A ${l2grep} "${str2grep}" $file | tail -n ${l2grep} | sed -n $(( ${floor}*(3*${n_at}+1)+ ${i}+1 ))p | awk -v col=$(( $k +1)) '{print $col}'`
		printf "%12.8f" ${array}
		done
	fi
	printf "\n"
done
