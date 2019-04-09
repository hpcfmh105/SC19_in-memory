#!/bin/bash


      #  if(strcmp(transport_string, "ADIOS_DISK_MPIIO") == 0){
      #  else if (strcmp(transport_string, "ADIOS_STAGING_DSPACES") == 0){
      #  else if (strcmp(transport_string, "ADIOS_STAGING_DIMES") == 0){
      #  else if (strcmp(transport_string, "ADIOS_STAGING_FLEXPATH") == 0){
      #  else if (strcmp(transport_string, "NATIVE_STAGING_DSPACES") == 0){
      #  else if (strcmp(transport_string, "NATIVE_STAGING_DIMES") == 0){


trans=ADIOS_STAGING_DSPACES


napp1=1024
napp2=512
napp3=64
nnode=128

#napp1=64
#napp2=32
#napp3=4
#nnode=8

#napp1=128
#napp2=64
#napp3=8
#nnode=16

#napp2=8192
#napp1=16384
#napp3=1024
#nnode=2048

#napp2=1024
#napp1=2048
#napp3=512
#nnode=258

#if [[ "$trans" == *"STAGING"* ]];then
#    DS_SERVER=${DATASPACES_DIR}/bin/dataspaces_server
#    MPI_DS="aprun  -n ${napp1} -o 0 $DS_SERVER -s ${napp1} -c $((${napp2}+${napp1}))"
#	sleep 10s
#fi


#16384
ss=1
dim1=256
dim2=256
#16384
while [ $dim1 -le 512 ]
do
  	echo "Welcome $nnode, $napp1, $napp2,$dim1, $dim2, $jobid times"
	sed "s/\${nnode}/${nnode}/g;s/\${napp1}/${napp1}/g;s/\${napp2}/${napp2}/g;s/\${trans}/${trans}/g;s/\${napp3}/${napp3}/g;s/\${dim1}/${dim1}/g;s/\${dim2}/${dim2}/g" laplace_titan_only.job  > testnewjob
  	jobid=`qsub -q debug  testnewjob`
  	echo "jobid: $jobid" 
	rc=0
    sleep 15s

    while [ $rc -eq 0 ]
    do
        check=`showq -u dhuang`
        echo $check | grep "Idle"  > /dev/null
        rc="$?"
        echo "rc: idle $jobid $rc"
        sleep 3s
    done
	rc=0
    sleep 15s
    while [ $rc -eq 0 ]
    do
        check=`showq -u dhuang`
        echo $check | grep "Running"  > /dev/null
        rc="$?"
        echo "rc: running $jobid $rc"
        sleep 3s
    done

	qdel $jobid

	echo "delete job: $jobid"
	sleep 5s
    if [ $(( $ss % 2 )) -eq 0 ]; then
    dim2=$(( $dim2 * 2 ))
    else
    dim1=$(( $dim1 * 2 ))
    fi
    ss=$(( $ss + 1 ))
done


#sed  's/${nnode}/3/g;s/${napp1}/2/g'  lammps_adios_titan_var.job  > testnewjob
