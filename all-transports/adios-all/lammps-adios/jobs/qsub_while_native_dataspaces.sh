#!/bin/bash


      #  if(strcmp(transport_string, "ADIOS_DISK_MPIIO") == 0){
      #  else if (strcmp(transport_string, "ADIOS_STAGING_DSPACES") == 0){
      #  else if (strcmp(transport_string, "ADIOS_STAGING_DIMES") == 0){
      #  else if (strcmp(transport_string, "ADIOS_STAGING_FLEXPATH") == 0){
      #  else if (strcmp(transport_string, "NATIVE_STAGING_DSPACES") == 0){
      #  else if (strcmp(transport_string, "NATIVE_STAGING_DIMES") == 0){


trans=NATIVE_STAGING_DSPACES


napp1=32
napp2=16
napp3=2
nnode=4

napp1=512
napp2=256
napp3=32
nnode=64

#napp1=256
#napp2=128
#napp3=16
#nnode=32
#napp1=1024
#napp2=512
#napp3=64
#nnode=128


#napp2=8192
#napp1=16384
#napp3=1024
#nnode=2048


#if [[ "$trans" == *"STAGING"* ]];then
#    DS_SERVER=${DATASPACES_DIR}/bin/dataspaces_server
#    MPI_DS="aprun  -n ${napp1} -o 0 $DS_SERVER -s ${napp1} -c $((${napp2}+${napp1}))"
#	sleep 10s
#fi
#napp1=4096
#napp2=2048
#napp3=256
#nnode=512
#napp2=4096
#napp1=8192
#napp3=512
#nnode=1024


#16384
while [ $napp1 -le 8192 ]
do
  	echo "Welcome $nnode, $napp1, $napp2, $jobid times"
	sed "s/\${nnode}/${nnode}/g;s/\${napp1}/${napp1}/g;s/\${napp2}/${napp2}/g;s/\${trans}/${trans}/g;s/\${napp3}/${napp3}/g" lammps_adios_titan_native_dataspaces_var.job  > testnewjobv2
  	jobid=`qsub   testnewjobv2`
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

  	napp1=$(( $napp1 * 2 ))
  	napp2=$(( $napp2 * 2 ))
	napp3=$(( $napp3 * 2 ))
  	nnode=$(( $nnode * 2 ))  
done


#sed  's/${nnode}/3/g;s/${napp1}/2/g'  lammps_adios_titan_var.job  > testnewjob
