#!/bin/bash


      #  if(strcmp(transport_string, "ADIOS_DISK_MPIIO") == 0){
      #  else if (strcmp(transport_string, "ADIOS_STAGING_DSPACES") == 0){
      #  else if (strcmp(transport_string, "ADIOS_STAGING_DIMES") == 0){
      #  else if (strcmp(transport_string, "ADIOS_STAGING_FLEXPATH") == 0){
      #  else if (strcmp(transport_string, "NATIVE_STAGING_DSPACES") == 0){
      #  else if (strcmp(transport_string, "NATIVE_STAGING_DIMES") == 0){


trans=ADIOS_STAGING_DSPACES


napp1=32
napp2=16
napp3=2
nnode=4

napp1=512
napp2=256
napp3=32
nnode=64

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
while [ $napp1 -le 8192 ]
do
  	echo "Welcome $nnode, $napp1, $napp2, $jobid times"
	sed "s/\${nnode}/${nnode}/g;s/\${napp1}/${napp1}/g;s/\${napp2}/${napp2}/g;s/\${trans}/${trans}/g;s/\${napp3}/${napp3}/g" lammps_adios_titan_dataspaces_var.job  > testnewjob
  	jobid=`qsub   testnewjob`
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
