#!/bin/sh

let concurrent=1000
#lcg-infosites --vo gridpp wms | grep ac.uk
export GLITE_WMS_WMPROXY_ENDPOINT=https://lcgwms05.gridpp.rl.ac.uk:7443/glite_wms_wmproxy_server

[ -e jobs.txt ] || touch jobs.txt
mapfile -t jobs < jobs.txt

while : ; do

    trap 'echo "unsafe to abort, wait for sleep"' 2

    echo "checking existing jobs"
    for job in "${jobs[@]}"
    do
	status=`glite-wms-job-status $job | grep Status: | cut -f 2 | cut -d: -f 2 | sed 's/\s//g'`
	case "$status" in 
	    Done*)
		echo "retrieving $job"
		glite-wms-job-output --dir /data/robinson/musun/grid $job
		;&
	    Cancelled|Aborted)
		echo $job $status >> complete.txt
		;&
	    "")
		jobs=(${jobs[@]/"$job"})
		;;
	esac
    done

    [ -e runjobs.stop ] || echo "topping up running jobs"
    while [ ${#jobs[@]} -lt $concurrent ] && [ ! -e runjobs.stop ]; do
	job=`glite-wms-job-submit -a test.jdl | grep https | grep 9000`
	echo "job $job submitted"
	if [ -z "$job" ]; then
	    echo "job submission failed"
	else
	    jobs+=($job)
	fi
    done
    for job in "${jobs[@]}"
    do
	echo "$job"
    done > jobs.txt

    if [ ${#jobs[@]} -eq 0 ]; then
	exit
    fi

    trap 2

    [ -e runjobs.abort ] && exit

    echo "sleep 10 minutes"
    sleep 60

done


