#!/bin/bash

if [ -z "$1" ]
then
echo "Missing parameter 'exec path'".
exit
fi
if [ -z "$2" ]
then
echo "Missing parameter 'scenario'".
exit
fi
if [ -z "$3" ]
then
echo "Missing parameter 'iterations'".
exit
fi
if [ -z "$4" ]
then
echo "Missing parameter 'max_nodes'".
exit
fi
if [ -z "$4" ]
then
echo "Missing parameter 'count_days'".
exit
fi

exec_path=$1
scenario=$2
iterations=$3
max_nodes=$4
count_days=$5

python "$exec_path"/parallel_r.py -t "$iterations" -n "$max_nodes" -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R notest "$scenario" "$count_days"

