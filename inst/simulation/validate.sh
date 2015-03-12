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
if [ -z "$5" ]
then
  echo "Missing parameter 'model'".
  exit
fi
if [ -z "$6" ]
then
  echo "Missing parameter 'free_mem'".
  exit
fi

exec_path=$1
scenario=$2
iterations=$3
max_nodes=$4
model=$5
free_mem=$6

if [ "$model" != "FMPModel" ]
then
  python "$exec_path"/parallel_r.py -p 15 -t "$iterations" -n "$max_nodes" -l 20.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/validate.R notest "$scenario" "$model"
fi
