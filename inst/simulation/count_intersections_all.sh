#!/bin/bash                                                                                                                            

scenarios=(A B C D E F)
scenarios2000=(Acombined Bcombined Ccombined Dcombined Ecombined Fcombined)
scenarios10days=(A10days B10days C10days D10days E10days F10days)
nodes=(1:50 1:50 1:50)
max_nodes=(60 60 60)
days=(1 1 10)

test="no-test"
if [ "$1" == "test" ]
then
  test="test"
fi

function count_intersections_single {
  echo $test $1 $2 $3 $4
  if [ "$test" != "test" ]
  then
    ./count_intersections.sh ~/git/RParallelScreen/ $1 $2 $3 $4
  fi
}

function count_intersections_scenario {
  _scenarios=$1[@]
  scenarios=("${!_scenarios}")
  nodes=$2
  max_nodes=$3
  numdays=$4

  for scenario in "${scenarios[@]}"
  do
    count_intersections_single $scenario $nodes $max_nodes $model
  done
}

count_intersections_scenario scenarios "${nodes[0]}" "${max_nodes[0]}" "${days[0]}"
count_intersections_scenario scenarios2000 "${nodes[0]}" "${max_nodes[0]}" "${days[1]}"
count_intersections_scenario scenarios10days "${nodes[0]}" "${max_nodes[0]}" "${days[2]}"
