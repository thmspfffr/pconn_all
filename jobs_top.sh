#!/bin/sh

# Submits n jobs to the torque queing system

for i in {1..50}
do
  let var1=40*$i;
  echo 'Start Job ' $i 'wait for: ' $var1 's'
  qsub -v var="$var1" jobs.sh
done