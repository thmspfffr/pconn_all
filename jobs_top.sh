#!/bin/sh

# Submits n jobs to the torque queing system

for i in {1..50}
do
  echo 'Start Job rev' $i
  qsub jobs.sh
  sleep 15
done

