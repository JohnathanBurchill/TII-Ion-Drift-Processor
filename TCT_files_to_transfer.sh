#!/bin/bash

if test "$#" -ne "3"; then
  echo "Usage: $0 referenceDate version dataset"
  exit 1
fi

testDate=$(date -d $1 +%s)
version=$2
dataset=$3

a=`ls /databases/TCT/${version}/${dataset}/*.ZIP`

for i in $a; do 
  f=`basename $i`
  d=`echo $f | cut -c 20-27` 
  fileDate=$(date -d $d +%s)
  if [ $fileDate -ge $testDate ]; then 
    ln -s $i /databases/TCT/from_uoc/$version/$dataset/$f
    echo Linking $f
  fi
done 

