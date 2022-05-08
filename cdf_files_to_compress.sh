#!/bin/bash

# checks number of TCT16 or TCT02 cdf files ot compress for specified satellite and version 
# http://stackoverflow.com/questions/28226229/how-to-loop-through-dates-using-bash

satellite=$1
version=$2
dataset=$3
numCpus=$4

if test "$#" -ne "4"; then
  echo "Usage: $0 satelliteLetter version dataset numCpus"
  exit 1
fi

compressionDir="/databases/TCT/${version}/${dataset}"
echo ${compressionDir}

filesToCompress=`find ${compressionDir} -name "SW_EXPT_EFI${satellite}_*.cdf"`
fileArray=()
for i in ${filesToCompress}; do
  ba=`basename $i .cdf`
  if [ ! -f ${compressionDir}/${ba}.ZIP ]; then
    fileArray+=("${i}")
  fi
done
numberToCompress="${#fileArray[@]}"
echo Number to compress: "${numberToCompress}"

exit 0



