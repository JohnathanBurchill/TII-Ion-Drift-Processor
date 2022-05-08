#!/bin/bash

# compresses the TCT16 or TCT02 data for specified satellite and version 
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

#for index in "${!fileArray[@]}"; do
#  echo $index ${fileArray[$index]}
#done

numFilesPerCpu=$((numberToCompress / (numCpus) + 1))
echo Compressing ${numberToCompress} files using ${numCpus} threads with ${numFilesPerCpu} per thread.


cpuNum=0
counter=0
start=0
stop=0

while [ $cpuNum -lt $numCpus ]; do
	xterm -geometry 80x1 -e "compress_TCT.sh ${compressionDir} ${satellite} ${numCpus} ${cpuNum}" &
	cpuNum=$((cpuNum + 1))
done



