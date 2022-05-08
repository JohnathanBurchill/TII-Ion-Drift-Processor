#!/bin/bash

# processes TCT16 and TCT02 data for specified satellite and version for a range of dates
# http://stackoverflow.com/questions/28226229/how-to-loop-through-dates-using-bash

satellite=$1
startDate=$2
stopDate=$3
calVersion=$4
exportVersion=$5
numCpus=$6

if test "$#" -ne "6"; then
  echo "Usage: $0 satelliteLetter startDate stopDate calVersion exportVersion numCpus"
  exit 1
fi

exportDir="/databases/TCT"

processingDate=$(date +%Y%m%dT%H%M%S)

dateToProcess="$startDate"
cpuNum=0

while [ $cpuNum -lt $numCpus ]; do
	cpuNum=$((cpuNum + 1))
	logFile="${exportDir}/${exportVersion}/TCT16/Swarm${satellite}_TCT16_${exportVersion}_processing_log_${processingDate}_cpu_${cpuNum}.log"
	xterm -geometry 80x1 -e "process_TCT.sh ${satellite} ${dateToProcess} ${stopDate} $numCpus $calVersion $exportVersion $exportDir $logFile" &
	dateToProcess=$(date -I -d "$dateToProcess + 1 day")
done



