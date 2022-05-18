#!/bin/bash

# TII Cross-Track Ion Drift Processor: parallel_process_TCT.sh

# Copyright (C) 2022  Johnathan K Burchill

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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



