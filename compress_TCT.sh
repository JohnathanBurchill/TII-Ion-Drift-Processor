#!/bin/bash

# TII Cross-Track Ion Drift Processor: compress_TCT.sh

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

# compresses TII cross-track flow CDFs for specified satellite and dirctory.
# http://stackoverflow.com/questions/28226229/how-to-loop-through-dates-using-bash

compressionDir=$1
satellite=$2
numCpus=$3
cpuNum=$4

if test "$#" -ne "4"; then
  echo "Usage: $0 compressionDirectory satelliteLetter numberOfCpus cpuNumber"
  exit 1
fi


# https://stackoverflow.com/questions/238073/how-to-add-a-progress-bar-to-a-shell-script
PROGRESS_BAR_WIDTH=40  # progress bar length in characters
draw_progress_bar() {
  # Arguments: current value, max value, unit of measurement (optional)
  local __value=$1
  local __max=$2
  local __unit=${3:-""}  # if unit is not supplied, do not display it

  # Calculate percentage
  if (( $__max < 1 )); then __max=1; fi  # anti zero division protection
  local __percentage=$(( 100 - ($__max*100 - $__value*100) / $__max ))

  # Rescale the bar according to the progress bar width
  local __num_bar=$(( $__percentage * $PROGRESS_BAR_WIDTH / 100 ))

  # Draw progress bar
  printf "["
  for b in $(seq 1 $__num_bar); do printf "#"; done
  for s in $(seq 1 $(( $PROGRESS_BAR_WIDTH - $__num_bar ))); do printf " "; done
  printf "] $__percentage%% ($__value / $__max $__unit)\r"
}


filesToCompress=`find ${compressionDir} -name "SW_EXPT_EFI${satellite}_*.cdf"`
fileArray=()
for i in ${filesToCompress}; do
  ba=`basename $i .cdf`
 if [ ! -f ${compressionDir}/${ba}.ZIP ]; then
  fileArray+=("${i}")
 fi
done

numberToCompress="${#fileArray[@]}"
index=${cpuNum}
numFilesPerCpu=$((numberToCompress / numCpus + 1))
counter=0

while [ "${index}" -lt "${numberToCompress}" ]; do
        f="${fileArray[${index}]}"
	base=`basename ${f} .cdf`
#	echo Compressing $base.cdf
	draw_progress_bar ${counter} ${numFilesPerCpu} "Compressing files"
	zip -q ${compressionDir}/${base}.ZIP ${f}
        index=$((index + numCpus))
	counter=$((counter + 1))
done


echo "["
