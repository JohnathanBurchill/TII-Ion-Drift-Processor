#!/bin/bash

# TII Cross-Track Ion Drift Processor: cdf_files_to_compress.txt

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



