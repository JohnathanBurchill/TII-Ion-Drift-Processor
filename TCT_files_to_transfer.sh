#!/bin/bash

# TII Cross-Track Ion Drift Processor: TCT_files_to_transfer.sh

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

