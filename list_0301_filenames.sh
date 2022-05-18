#!/bin/bash

# TII Cross-Track Ion Drift Processor: list_0301_filenames.sh

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

(cd /storage/TCT/0301/TCT16 && ls SW_EXPT*.cdf > ../SwarmABC_filenames16_0301.txt)
(cd /storage/TCT/0301/TCT02 && ls SW_EXPT*.cdf > ../SwarmABC_filenames02_0301.txt)

