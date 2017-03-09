#!/usr/bin/env python

"""
Script that takes a command and expands any widlcards before
running it. 
"""
# This file is part of 'python-fmask' - a cloud masking module
# Copyright (C) 2015  Neil Flood
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
from __future__ import print_function, division

import os
import sys
import glob
from rios.parallel.jobmanager import find_executable

if len(sys.argv) < 2:
    raise SystemExit("missing arguments to expand")

outList = [] # will have expanded wildcards
exeName = sys.argv[1]
if exeName.endswith('.py'):
    # look up the path of the script and execute with python
    outList.append(sys.executable)
    scriptPath = find_executable(exeName)
    if scriptPath is None:
        raise SystemExit("Can't find %s" % exeName)
        
    outList.append(scriptPath)
else:
    # normal exe - just append it
    outList.append(exeName)
    
# go through the remaining args and expand them if we can
for arg in sys.argv[2:]:
    expanded = glob.glob(arg)
    if len(expanded) == 0:
        # not a file, just add the original string
        outList.append(arg)
    else:
        outList.extend(expanded)
        
cmd = ' '.join(outList)
#print(cmd)
retval = os.system(cmd)
sys.exit(retval)
