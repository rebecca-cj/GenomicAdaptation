## Calculate mean covariance matrix from X matrices (for Bayenv2) ##

#  Copyright 2016 Rebecca Jordan
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


# To be used in conjunction with meanCovMat.sh
# Input:  Output from meanCovMat.sh = series of files, each containing a single output covariance matrix; 
#	  file names = "file_prefix.X" where X is a series from 1 to the number of matrix files
# Usage: meanCovMat.py <file_prefix> <number_of_matrices_to_be_averages>

import sys
import numpy as np

prefix=sys.argv[1]
print prefix
reps=sys.argv[2]
print reps

covmat={}
for i in range(1,int(reps)+1):
	covmat[i]=np.loadtxt(prefix+'.'+str(i))

tmp=0
for i in range(1,int(reps)+1):
	tmpX=np.add(tmp,covmat[i])
	tmp=tmpX
meanmatrix=tmp/int(reps)

np.savetxt(prefix+'.mean',meanmatrix,delimiter='\t')
