#! /usr/bin/python

"""
Plot results stemming from the running of fouvol.py. See for example,
https://github.com/variationalform/fouvol

Copyright (c) 2020, Simon Shaw
(https://github.com/variationalform, https://www.brunel.ac.uk/people/simon-shaw).
The moral right of the author has been asserted.

These codes are free software; you can redistribute them and/or
modify them under the terms of the GNU General Public License Version 3 - the terms
of which should accompany this source code.
"""

# sudo apt-get install python-tk
# sudo apt-get install python-matplotlib
# sudo apt-get install python-scipy

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import os, math
import itertools
import commands

# we need to find the starting value of Nlim as Lmin from the first data line
Lmin = int(commands.getoutput("head -3 errortable.raw | tail -1 | tr '&' ';' | tr -s ' ' | cut -d';' -f2 | cut -d' ' -f2"))

# this controls whether L doubles or quadruples etc in the legends
DL=2  # a change here must be mirrored in bigrun.sh, compare.py

# Ref: http://www.folkstalk.com/2013/03/sed-remove-lines-file-unix-examples.html
os.system("cat errortable.raw | sed '1,3d' | sed '$d' | tr '&' ' ' | sed 's/\\\\//g' > errortable.dat")
os.system("cat timestable.raw | sed '1,3d' | sed '$d' | tr '&' ' ' | sed 's/\\\\//g' > timestable.dat")

# Ref: https://python4astronomers.github.io/files/asciifiles.html
# Open files, assume they have the same structure
fe = open('errortable.dat', 'r')
ft = open('timestable.dat', 'r')
Nc = len(fe.readline().split())
# count the remaining lines -  there's probably a better way to this but I am in a hurry!
Nl=1;
for line in fe:
  Nl = Nl+1

fe.seek(0,0) # return to beginning of file
errvals = np.zeros((Nl,Nc))
timvals = np.zeros((Nl,Nc))

# Loop over lines and extract variables of interest
linecount = 0
for line in fe:
  line = line.strip()
  columns = line.split()
  for col in range(0,Nc):
    errvals[linecount,col] = float(columns[col])
  linecount = linecount + 1
linecount = 0
for line in ft:
  line = line.strip()
  columns = line.split()
  for col in range(0,Nc):
    timvals[linecount,col] = float(columns[col])
  linecount = linecount + 1
fe.close()
ft.close()
# this causes log(0) errors, and also variatons across runs
#xlo=0; xhi=int(1.2*max(errvals[:,0]));
xlo=1; xhi=1000000; # non-zeros are saferfor log plots

# Ref: https://stackoverflow.com/questions/13091649/unique-plot-marker-for-each-plot-in-matplotlib
marker = itertools.cycle(( '+','s','o','*','v','^','<','>','x','D','d','.'))
# set up the error plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 18})
plt.xlabel(r'$N_t$',fontsize=25)
plt.ylabel(r'error at time $T$',fontsize=25)
plt.tick_params(labelsize=25)

for col in range(1,Nc):
  print 'plotting these:', errvals[:,0], errvals[:,col]
  plt.loglog(errvals[:,0], errvals[:,col], \
    marker = marker.next(), linestyle='-', color='black', \
    linewidth=2, markersize=8, \
    label="$L = "+str(Lmin*2**(DL*(col-1)))+"$" )
#    label="$N_\mathrm{lim} = "+str(Lmin*2**(col-1))+"$" )

plt.legend(loc="lower left")
#plt.legend(loc="upper right")
plt.xlim(xlo,xhi)
plt.ylim(0.000000000001,10.0)
plt.title(r'errors at time $T$',fontsize=25)
plt.tight_layout()
plt.savefig('errors.png', format='png', dpi=750)
#plt.savefig('errors.jpg', format='jpg')  # jpg creation creates errors for me on some machines.
plt.savefig('errors.eps', format='eps', dpi=1000)
plt.grid(True)
plt.clf()

# set up the time plot
xlo=10; xhi=1000000;
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 18})
plt.xlabel(r'$N_t$',fontsize=25)
plt.ylabel(r'time (s)',fontsize=25)
plt.tick_params(labelsize=25)

for col in range(1,Nc):
  print 'plotting these:', timvals[:,0], timvals[:,col]
  plt.loglog(timvals[:,0], timvals[:,col], \
    marker = marker.next(), linestyle='-', color='black', \
    linewidth=2, markersize=8, \
    label="$L = "+str(Lmin*2**(DL*(col-1)))+"$" )
#    label="$N_\mathrm{lim} = "+str(Lmin*2**(col-1))+"$" )

plt.legend(loc="upper left")
plt.xlim(xlo,xhi)
plt.ylim(0.1,100000.0)
plt.title(r'wall clock times (s)',fontsize=25)
plt.tight_layout()
plt.savefig('timings.png', format='png', dpi=750)
#plt.savefig('timings.jpg', format='jpg')
plt.savefig('timings.eps', format='eps', dpi=1000)
plt.grid(True)

exit(0)
