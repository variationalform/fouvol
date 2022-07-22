#! /usr/bin/python

"""
Used to compare results stemming from the running of fouvol.py. See for example,
https://github.com/variationalform/fouvol

Copyright (c) 2020, Simon Shaw
(https://github.com/variationalform, https://www.brunel.ac.uk/people/simon-shaw).
The moral right of the author has been asserted.

These codes are free software; you can redistribute them and/or
modify them under the terms of the GNU General Public License Version 3 - the terms
of which should accompany this source code.
"""

# https://matplotlib.org/1.4.1/api/markers_api.html
# sudo apt-get install python-tk
# sudo apt-get install python-matplotlib
# sudo apt-get install python-scipy

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import os, sys, getopt, time, math
import itertools
import commands

# help function - TO DO: need to decide and specify defaults
def usage():
#  set_log_active(True)
  print("-h   or --help")
  print("-r s or --reference s   directory containing reference data")
  print("-c s or --comparison s  directory containing reference data")
  print("-X s or --gfx s         to block (gfx>0) or hold (gfx<0) screen gfx")
  print(" ")
  os.system('date +%Y_%m_%d_%H-%M-%S')
  print (time.strftime("%d/%m/%Y at %H:%M:%S"))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- --   C O D E    B E G I N S   -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# parse the command line
try:
  opts, args = getopt.getopt(sys.argv[1:], "hr:c:X:",
                   [
                    "help",           # obvious
                    "reference=",     # reference data directory for comparison
                    "comparison=",    # data directory to compare reference to
                    "gfx=",           # to block/hold screen graphics
                    ])

except getopt.GetoptError as err:
  # print help information and exit:
  print(err) # will print something like "option -a not recognized"
  usage()
  sys.exit(2)

# defaults definitions
beingloud  = 20            # level of verbosity
reference  = './'
comparison = './'
gfx        = 0             # screen gfx off by default; otherwise block/hold 
# this controls whether L doubles or quadruples etc in the legends
DL=2  # a change here must be mirrored in bigrun.sh, plotter.py


for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit()
  elif o in ("-r", "--reference"):
    reference = str(a)
    if beingloud > 19: print('using reference  data in = %s;' % reference)
  elif o in ("-c", "--comparison"):
    comparison = str(a)
    if beingloud > 19: print('using comparison data in = %s;' % comparison)
  elif o in ("-X", "--gfx"):
    gfx = float(a)
    if beingloud > 19: print('screen gfx switched on, flag = %f;' % gfx)
  else:
    assert False, "unhandled option"

# rely on the pre-creation of error table.dat and timestable.dat from plotter.py
#os.system("cat errortable.raw | sed '1,3d' | sed '$d' | tr '&' ' ' | sed 's/\\\\//g' > errortable.dat")
#os.system("cat timestable.raw | sed '1,3d' | sed '$d' | tr '&' ' ' | sed 's/\\\\//g' > timestable.dat")

# - - - - - - - - - - - - - - - - 
# read the reference data
# - - - - - - - - - - - - - - - - 
cwd = os.getcwd()
print 'Starting in cwd = ', cwd
os.chdir(reference); # os.system('ls')
# Open files, assume they have the same structure
rfe = open('errortable.dat', 'r')
rft = open('timestable.dat', 'r')
Nc = len(rfe.readline().split())
# count the remaining lines
Nl=1;
for line in rfe:
  Nl = Nl+1

rfe.seek(0,0) # return to beginning of file
rerrvals = np.zeros((Nl,Nc))
rtimvals = np.zeros((Nl,Nc))

# Loop over lines and extract variables of interest in the reference data
linecount = 0
for line in rfe:
  line = line.strip()
  columns = line.split()
  for col in range(0,Nc):
    rerrvals[linecount,col] = float(columns[col])
  linecount = linecount + 1
linecount = 0
for line in rft:
  line = line.strip()
  columns = line.split()
  for col in range(0,Nc):
    rtimvals[linecount,col] = float(columns[col])
  linecount = linecount + 1
rfe.close()
rft.close()
xlo=0; xhi=int(1.2*max(rerrvals[:,0]));

# - - - - - - - - - - - - - - - - 
# read the comparison data
# - - - - - - - - - - - - - - - - 
print 'captured data in', reference
os.chdir(cwd); os.chdir(comparison); # os.system('ls')
cwd = os.getcwd()
print 'Now in cwd = ', cwd
print '...'
# Open files, assume they have the same structure
cfe = open('errortable.dat', 'r')
cft = open('timestable.dat', 'r')

# we need to find the starting value of Nlim as Lmin from the first data line
Lmin = int(commands.getoutput("head -3 errortable.raw | tail -1 | tr '&' ';' | tr -s ' ' | cut -d';' -f2 | cut -d' ' -f2"))

Nc = len(cfe.readline().split())
# count the remaining lines
Nl=1;
for line in cfe:
  Nl = Nl+1

cfe.seek(0,0) # return to beginning of file
cerrvals = np.zeros((Nl,Nc))
ctimvals = np.zeros((Nl,Nc))

# Loop over lines and extract variables of interest in the reference data
linecount = 0
for line in cfe:
  line = line.strip()
  columns = line.split()
  for col in range(0,Nc):
    cerrvals[linecount,col] = float(columns[col])
  linecount = linecount + 1
linecount = 0
for line in cft:
  line = line.strip()
  columns = line.split()
  for col in range(0,Nc):
    ctimvals[linecount,col] = float(columns[col])
  linecount = linecount + 1
cfe.close()
cft.close()
xlo=0; xhi=int(1.2*max(rerrvals[:,0]));

print 'captured data in', comparison
print '...'

# plot the time-error curves
marker = itertools.cycle(( '+','s','o','*','v','^','<','>','x','D','d','.'))
# set up the error plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 16})
plt.xlabel(r'error',fontsize=25)
plt.ylabel(r'execution time (s)',fontsize=25)
plt.tick_params(labelsize=25)

#print 'plotting reference:', rerrvals[:,1], rtimvals[:,1]
plt.loglog(rerrvals[:,1], rtimvals[:,1], \
    marker = marker.next(), linestyle='--', color='black', \
    linewidth=2, markersize=8,
    label="reference" )
# Nc has the number of reference columns
for col in range(1,Nc):
  #print 'plotting comparison:', cerrvals[:,col], ctimvals[:,col]
  plt.loglog(cerrvals[:,col], ctimvals[:,col], \
    marker = marker.next(), linestyle='-', color='black', \
    linewidth=2, markersize=8,
    label="$L = "+str(Lmin*2**(DL*(col-1)))+"$" )

plt.xlim(0.0000001,1000)
plt.ylim(0.1,100000)
plt.legend(loc="upper right")
plt.title(r'error/time comparison',fontsize=25)
plt.tight_layout()
plt.savefig('compare.png', format='png', dpi=750)
#plt.savefig('compare.jpg', format='jpg')
plt.savefig('compare.eps', format='eps', dpi=1000)
plt.grid(True)
if gfx > 0:
  plt.show(block=False); time.sleep(gfx); plt.close()
elif gfx < 0:
  plt.show()
plt.clf()

exit(0)
