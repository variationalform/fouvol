#!/usr/bin/python


"""
fouvol.py. See for example, https://github.com/variationalform/fouvol

Copyright (c) 2020, Simon Shaw
(https://github.com/variationalform, https://www.brunel.ac.uk/people/simon-shaw).
The moral right of the author has been asserted.

These codes are free software; you can redistribute them and/or
modify them under the terms of the GNU General Public License Version 3 - the terms
of which should accompany this source code.
"""

import os, sys, getopt, time, math
import numpy as np
import matplotlib.pyplot as plt
from fouker import *

# exact solution, u(t) - assumes that tbar = 0
def ut(t,alpha):
  return t**(-alpha)

# RHS function, f(t) - assumes that tbar = 0
def ft(t,varphi0,alpha):
  return t**(-alpha) + varphi0*math.pi*(-alpha)*t/math.sin(-math.pi*alpha)

# help function - TO DO: need to decide and specify defaults
def usage():
#  set_log_active(True)
  print("-h   or --help")
  print("-v n or --garrulous n   to specify verbosity level")
  print("-m n or --smoothness n  to specify C^m for Hermite interpolants")
  print("-s n or --solver n      to specify product solver (1=rectangle; 2=trapezoidal; ... ")
  print("-L n or --Nlim n        to specify \\sum_{n=-N)^N c_n exp(inx) for Fourier Series")
  print("-N n or --Nt n          to specify number of time steps")
  print("-T t or --Time t        to specify final time")
  print("-1 t or --T1 t          to specify N1=ceiling(T1/dt) for the quadrature interval")
  print("-E e or --preverr e     to specify a previous error to compare with")
  print("-a a or --alpha a       to specify \\alpha in \\varphi(t) = \\varphi_0 t^\\alpha")
  print("-X s or --gfx s         to block (gfx>0) or hold (gfx<0) screen gfx")
  print("-P   or --makeplots     if given then plots of phi, Fourier proxy and U(t)")
  print("-F   or --FSerrors      Compute a table of Fourier Series errors and quit (needs Nt large for accurate T1)")
  print("-I   or --invinterp     if given then matrix solve and not divided differences for Hermite interpolants.")
  print(" ")
  os.system('date +%Y_%m_%d_%H-%M-%S')
  print (time.strftime("%d/%m/%Y at %H:%M:%S"))


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- --   C O D E    B E G I N S   -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# parse the command line
try:
  opts, args = getopt.getopt(sys.argv[1:], "hv:m:s:L:N:T:1:E:a:X:PFI",
                   [
                    "help",           # obvious
                    "garrulous=",     # level of verbosity
                    "smoothness=",    # interpolant continuity
                    "solver=",        # specify solver
                    "Nlim=",          # two-sided sum limits for Fourier Series
                    "Nt=",            # number of time steps
                    "Time=",          # final time, T
                    "T1=",            # T1, the quadrature interval - used to obtain N1
                    "preverr=",       # a previous error to compare with
                    "alpha=",         # in \varphi(t) = \varphi_0 t^\alpha
                    "gfx=",           # to block/hold screen graphics
                    "makeplots",      # if given then create plots of phi, Fourier proxy and U(t)
                    "FSerrors",       # Compute a table of Fourier Series errors and quit
                    "invinterp",      # if true use matrix solve and not divided differences for Hermite interpolants.
                    ])

except getopt.GetoptError as err:
  # print help information and exit:
  print(err) # will print something like "option -a not recognized"
  usage()
  sys.exit(2)

# defaults definitions
beingloud     =  0              # level of verbosity
smoothness    =  1              # interpolant continuity
solver        =  1              # default to product rectangle rule
Nlim          =  4              # two-sided sum limits for Fourier Series
Nt            =  4              # number of time steps over (0,T)
N1            =  Nt             # number of time steps over (0,T1)
varphi0       =  0.2
alpha         = -0.4
T1            =  0.5
T             =  math.pi-T1     # final time - WE NEED FOURIER SERIES OVER (-L,L) NOT (-pi,pi)
Tx            =  T+T1
tbar          =  0.0;
Nvals         =  10000          # for smooth plots
gfx           =  0              # screen gfx off by default; otherwise block/hold 
makeplots     =  0              # if true create plots of phi, Fourier proxy and U(t)
getFSerrors   =  0              # get a table of Fourier seies errors and quit
FSerror       =  0.0;
#prev_error    = 0              # previous error (estimate rate of convergence)
DDuse         = 1               # if true use divided differences for Hermite interpolants


for o, a in opts:
  if o in ("-v","--garrulous"):
    beingloud = int(a)
    if beingloud > 19:
      print 'Command Line: using: '
      print('loud level %d;' % beingloud),
  elif o in ("-h", "--help"):
    usage()
    sys.exit()
  elif o in ("-m", "--smoothness"):
    smoothness = int(a)
    if beingloud > 19: print('smoothness = %d;' % smoothness),
  elif o in ("-s", "--solver"):
    solver = int(a)
    if beingloud > 19: print('solver = %d;' % solver),
  elif o in ("-L", "--Nlim"):
    Nlim = int(a)
    if beingloud > 19: print('Fourier sum limits +/- = %d;' % Nlim),
  elif o in ("-N", "--Nt"):
    Nt = int(a)
    if beingloud > 19: print('Number of time steps = %d;' % Nt),
  elif o in ("-T", "--Time"):
    T = float(a)
    if beingloud > 19: print('final time = %f;' % T),
  elif o in ("-1", "--T1"):
    T1 = float(a)
    if beingloud > 19: print('T1 = %f;' % T1),
  elif o in ("-a", "--alpha"):
    alpha = float(a)
    if beingloud > 19: print('alpha = %f;' % alpha),
  elif o in ("-X", "--gfx"):
    gfx = float(a)
    if beingloud > 19: print('screen gfx switched on, flag = %f;' % gfx),
  elif o in ("-P", "--makeplots"):
    makeplots = 1
    if beingloud > 19: print('creating plots of phi, Fourier proxy and U(t);'),
  elif o in ("-F", "--FSerrors"):
    getFSerrors = 1
    if beingloud > 19: print('Tabulating Fourier Series error;'),
  elif o in ("-I", "--invinterp"):
    DDuse = 0
    if beingloud > 19: print('not using divided differences;'),
  else:
    assert False, "unhandled option"

# sanity checks
dt = T/Nt
N1=int(math.ceil(T1/dt))
# crude hack to get into bigrun.sh tables
#os.system('echo '+N1+' > N1_value.txt')
T1 = N1*dt
Tx = T + T1
if beingloud > 19:
  print "\nSanity adjustments give ", 
  print('T = %f; Nt = %d; dt = %f; T1 = %f; N1 = %d; Tx = %f' % (T,Nt,dt,T1,N1,Tx) )

if beingloud > 19:
  # print starting time
  os.system('date +%Y_%m_%d_%H-%M-%S')
  print (time.strftime("%d/%m/%Y at %H:%M:%S"))
  print(" ")
  print('Command line parsing complete... Defaults not yet set and specified')  # TO DO:

# don't spend time on this for the basic rules
if solver > 2:
  old_way=0
  # obtain the Fourier series coefficients and polynomial coefficients
  if old_way == 1:
    cnR, pL, pR = get_FourierCoefficients(tbar,T1,T,Tx,varphi0,alpha,smoothness,Nlim)
  else:
    if DDuse == 0:
      cnR, pL, pR = get_FourierCoefficientsNew(tbar,T1,T,Tx,varphi0,alpha,smoothness,Nlim,DDuse,beingloud)
    else:
      cnR, zL, DDmatL, zR, DDmatR = get_FourierCoefficientsNew(tbar,T1,T,Tx,varphi0,alpha,smoothness,Nlim,DDuse,beingloud)

  if getFSerrors:
    get_FSConvergenceResults(varphi0, alpha, tbar, T, T1, Tx, Nvals, cnR, Nlim, smoothness, gfx, old_way,beingloud)
    exit(0)

  FSerror = get_partialerror(varphi0,alpha,Nvals,T1,T,tbar,cnR,Nlim,Tx)
  if beingloud > 19:
    print 'FSerror = ', FSerror

  if makeplots > 0:
    fig = plt.figure(figsize=(8, 7))
    # these should come before plotting otherwise the fonts get messed up
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rcParams.update({'font.size': 18})
    plot_varphi(varphi0,alpha,Nvals,T,tbar)
    plot_partvarphi(varphi0,alpha,Nvals,T,T1,tbar)
    plot_FS(Tx,Nvals,cnR,Nlim)
    if old_way == 1:
      plot_leftHermite(Nvals,T1,pL)
      plot_rightHermite(Nvals,T,Tx,pR)
    else:
      if DDuse == 0:
        plot_basicHermite(Nvals,0,T1,pL,2*smoothness+1,'-b',r'$p_L(t)$')
        plot_basicHermite(Nvals,T,Tx,pR,2*smoothness+1,'-r',r'$p_R(t)$')
      else:
        plot_DDHermite(Nvals,0,T1,zL,DDmatL,smoothness,'-b',r'$p_L(t)$')
#        plot_basicHermite(Nvals,0,T1,pL,2*smoothness+1,'-b',r'$p_L(t)$')
        plot_DDHermite(Nvals,T,Tx,zR,DDmatR,smoothness,'-r',r'$p_R(t)$')
#        plot_basicHermite(Nvals,T,Tx,pR,2*smoothness+1,'-r',r'$p_R(t)$')
    plt.xlabel(r'time, $t$',fontsize=25)
    plt.ylabel(r'$\varphi(t)$',fontsize=25)
    plt.tick_params(labelsize=25)
    plt.legend(loc="upper right")
    plt.xlim(0, Tx);
    plt.ylim(0, 1.3*get_LeftHermiteInterceptValue(varphi0,T1,tbar,alpha))
    plt.title(r'$\varphi(t)$ with FSerror = '+str(FSerror),fontsize=25)
    plt.tight_layout()
    plt.savefig('varphirep.png', format='png', dpi=750)
#    print 'WARNING: fouvol.py - jpg needs uncommenting for varphirep '
#    plt.savefig('varphirep.jpg', format='jpg')
    plt.savefig('varphirep.eps', format='eps', dpi=1000)
    if gfx > 0:
      plt.show(block=False); time.sleep(gfx); plt.close()
    elif gfx < 0:
      plt.show()
    plt.clf()
    
    Nper=4
    fig = plt.figure(figsize=(10, 5))
    plot_varphi(varphi0,alpha,Nvals,T,tbar)
    plot_partvarphi(varphi0,alpha,Nvals,T,T1,tbar)
    if old_way == 1:
      plot_leftHermite(Nvals,T1,pL)
      plot_rightHermite(Nvals,T,Tx,pR)
    else:
      if DDuse == 0:
        plot_basicHermite(Nvals,0,T1,pL,2*smoothness+1,'-b',r'$p_L(t)$')
        plot_basicHermite(Nvals,T,Tx,pR,2*smoothness+1,'-r',r'$p_R(t)$')
      else:
        plot_DDHermite(Nvals,0,T1,zL,DDmatL,smoothness,'-b',r'$p_L(t)$')
#       plot_basicHermite(Nvals,0,T1,pL,2*smoothness+1,'-b',r'$p_L(t)$')
        plot_DDHermite(Nvals,T,Tx,zR,DDmatR,smoothness,'-r',r'$p_R(t)$')
#        plot_basicHermite(Nvals,T,Tx,pR,2*smoothness+1,'-r',r'$p_R(t)$')
    plot_long_FS(Nper,Tx,2*Nper*Nvals,cnR,Nlim,varphi0,alpha,T,tbar)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rcParams.update({'font.size': 18})
    plt.xlabel(r'time, $t$',fontsize=25)
    plt.ylabel(r'$\varphi(t)$, $\psi(t)$ and $F\psi(t)$',fontsize=25)
    plt.tick_params(labelsize=25)
    plt.legend(loc="upper right")
    plt.xlim(-Nper*Tx, Nper*Tx);
    plt.ylim(0,0.7)
    #plt.ylim(0, 1.3*get_LeftHermiteInterceptValue(varphi0,T1,tbar,alpha))
    plt.title(r'$\varphi(t)$, the splines, and the periodic extension',fontsize=25)
    plt.tight_layout()
    plt.savefig('varphirepextended.png', format='png', dpi=750)
#    print 'WARNING: fouvol.py - jpg needs uncommenting for varphirepextended '
#    plt.savefig('varphirepextended.jpg', format='jpg')
    plt.savefig('varphirepextended.eps', format='eps', dpi=1000)
    if gfx > 0:
      plt.show(block=False); time.sleep(gfx); plt.close()
    elif gfx < 0:
      plt.show()
    plt.clf()

# begin solution process, store solution here - avoid solving for u(0)
uvals = np.zeros(Nt+1)  # don't need the ,1) here - it makes abs() return [xxx] not xxx WHERE ELSE?
cheat = 0
# set up the complex Fourier history - whether or not it is needed
phistarR = np.zeros(2*Nlim+1) 
phistarI = np.zeros(2*Nlim+1) 

# product rectangle rule
if solver == 1:
  # set up the denominator for product rectangle rule
  if tbar > 0:
    denom = 1.0 + varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + varphi0/(alpha+1)*( dt**(alpha+1) )
  # begin time stepping
  for i in range(1,Nt+1):
    # get history, list comprehension seems around 5% faster
    hist = sum([varphi0/(alpha+1)*(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*uvals[j] for j in range(1,i)])
    uvals[i] = ( ft(i*dt,varphi0,alpha) - hist ) / denom
    if cheat and i<Nt/4:
      uvals[i] = ut(i*dt,alpha)

# product trapezoidal rule for all except first interval
elif solver == 2:
  # set up the denominator for product rectangle rule for first interval
  if tbar > 0:
    denom = 1.0 + varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + varphi0/(alpha+1)*( dt**(alpha+1) )
  # product rectangle rule on first interval, convenience setting for u[0] (used in history sum)
  uvals[1] = ft(dt,varphi0,alpha) / denom
  uvals[0] = uvals[1]
  # set up the denominator for product trapezoidal rule
  if tbar > 0:
    denom = 1.0 + 0.5*varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + 0.5*varphi0/(alpha+1)*( dt**(alpha+1) )
  # begin time stepping
  for i in range(2,Nt+1):
    hist = sum([0.5*varphi0/(alpha+1)*(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*(uvals[j]+uvals[j-1]) \
                   for j in range(1,i)])
    hist = hist + 0.5*varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) ) * uvals[i-1]
    uvals[i] = ( ft(i*dt,varphi0,alpha) - hist ) / denom
    if cheat and i<Nt/4:
      uvals[i] = ut(i*dt,alpha)

# Fourier, with product rectangle rule
elif solver == 3:
  # set up the denominator for product rectangle rule
  if tbar > 0:
    denom = 1.0 + varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + varphi0/(alpha+1)*( dt**(alpha+1) )
  # begin time stepping
  for i in range(1,Nt+1):
    # update Fourier history variables from previous time step
    if i > N1:
      for j in range(-Nlim,0):
        Rpart = phistarR[j+Nlim]; Ipart = phistarI[j+Nlim];
        costerm = math.cos(math.pi*j*dt/Tx); sinterm = math.sin(math.pi*j*dt/Tx)
        cosT1   = math.cos(math.pi*j*T1/Tx); sinT1   = math.sin(math.pi*j*T1/Tx)
        # complex attenuation
        phistarR[j+Nlim] = costerm*Rpart - sinterm*Ipart
        phistarI[j+Nlim] = sinterm*Rpart + costerm*Ipart
        # product integration improvement - the old way used the rectangle rule
        update = cosT1*sinterm-sinT1*(1.0-costerm)
        coeff = cnR[j+Nlim]*uvals[i-N1]*Tx/j/math.pi
        phistarR[j+Nlim] = phistarR[j+Nlim] + coeff*update
        update = cosT1*(1.0-costerm)+sinT1*sinterm
        phistarI[j+Nlim] = phistarI[j+Nlim] + coeff*update
      # the special case where j=0 - no 'attenuation' is necessary
      phistarR[0+Nlim] = phistarR[0+Nlim] + cnR[0+Nlim]*uvals[i-N1]*dt

    # quadrature over j={i-N1,...,i-2,i-1}
    jlo = int(max(1,i+1-N1))
    hist = varphi0/(alpha+1)*sum([(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*uvals[j] for j in range(jlo,i)])
    # Fourier contribution for far history
    if i > N1:
      hist = hist + phistarR[0+Nlim]
      hist = hist + 2*sum([phistarR[j+Nlim] for j in range(-Nlim,0)])
    # solve...
    uvals[i] = ( ft(i*dt,varphi0,alpha) - hist ) / denom
    if cheat and i<Nt/4:
      uvals[i] = ut(i*dt,alpha)

# Fourier, with product trapezoidal rule
elif solver == 4:
  # set up the denominator for product rectangle rule for first interval
  if tbar > 0:
    denom = 1.0 + varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + varphi0/(alpha+1)*( dt**(alpha+1) )
  # product rectangle rule on first interval, convenience setting for u[0] (used in history sum)
  uvals[1] = ft(dt,varphi0,alpha) / denom
  uvals[0] = uvals[1]
  # set up the denominator for product trapezoidal rule
  if tbar > 0:
    denom = 1.0 + 0.5*varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + 0.5*varphi0/(alpha+1)*( dt**(alpha+1) )
  # begin time stepping
  for i in range(2,Nt+1):
    # update Fourier history variables from previous time step
    if i > N1:
      for j in range(-Nlim,0):
        Rpart = phistarR[j+Nlim]; Ipart = phistarI[j+Nlim];
        costerm = math.cos(math.pi*j*dt/Tx); sinterm = math.sin(math.pi*j*dt/Tx)
        cosT1   = math.cos(math.pi*j*T1/Tx); sinT1   = math.sin(math.pi*j*T1/Tx)
        # complex attenuation
        phistarR[j+Nlim] = costerm*Rpart - sinterm*Ipart
        phistarI[j+Nlim] = sinterm*Rpart + costerm*Ipart
        # product integration
        update = cosT1*sinterm-sinT1*(1.0-costerm)
        coeff = cnR[j+Nlim]*(uvals[i-N1]+uvals[i-1-N1])/2.0*Tx/j/math.pi
        phistarR[j+Nlim] = phistarR[j+Nlim] + coeff*update
        update = cosT1*(1.0-costerm)+sinT1*sinterm
        phistarI[j+Nlim] = phistarI[j+Nlim] + coeff*update
      # the special case where j=0 - no 'attenuation' is necessary
      phistarR[0+Nlim] = phistarR[0+Nlim] + cnR[0+Nlim]*(uvals[i-N1]+uvals[i-1-N1])/2.0*dt

    # quadrature over j={i-N1,...,i-2,i-1}
    jlo = int(max(1,i+1-N1))
    hist = sum([0.5*varphi0/(alpha+1)*(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*(uvals[j]+uvals[j-1]) \
           for j in range(jlo,i)])
    hist = hist + 0.5*varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )*uvals[i-1]
    # Fourier contribution for far history
    if i > N1:
      hist = hist + phistarR[0+Nlim]
      hist = hist + 2*sum([phistarR[j+Nlim] for j in range(-Nlim,0)])
    # solve...
    uvals[i] = ( ft(i*dt,varphi0,alpha) - hist ) / denom
    if cheat and i<Nt/4:
      uvals[i] = ut(i*dt,alpha)

# Sparse quad, rect, trap, Simp and Boole
# FAOCQ
# McLean's panel method
#elif solver == ?:
  # set up the denominator for product rectangle rule for first interval

# calc max error
error = abs(ut(T,alpha)-uvals[Nt])
print 'FSerror = ', FSerror
print '|u(T)-U| = ', error

# plot the solution over [dt,T] - using list comprehension
def plot_ut(uvals,T,Nt,gfx):
  tvals = [float(0.0+i*(T-0.0)/Nt) for i in range(1,Nt+1)]
  ut    = [ uvals[i] for i in range(1,Nt+1) ]
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  plt.plot(tvals, ut, "-k", label=r'$U(t)$')
  plt.rcParams.update({'font.size': 18})
  plt.xlabel(r'time, $t$',fontsize=25)
  plt.ylabel(r'$U(t)$',fontsize=25)
  plt.tick_params(labelsize=25)
  plt.xlim(0, T);
  plt.ylim(0, 1.1*max(abs(uvals)))
  plt.title(r'approximation, $U(t)$',fontsize=25)
  plt.tight_layout()
  plt.savefig('Ut.png', format='png', dpi=750)
#  print 'WARNING: fouvol.py - Ut.jpg needs uncommenting'
#  plt.savefig('Ut.jpg', format='jpg')
  plt.savefig('Ut.eps', format='eps', dpi=1000)
  if gfx > 0:
    plt.show(block=False); time.sleep(gfx); plt.close()
  elif gfx < 0:
    plt.show()
  plt.clf()

# plot the solution
if makeplots > 0:
  plot_ut(uvals,T,Nt,gfx)

if beingloud > 19:
  # print finish time
  print('Finished at...')
  os.system('date +%Y_%m_%d_%H-%M-%S')
  print (time.strftime("%d/%m/%Y at %H:%M:%S"))

