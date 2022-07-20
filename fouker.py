#!/usr/bin/python

"""

"""

from scipy.integrate import quad
import scipy
import os, sys, getopt, time, math
import numpy as np
import matplotlib.pyplot as plt
import itertools

# plot \varphi(t) over [epsilon,T] - using list comprehension
def plot_varphi(varphi0,alpha,Nvals,T,tbar):
  tvals = [float(0.0+i*(T-0.0)/Nvals) for i in range(1,Nvals+1)]
  phivals = [ varphi0*(time+tbar)**alpha for time in tvals ]
  plt.plot(tvals, phivals, ":g", linewidth=2, markersize=0) # no label, that appears in the plot below

# plot \varphi(t) over [T1,T]
def plot_partvarphi(varphi0,alpha,Nvals,T,T1,tbar):
  tvals = [float(T1+i*(T-T1)/Nvals) for i in range(0,Nvals+1)]
  phivals = [ varphi0*(time+tbar)**alpha for time in tvals ]
  plt.plot(tvals, phivals, "-g", label=r'$\varphi(t)$', linewidth=2, markersize=0)

# plot pC(t) - plot the basic Hermite polynomial over \tau\in[0,1]
def plot_basicHermite(Nvals,a,b,pC,degree,line_str,label_str):
  tvals = [float(a+i*(b-a)/Nvals) for i in range(0,Nvals+1)]
  pCvals = [ get_MappedHornerPolyValue(pC,a,b,time,degree) for time in tvals ]
  plt.plot(tvals, pCvals, line_str, label=label_str, linewidth=2, markersize=0)

# plot the Hermite interpolant using DD's 
def plot_DDHermite(Nvals,a,b,z,DDmat,smoothness,line_str,label_str):
  tvals = [float(a+i*(b-a)/Nvals) for i in range(0,Nvals+1)]
  pCvals = [ get_DDHermiteValue(time,DDmat,smoothness,z) for time in tvals ]
#  pCvals = [ get_MappedHornerPolyValue(pC,a,b,time,degree) for time in tvals ]
  plt.plot(tvals, pCvals, line_str, label=label_str, linewidth=2, markersize=0)

# plot pL(t) - the left side Hermite polynomial over [0,T1]
def plot_leftHermite(Nvals,T1,pL):
  tvals = [float(0.+i*(T1-0.)/Nvals) for i in range(0,Nvals+1)]
  pLvals = [ np.polyval(np.flipud(pL),time) for time in tvals ]
  plt.plot(tvals, pLvals, "-b", label=r'$p_L(t)$', linewidth=2, markersize=0)

# plot pR(t) - the right side Hermite polynomial over [T,Tx]
def plot_rightHermite(Nvals,T,Tx,pR):
  tvals = [float(T+i*(Tx-T)/Nvals) for i in range(0,Nvals+1)]
  pRvals = [ np.polyval(np.flipud(pR),time) for time in tvals ]
  plt.plot(tvals, pRvals, "-r", label=r'$p_R(t)$', linewidth=2, markersize=0)

# plot the Fourier series
def plot_FS(Tx,Nvals,cnR,Nlim):
  tvals = [float(0.+i*(Tx-0.)/Nvals) for i in range(0,Nvals+1)]
  FSvals = [ sum([cnR[k+Nlim]*math.cos(math.pi*k*tvals[i]/Tx) for k in range(-Nlim,Nlim+1,1)]) for i in range(0,Nvals+1,1)]
  plt.plot(tvals, FSvals, "-k", label=r'$F\psi(t)$', linewidth=2, markersize=0)
  return FSvals

# plot the Fourier series and periodic extension over -Nper to Nper periods
def plot_long_FS(Nper,Tx,Nvals,cnR,Nlim,varphi0,alpha,T,tbar):
  tvals = [float(-Nper*Tx+i*(2*Nper*Tx-0.)/Nvals) for i in range(0,Nvals+1)]
  FSvals = [ sum([cnR[k+Nlim]*math.cos(math.pi*k*tvals[i]/Tx) for k in range(-Nlim,Nlim+1,1)]) for i in range(0,Nvals+1,1)]
#  fig = plt.figure(figsize=(8, 6))
#  ax = plt.axes()  #  ax = plt.gca()    # get the axis handle
#  ax.set_ylim([0.0, 0.7])
  plt.plot(tvals, FSvals, "-k", label=r'$F\psi(t)$', linewidth=1, markersize=0)
  #plt.ylim(0,0.7) # plt.tight_layout()
  #plt.ylim(top=0.7)
#  ax.set_aspect(80) # sets the height to width ratio
  return FSvals

# get a vector of Fourier series values at a vector of times
def get_FS(Tx,Nvals,cnR,Nlim):
  tvals = [float(0.+i*(Tx-0.)/Nvals) for i in range(0,Nvals+1)]
  FSvals = [ sum([cnR[k+Nlim]*math.cos(math.pi*k*tvals[i]/Tx) for k in range(-Nlim,Nlim+1,1)]) for i in range(0,Nvals+1,1)]
  return FSvals

# get the Hermite coefficients for \tau\in[0,1] with RHS scaling
def get_basicHermite(n,f,fscale,loudness):
  A = np.zeros((2*n+2,2*n+2)); pC=np.zeros((2*n+2,1)); coeff = np.ones((2*n+2,1))
  for row in range(0,2*n+1,2):
    apow = 1;
    f[row  ] = f[row  ]*fscale**(row/2)
    f[row+1] = f[row+1]*fscale**(row/2)
    for col in range(row/2,2*n+2,1):
      A[row  ,col] = coeff[col]*apow;
      A[row+1,col] = coeff[col];
      coeff[col]   = coeff[col] * (col-row/2);
      apow = 0.0
  pC = np.linalg.solve(A, f)
  #Check that the solution is correct:
  if loudness > 0:
    residnorm = np.linalg.norm(f - np.dot(A, pC), np.inf)
    if residnorm > 0.01:
      print 'WARNING: Hermite linear solve has large residual, allclose may be FALSE'
      print np.allclose(np.dot(A, pC), f)
      print residnorm
  return pC

## get a C^m Hermite interpolant divided differences matrix
#def get_DDHermite(a,b,m,fa,fb):
  ## an empty matrix, and a vector of duplicated x values
  #DDmat = np.zeros((2*m+2,2*m+2))
  #z     = np.zeros(2*m+2)
  ## populate the duplicated values
  #for k in range(0,m+1):
    #z[k] = a; z[m+1+k] = b;
  ## intialize the matrix
  #for k in range(0,m+1):
    #DDmat[k,k] = fa;
    #DDmat[k+m+1,k+m+1] = fb;   
  ## now fill the matrix up, do the diagonals for duplicated points first
  #for n in range(1,m+1):        # update super-diagonal 1, 2, ..., m
    #for k in range(m-n+1, 2*m-1):    # update rows 0, 1, ..., 2m+1-n 
      #DDmat[k,k+n] = ( DDmat[k+1,k+n]-DDmat[k,k+n-1] ) / ( z[k+n]-z[k] )
  ## finish off with the remaining super diagonals
  #for n in range(m+1,2*m+2):        # update super-diagonal m+1, ..., 2m+1
    #for k in range(0, 2*m+2-n):     # update rows 0, 1, ..., 2m+1-n
      #DDmat[k,k+n] = ( DDmat[k+1,k+n]-DDmat[k,k+n-1] ) / ( z[k+n]-z[k] )
##  print z
##  print ' '
##  print DDmat
  #return DDmat, z

# get C^m Hermite interpolant divided differences matrix with m derivs at a and b
def get_DDHermite(a,b,m,fa,fb,da,db):
  # an empty matrix, and a vector of duplicated x values
  DDmat = np.zeros((2*m+2,2*m+2))
  z     = np.zeros(2*m+2)
  # populate the duplicated values
  for k in range(0,m+1):
    z[k] = a; z[m+1+k] = b;
  # intialize the matrix
  for k in range(0,m+1):
    DDmat[k,k] = fa;
    DDmat[k+m+1,k+m+1] = fb;
  # now fill the matrix up, do the diagonals for duplicated points first
  nfac = 1
  for n in range(1,m+1):        # update super-diagonal zero blocks 
    nfac = n*nfac
    for k in range(0,m+1-n):    # update rows from top and bottom 
      DDmat[k,  k+n]       = da[n]/nfac
      DDmat[m+k+1,m+k+n+1] = db[n]/nfac
  # finish off with the remaining super diagonal sub block, bottom up
  for k in range(m,-1,-1):        # update upper right square block
    for n in range(m+1, 2*m+2):     # update rows 0, 1, ..., 2m+1-n
      DDmat[k,n] = ( DDmat[k+1,n]-DDmat[k,n-1] ) / ( b-a ) # denom is always b-a
  #print z
  #print ' '
  #print DDmat
  return DDmat, z

# get the value of the Divided difference version of the Hermite at t
def get_DDHermiteValue(t,DDmat,m,z):
  total = 0.0
  coeff = 1.0
  for k in range(0,2*m+2):
    total = total+DDmat[0,k]*coeff
    coeff = coeff*(t-z[k])
  return total

## get the left side Hermite coefficients
#def get_leftHermite(a,b,n,fL):
  #A = np.zeros((2*n+2,2*n+2)); pL=np.zeros((2*n+2,1)); coeff = np.ones((2*n+2,1))
  #for row in range(0,2*n+1,2):
    #apow = 1; bpow = 1;
    #for col in range(row/2,2*n+2,1):
      #A[row  ,col] = coeff[col] * apow;
      #A[row+1,col] = coeff[col] * bpow;
      #coeff[col]   = coeff[col] * (col-row/2);
      #apow = a*apow; bpow = b*bpow;
  #pL = np.linalg.solve(A, fL)
  ##Check that the solution is correct:
  #print np.allclose(np.dot(A, pL), fL)
  #return pL

## get the right side Hermite coefficients
#def get_rightHermite(a,b,n,fR):
  #A = np.zeros((2*n+2,2*n+2)); pR=np.zeros((2*n+2,1)); coeff = np.ones((2*n+2,1))
  #for row in range(0,2*n+1,2):
    #apow = 1; bpow = 1;
    #for col in range(row/2,2*n+2,1):
      #A[row  ,col] = coeff[col] * apow;
      #A[row+1,col] = coeff[col] * bpow;
      #coeff[col]   = coeff[col] * (col-row/2);
      #apow = a*apow; bpow = b*bpow;
  #pR = np.linalg.solve(A, fR)
  ##Check that the solution is correct:
  #print np.allclose(np.dot(A, pR), fR)
  #return pR

# calculate error over the part that matters
def get_partialerror(varphi0,alpha,Nvals,T1,T,tbar,cnR,Nlim,Tx):
  error = 0
  tvals = np.zeros(Nvals+1)
  FSvals = np.zeros(Nvals+1)
  phivals = np.zeros(Nvals+1)
  for i in range(0,Nvals+1,1):
    t        = float(T1+i*(T-T1)/Nvals)
    tvals[i] = t
    FSvals[i] = 0
    # earlier results had range(-Nlim,Nlim,1)
    for k in range(-Nlim,Nlim+1,1):
      FSvals[i] = FSvals[i] + cnR[k+Nlim]*math.cos(math.pi*k*t/Tx)
    phivals[i] = varphi0*(t+tbar)**alpha
    error = max(abs(FSvals[i]-phivals[i]), error)
  return error

def get_HermiteInterceptValues(lbeta,rbeta,varphi0,T1,T,Tx,tbar,alpha):
  lval = varphi0*(T1+tbar)**alpha - T1*lbeta*alpha*varphi0*(T1+tbar)**(alpha-1)
  rval = varphi0*(T+tbar)**alpha+(1-rbeta)*(Tx-T)*alpha*varphi0*(T+tbar)**(alpha-1)
  return lval, rval

def get_LeftHermiteInterceptValue(varphi0,T1,tbar,alpha):
  return 1.4*varphi0*(T1+tbar)**alpha

def get_RightHermiteInterceptValue(varphi0,T,tbar,alpha):
  return 0.95*varphi0*(T+tbar)**alpha

## obsolete - and needs updating to non 2\pi period, and intercept value implementation needs checking
#def get_FourierCoefficients(tbar,T1,T,Tx,varphi0,alpha,smoothness,Nlim):
  #lval, rval = get_HermiteInterceptValues(0.5,0.5,varphi0,T1,T,Tx,tbar,alpha)
  ##calculate the Hermite at the left end
  #a=0.; b=T1; #n=smoothness;
  #fL=np.zeros((2*smoothness+2,1)); 
  #fL[0] = lval
  #fL[1] = varphi0*(T1+tbar)**alpha
  #for k in range(2,2*smoothness+2,2):
    #fL[k]   = 0.0
    #fL[k+1] = fL[k-1] * (alpha-(k-2.0)/2.0) * (T1+tbar)**(-1)
  #pL = get_leftHermite(a, b, smoothness, fL)

  ##calculate the Hermite at the right end
  #a=T; b=Tx; # n=2;
  #fR=np.zeros((2*smoothness+2,1)); 
  #fR[0] = varphi0*(T+tbar)**alpha
  #fR[1] = rval
  #for k in range(2,2*smoothness+2,2):
    #fR[k]   = fR[k-2] * (alpha-(k-2.0)/2.0) * (T+tbar)**(-1)
    #fR[k+1] = 0.0
  #pR = get_rightHermite(a, b, smoothness, fR)

  ## calculate Fourier series
  #cnR= np.zeros((2*Nlim+1,1)) 
  ## basic assumption: Im part = 0, Re part symmetric about n=0. Save time.
  #for n in range(0,Nlim+1):
  ##  print realpart % realpart[1] has error
    #realpart = quad(lambda t: np.polyval(np.flipud(pL),t) * math.cos(n*t), 0, T1, limit=1000)
    #cnR[n+Nlim] = 1.0/math.pi*realpart[0]
    #realpart = quad(lambda t: varphi0*(t+tbar)**alpha * math.cos(n*t), T1, T, limit=1000)
    #cnR[n+Nlim] = cnR[n+Nlim] + 1.0/math.pi*realpart[0]
    #realpart = quad(lambda t: np.polyval(np.flipud(pR),t) * math.cos(n*t), T, Tx, limit=1000)
    #cnR[n+Nlim] = cnR[n+Nlim] + 1.0/math.pi*realpart[0]

  #for n in range(-Nlim,0):
    #cnR[n+Nlim] = cnR[-n+Nlim]

  #return cnR, pL, pR

def get_MappedHornerPolyValue(pC,a,b,t,degree):
  tau = (t-a)/(b-a)
  value = pC[degree]
  for k in range(degree-1,-1,-1):
    value = pC[k] + tau*value
  return value

def get_FourierCoefficientsNew(tbar,T1,T,Tx,varphi0,alpha,smoothness,Nlim,DDuse,loudness):
  lval, rval = get_HermiteInterceptValues(0.9,0.5,varphi0,T1,T,Tx,tbar,alpha)
  #calculate the Hermite at the left end
  if DDuse == 1:
    # get the function values
    fa = lval; fb = varphi0*(T1+tbar)**alpha;
    # get the derivative values (all zero on left), d[1] = first deriv etc
    da = np.zeros(1+smoothness); db = np.zeros(1+smoothness)
    deriv = alpha*fb/(T1+tbar) # init to first deriv at T1
    for n in range(1,1+smoothness):
      db[n] = deriv
      deriv = (alpha-n)*deriv/(T1+tbar)
    DDmatL, zL = get_DDHermite(0,T1,smoothness,fa,fb,da,db)
  else:    
    a=0.; b=1.
    fL=np.zeros((2*smoothness+2,1)); 
    fL[0] = lval
    fL[1] = varphi0*(T1+tbar)**alpha
    for k in range(2,2*smoothness+2,2):
      fL[k]   = 0.0
      fL[k+1] = fL[k-1] * (alpha-(k-2.0)/2.0) * (T1+tbar)**(-1)
    pL = get_basicHermite(smoothness,fL,T1-0.,loudness)
  
  #mtmp = 6
  #da = np.arange(mtmp); db = 1+np.arange(mtmp)
  #DDmatL, zL = get_DDHermite(0,0.5,mtmp,2,4,da,db)
  #print get_DDHermiteValue(0.0,DDmatL,2,zL)
  #print get_DDHermiteValue(0.5,DDmatL,2,zL)
  #exit()

  #calculate the Hermite at the right end
  if DDuse == 1:
    # get the function values
    fa = varphi0*(T+tbar)**alpha; fb = rval
    # get the derivative values (all zero on right), d[1] = first deriv etc
    da = np.zeros(1+smoothness); db = np.zeros(1+smoothness)
    deriv = alpha*fb/(T+tbar) # init to first deriv at T1
    for n in range(1,1+smoothness):
      da[n] = deriv
      deriv = (alpha-n)*deriv/(T+tbar)
    DDmatR, zR = get_DDHermite(T,T+Tx,smoothness,fa,fb,da,db)
  else:    
    a=0.; b=1.
    fR=np.zeros((2*smoothness+2,1)); 
    fR[0] =     varphi0*(T+tbar)**alpha
    fR[1] = rval
    for k in range(2,2*smoothness+2,2):
      fR[k]   = fR[k-2] * (alpha-(k-2.0)/2.0) * (T+tbar)**(-1)
      fR[k+1] = 0.0
    pR = get_basicHermite(smoothness,fR,Tx-T,loudness)
  
  # calculate Fourier series
  cnR= np.zeros((2*Nlim+1,1)) 
  # basic assumption: Im part = 0, Re part symmetric about n=0. Save time.
  for n in range(0,Nlim+1):
  #  print realpart % realpart[1] has error
    if DDuse == 1:
      realpart = quad(lambda t: get_DDHermiteValue(t,DDmatL,smoothness,zL) * math.cos(math.pi*n*t/Tx), 0, T1, limit=80000)
    else:
      realpart = quad(lambda t: get_MappedHornerPolyValue(pL,0,T1,t,2*smoothness+1) * math.cos(math.pi*n*t/Tx), 0, T1, limit=80000)
    cnR[n+Nlim] = 1.0/Tx*realpart[0]
    realpart = quad(lambda t: varphi0*(t+tbar)**alpha * math.cos(math.pi*n*t/Tx), T1, T, limit=80000)
    cnR[n+Nlim] = cnR[n+Nlim] + 1.0/Tx*realpart[0]
    if DDuse == 1:
      realpart = quad(lambda t: get_DDHermiteValue(t,DDmatR,smoothness,zR) * math.cos(math.pi*n*t/Tx), T, Tx, limit=80000)
    else:
      realpart = quad(lambda t: get_MappedHornerPolyValue(pR,T,Tx,t,2*smoothness+1) * math.cos(math.pi*n*t/Tx), T, Tx, limit=80000)
    cnR[n+Nlim] = cnR[n+Nlim] + 1.0/Tx*realpart[0]

  for n in range(-Nlim,0):
    cnR[n+Nlim] = cnR[-n+Nlim]

  if DDuse == 1:
    return cnR, zL, DDmatL, zR, DDmatR
  else:
    return cnR, pL, pR

# e.g. time ./fouvol.py  -v 20 -a -0.51 --T1 0.2 -T 1 -F --Nt 1000000 -s 3
def get_FSConvergenceResults(varphi0, alpha, tbar, T, T1, Tx, Nvals, cnR, Nlim, smoothness, gfx, old_way, loudness):
  max_m=10; max_L=12
#  max_m=4; max_L=4
#  max_m=6; max_L=6
  errors = np.zeros((max_L,max_m)) 
  ratios = np.zeros((max_L,max_m)) 
  
  print "Wait, computing errors... (Did you set Nt large enough to allow N1 to be properly defined?)"
  for i in range(0,max_L):
    Nlim = 2**i
    for m in range(1,max_m):
      print 'Currently on Nlim = ', Nlim, ' m = ', m
      smoothness = m
      if old_way == 1:
        cnR, pL, pR = get_FourierCoefficients(tbar,T1,T,Tx,varphi0,alpha,smoothness,Nlim)
      else:
        cnR, pL, pR = get_FourierCoefficientsNew(tbar,T1,T,Tx,varphi0,alpha,smoothness,Nlim,DDuse,loudness)
      FSerror = get_partialerror(varphi0,alpha,Nvals,T1,T,tbar,cnR,Nlim,Tx)
      errors[i,m] = FSerror

  for i in range(1,max_L):
    Nlim = 2**i
    for m in range(1,max_m):
      ratios[i,m] = errors[i-1,m]/errors[i,m]

  print '{\\tiny\\begin{verbatim}'
  print sys.argv
  print '\\end{verbatim}}'
  print '\\begin{tabular}{|l|',
  for m in range(1,max_m):
    print 'l|',
  print '}\hline'

  print '$L\\mid m$',
  for m in range(1,max_m):
    print "&    ${0:7d}$         ".format(m),
  print "  \\\\\\hline"
  for i in range(0,max_L):
    Nlim = 2**i
    print "${0:6d}$ ".format(Nlim), 
    for m in range(1,max_m):
      #(sign, digits, exponent) = Decimal(errors[i,m]).as_tuple()
      #print sign, digits, exponent
      if i>0:
        print "& $"+str('{:8.1e}'.format(errors[i,m])).replace("e","(")+")_{{({:5.1f})}}$".format(ratios[i,m]),
      else:
        print "& $"+str('{:8.1e}'.format(errors[i,m])).replace("e","(")+")$          ",
    print '  \\\\'
  print '\\hline\\end{tabular}'
  
  # create graphics
  marker = itertools.cycle((',', '+', '.', 'o', '*')) 
  Nlims = [float(2**i) for i in range(0,max_L)]
  for m in range(1,max_m):
    plt.loglog(Nlims, errors[:,m], \
    marker = marker.next(), linestyle='-', color='black', linewidth=4, markersize=18, \
    label=r'$'+str(m)+'$' )
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  plt.rcParams.update({'font.size': 18})
  plt.xlabel(r'$N$',fontsize=25)
  plt.ylabel(r'$\Vert(I-F_N)\varphi\Vert_{L_\infty(T_1,T)}$',fontsize=25)
  plt.tick_params(labelsize=25)
  plt.legend(loc="lower left")
  plt.title(r'Fourier Series error (approximate)',fontsize=25)
  plt.tight_layout()
  plt.savefig('FSerror.png', format='png', dpi=750)
  plt.savefig('FSerror.jpg', format='jpg')
  plt.savefig('FSerror.eps', format='eps', dpi=1000)
  if gfx > 0:
    plt.show(block=False); time.sleep(gfx); plt.close()
  elif gfx < 0:
    plt.show()
  plt.clf()
  
