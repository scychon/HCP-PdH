from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from drudenoseplugin import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import mdtraj as md
import multiprocessing
from multiprocessing import Process, Queue


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits import mplot3d
from matplotlib.ticker import FormatStrFormatter
from matplotlib.font_manager import FontProperties
import os
import scipy.spatial as ss
import numpy as np

font = {'weight' : 'bold',
        'size'   : 22}
plt.rc('font',**font)
fontP = FontProperties()
#fontP.set_family('normal')
fontP.set_weight('bold')
fontP.set_size('18')

plt.ion()
plt.show()
plt.clf()

temperature=300*kelvin

rc_pdpd = 5.35
rc_pdh = 4.95
rc_pdh = 5.35
rc_hh = 5.35
rs = 0.3
rho0_pd = 10.261

strdir = '/data/cson/PdH/init_coord/'
# Read nanoparticle postions
pdpos = np.genfromtxt(strdir+'nanoptl_all.xyz',skip_header=2,skip_footer=1)[:,1:4]
rcom = np.average(pdpos,axis=0)
npd = len(pdpos)

x=pdpos[:,0]
y=pdpos[:,1]
z=pdpos[:,2]

bins = np.arange(51)
centers=(bins[1:]+bins[:-1])/2

histpd,bins=np.histogram(np.linalg.norm(pdpos-rcom,axis=1),bins=bins)

hpos0 = np.genfromtxt('first.xyz',skip_header=2)[:,1:4]
rcomh0 = np.average(hpos0,axis=0)
nh = len(hpos0)

hpose = np.genfromtxt('last.xyz',skip_header=2)[:,1:4]

histh0,bins=np.histogram(np.linalg.norm(hpos0-rcom,axis=1),bins=bins)
histhe,bins=np.histogram(np.linalg.norm(hpose-rcom,axis=1),bins=bins)
plt.scatter(centers,histpd/(centers*centers),label='Pd')
plt.scatter(centers,histh0/(centers*centers),label='H0')
plt.scatter(centers,histhe/(centers*centers),label='He')
plt.legend()

plt.scatter(centers,(histhe-histpd)/(centers*centers),label='He')
plt.scatter(centers,(histh0-histpd)/(centers*centers),label='H0')


def fpdu(rho):
    rhop=rho/50
    fpduval = 295878.9003038662 *(rhop-0.20581955357385892)*(rhop-0.081228755904399)*rhop*(rhop+0.05298811034615951) \
   *(rhop*(rhop-2.4242616904962846) + 1.4791899886249564) * (rhop*(rhop-2.1376274623740064) + 1.2169215689822592) \
   *(rhop*(rhop-1.6486007989726832) + 0.8159825255339774) * (rhop*(rhop-1.0749204110338482) + 0.42007491336688396) \
   *(rhop*(rhop-0.5128056047933808) + 0.12468685331167456)
    return fpduval

fpdup0 = 0.2515463239137249
def fpd(rho):
    return fpdu(rho)-fpdup0*rho


def rhoa(r5):
    rhoaval = r5*(r5*(r5-2.7267629107325706) + 1.8716766113599643) * (r5*(r5-2.50290548851635) + 1.668549182690922) \
                *(r5*(r5-2.0924467509943674) + 1.3150372774478005) * (r5*(r5-1.564328475106985) + 0.8987511149780485) \
                *(r5*(r5-1.009780903403673) + 0.5124363774128722) * (r5*(r5-0.5304054524800665) + 0.2169886022464641) \
                *(r5*(r5-0.1356566408715063) + 0.035852347523891395)
    return rhoaval

def rc_scale(r,rc,rs):
    condcut = r<rc
    condscale = r>rc-rs
    condfull = r<=rc-rs
    scale = condfull+condcut*condscale*(1+np.cos(np.pi*(r-rc+rs)/rs))/2
    return scale

def rhoa_pd(r):
    r5=r/5
    rhoapdval = -0.02972698211669922 + r5*(0.6676807403564453 + r5*(-255.8965835571289 + r5*(14673.409149169922-(2.597301181336601e7)*rhoa(r5))))
    return rhoapdval*rc_scale(r,rc_pdpd,rs)


def fhu(rho):
    (ah, bh, ch, dh, eh) = (9.99780, 60.0155, 0.000197047, 1.18860, 0.0540638)
    rhoeh = rho+eh
    rhoehdh = rhoeh**dh
    fhuval = -ch*rhoehdh*((rhoeh**2)/(2+dh) - (ah+bh)/(1+dh)*rhoeh + (ah*bh)/dh)
    return fhuval

def fhup(rho):
    (ah, bh, ch, dh, eh) = (9.99780, 60.0155, 0.000197047, 1.18860, 0.0540638)
    rhoeh = rho+eh
    fhupval = -ch*(rhoeh**(1+dh) - (ah+bh)*rhoeh**dh + (ah*bh)*rhoeh**(dh-1))
    return fhupval

rho0_h=7.98909
fhup0=fhup(rho0_h)
fhup0=-0.0296601759808718

def fh(rho):
    return fhu(rho)-fhup0*rho

def rhoa_h(r):
    CH = 11.0025
    delta_h = 1.30927
    return  CH * np.exp(-delta_h*r)*rc_scale(r,rc_hh,rs)

def pi_pdu(r):
    r5=r/5
    pipdval = -79415.24035137112 *(r5-1.0699996145674568)*(r5-1.06015072612581) \
             *(r5-0.42433991011376526)*(r5+0.06169160085238687) \
             *(r5*(r5-2.0586473420376348) + 1.0683922574015199)*(r5*(r5-1.6696359816422877) + 0.7337878627470482) \
             *(r5*(r5-1.1690370066230809) + 0.3909805777737639)*(r5*(r5-0.2635598721249787) + 0.033551116514910245)
    return pipdval

def pi_pd(r):
    pival = pi_pdu(r)+2*fpdup0*rhoa_pd(r)
    return pival*rc_scale(r,rc_pdpd,rs)

def pi_pdh(r):
    (alpha,beta,r0,Dval) = (4.82613,2.13158,1.50964,0.2494540)
    pival = Dval*(beta*np.exp(-alpha*(r-r0))-alpha*np.exp(-beta*(r-r0)))
    return pival*rc_scale(r,rc_pdh,rs)

def pi_hh(r):
    (alpha,beta,r0,Dval) = (3.67263, 1.47797,2.51980,0.0661496)
    pival = Dval*(beta*np.exp(-alpha*(r-r0))-alpha*np.exp(-beta*(r-r0)))
    return pival*rc_scale(r,rc_hh,rs)

plt.clf()
xlist = np.arange(0.01,5,.01)
plt.plot(xlist,rhoa_pd(xlist),'-',label='Pd')
plt.plot(xlist,rhoa_h(xlist),'--',label='H')
plt.legend()
plt.xlim(0,5)
plt.ylim(0,5)
plt.yticks(np.arange(0,5.5,1))
plt.ylabel(r'$\rho$')
plt.xlabel(r'$r$ $(\AA)$')


plt.figure()
plt.clf()
xlist = np.arange(0.1,50,.1)
plt.plot(xlist,fpd(xlist),'-',label='Pd')
plt.plot(xlist,fh(xlist),'--',label='H')
plt.legend()
plt.xlim(0,50)
plt.ylim(-4,12)
plt.yticks(np.arange(-4,13,2))
plt.ylabel(r'$F$ $(eV)$')
plt.xlabel(r'$\rho$')

plt.figure()
plt.clf()
xlist = np.arange(0.01,5,.01)
plt.plot(xlist,pi_pd(xlist),'-',label='Pd')
plt.plot(xlist,pi_pdh(xlist),'--',label='PdH')
plt.plot(xlist,pi_hh(xlist),'--',label='H')
plt.legend()
plt.xlim(1,5)
plt.ylim(-1,1)
#plt.yticks(np.arange(-1,1,1))
plt.ylabel(r'$\phi (eV)$')
plt.xlabel(r'$r$ $(\AA)$')



def fcc_pdhoct(aval):
    pdpos = []
    hpos = []
    for i in [0,1]:
      for j in [0,1]:
        for k in [0,1]:
          xval = aval/2*np.asarray((i,j,k))
          if (i+j+k)%2==0:
            pdpos.append(xval)
          elif (i+j+k)%2==1:
            hpos.append(xval)
    return (pdpos,hpos)

def ener_fcc_pdhoct(aval,hfrac):
    (pdpos,hpos) = fcc_pdhoct(aval)
    rhopd = np.zeros(4)
    rhoh = np.zeros(4)
    ftot = fpdval = fhval = 0
    phitot = phi_pdpd = phi_pdh = phi_hh = 0
    for jj in range(len(pdpos)):
        pdval = pdpos[jj]
        hval = hpos[jj]
        for i in range(-2,4):
          for j in range(-2,4):
            for k in range(-2,4):
              xval = aval/2*np.asarray((i,j,k))
              drpd = np.linalg.norm(pdval-xval)
              drh = np.linalg.norm(hval-xval)
              if (i+j+k)%2==0:
                if drpd < 5.35 and drpd > 0.1:
                    rhopd[jj]+=rhoa_pd(drpd)
                    phi_pdpd += pi_pd(drpd)
                if drh < 5.35 and drh > 0.1:
                    rhoh[jj]+=rhoa_pd(drh)
                    phi_pdh += hfrac*pi_pdh(drh)
              elif (i+j+k)%2==1:
                if drpd < 5.35 and drpd > 0.1:
                    rhopd[jj]+=hfrac*rhoa_h(drpd)
                    phi_pdh += hfrac*pi_pdh(drpd)
                if drh < 5.35 and drh > 0.1:
                    rhoh[jj]+=hfrac*rhoa_h(drh)
                    phi_hh += hfrac**2 * pi_hh(drh)
        fpdval += fpd(rhopd[jj])
        fhval += hfrac*fh(rhoh[jj])
    ftot = fpdval + fhval
    phi_pdpd *= 0.5
    phi_pdh *= 0.5
    phi_hh *= 0.5
    phitot = phi_pdpd + phi_pdh + phi_hh
    Etot = ftot + phitot
    return (Etot, ftot, phitot, fpdval, fhval, phi_pdpd, phi_pdh, phi_hh, rhopd,rhoh)

def fcc_znblend(aval):
    pdpos = []
    hpos = []
    hdiff = aval/4*np.asarray((1,1,1))
    for i in [0,1]:
      for j in [0,1]:
        for k in [0,1]:
          xval = aval/2*np.asarray((i,j,k))
          if (i+j+k)%2==0:
            pdpos.append(xval)
            hpos.append(xval+hdiff)
    return (pdpos,hpos)

def ener_fcc_znblend(aval,hfrac):
    (pdpos,hpos) = fcc_znblend(aval)
    hdiff = aval/4*np.asarray((1,1,1))
    rhopd = np.zeros(4)
    rhoh = np.zeros(4)
    ftot = fpdval = fhval = 0
    phitot = phi_pdpd = phi_pdh = phi_hh = 0
    for jj in range(len(pdpos)):
        pdval = pdpos[jj]
        hval = hpos[jj]
        for i in range(-2,4):
          for j in range(-2,4):
            for k in range(-2,4):
              if (i+j+k)%2==0:
                xval = aval/2*np.asarray((i,j,k))
                drpdpd = np.linalg.norm(pdval-xval)
                drhpd = np.linalg.norm(hval-xval)
                drpdh = np.linalg.norm(pdval-xval-hdiff)
                drhh = np.linalg.norm(hval-xval-hdiff)
                if drpdpd < 5.35 and drpdpd > 0.1:
                    rhopd[jj]+=rhoa_pd(drpdpd)
                    phi_pdpd += pi_pd(drpdpd)
                if drhpd < 5.35 and drhpd > 0.1:
                    rhoh[jj]+=rhoa_pd(drhpd)
                    phi_pdh += hfrac*pi_pdh(drhpd)
                if drpdh < 5.35 and drpdh > 0.1:
                    rhopd[jj]+=hfrac*rhoa_h(drpdh)
                    phi_pdh += hfrac*pi_pdh(drpdh)
                if drhh < 5.35 and drhh > 0.1:
                    rhoh[jj]+=hfrac*rhoa_h(drhh)
                    phi_hh += hfrac**2 * pi_hh(drhh)
        fpdval += fpd(rhopd[jj])
        fhval += hfrac*fh(rhoh[jj])
    ftot = fpdval + fhval
    phi_pdpd *= 0.5
    phi_pdh *= 0.5
    phi_hh *= 0.5
    phitot = phi_pdpd + phi_pdh + phi_hh
    Etot = ftot + phitot
    return (Etot, ftot, phitot, fpdval, fhval, phi_pdpd, phi_pdh, phi_hh, rhopd,rhoh)

###################
# HCP Lattice model
###################
vec0=np.array((1,0,0))
vec1 = np.array((np.cos(np.pi*2/3),np.sin(np.pi*2/3),0))
vec2 = np.array((0,0,1))

def ener_hcp(aval,cval,hfrac,xlist):
    (x0,x1,x2,x3)= xlist
    pdpos = [x0,x1]
    hpos = [x2,x3]
    rhopd = np.zeros(2)
    rhoh = np.zeros(2)
    ftot = fpdval = fhval = 0
    phitot = phi_pdpd = phi_pdh = phi_hh = 0
    for jj in range(len(pdpos)):
        pdval = pdpos[jj]
        hval = hpos[jj]
        for i in range(-4,5):
          for j in range(-4,5):
            for k in range(-4,5):
              xorg = (vec0*i + vec1*j)*aval + (vec2*k)*cval
              for kk in range(2):
                xval = xorg+pdpos[kk]
                drpdpd = np.linalg.norm(pdval-xval)
                drhpd = np.linalg.norm(hval-xval)
                xval = xorg+hpos[kk]
                drpdh = np.linalg.norm(pdval-xval)
                drhh = np.linalg.norm(hval-xval)
                if drpdpd < 5.35 and drpdpd > 0.1:
                    rhopd[jj]+=rhoa_pd(drpdpd)
                    phi_pdpd += pi_pd(drpdpd)
                if drhpd < 5.35 and drhpd > 0.1:
                    rhoh[jj]+=rhoa_pd(drhpd)
                    phi_pdh += hfrac*pi_pdh(drhpd)
                if drpdh < 5.35 and drpdh > 0.1:
                    rhopd[jj]+=hfrac*rhoa_h(drpdh)
                    phi_pdh += hfrac*pi_pdh(drpdh)
                if drhh < 5.35 and drhh > 0.1:
                    rhoh[jj]+=hfrac*rhoa_h(drhh)
                    phi_hh += hfrac**2 * pi_hh(drhh)
        fpdval += fpd(rhopd[jj])
        fhval += hfrac*fh(rhoh[jj])
    ftot = fpdval + fhval
    phi_pdpd *= 0.5
    phi_pdh *= 0.5
    phi_hh *= 0.5
    phitot = phi_pdpd + phi_pdh + phi_hh
    Etot = ftot + phitot
    return (Etot, ftot, phitot, fpdval, fhval, phi_pdpd, phi_pdh, phi_hh, rhopd,rhoh)

def ener_hcp_oct(aval,cval,hfrac):
    x0 = np.array((0,0,0))
    x1 = np.array((0,np.sqrt(1/3)*aval,.5*cval))
    x2 = np.array((aval/2,np.sqrt(1/12)*aval,cval/4))
    x3 = np.array((aval/2,np.sqrt(1/12)*aval,cval*3/4))
    poslist = (x0,x1,x2,x3)
    return ener_hcp(aval,cval,hfrac,poslist)

def ener_hcp_tetra(aval,cval,hfrac):
    x0 = np.array((0,0,0))
    x1 = np.array((0,np.sqrt(1/3)*aval,.5*cval))
    x2 = np.array((0,0,cval*5/8))
    x3 = np.array((0,np.sqrt(1/3)*aval,cval/8))
    poslist = (x0,x1,x2,x3)
    return ener_hcp(aval,cval,hfrac,poslist)

####################
#run density calc
plt.figure()
xlist = np.linspace(3.887,4.5,100)
xlist = np.arange(4.0,4.3,0.005)
covera_hcp = sqrt(8/3)
xlist = np.linspace(3.5,4.5,100)
listfccoc = np.zeros((len(xlist),8))
listfccte = np.zeros((len(xlist),8))
listhcpoc = np.zeros((len(xlist),8))
listhcpte = np.zeros((len(xlist),8))
fhval = 0.2
for i in range(len(xlist)):
    listfccoc[i] = ener_fcc_pdhoct(xlist[i],fhval)[0:8]
    listfccte[i] = ener_fcc_znblend(xlist[i],fhval)[0:8]
    listhcpoc[i] = ener_hcp_oct(xlist[i]/sqrt(2),covera_hcp*xlist[i]/sqrt(2),fhval)[0:8]
    listhcpte[i] = ener_hcp_tetra(xlist[i]/sqrt(2),covera_hcp*xlist[i]/sqrt(2),fhval)[0:8]

plt.clf()
plt.plot(xlist/sqrt(2),listfccoc[:,0]/4,label='FCC_OC')
plt.plot(xlist/sqrt(2),listfccte[:,0]/4,label='FCC_TE')
plt.plot(xlist/sqrt(2),listhcpoc[:,0]/2,label='HCP_OC')
plt.plot(xlist/sqrt(2),listhcpte[:,0]/2,label='HCP_TE')
plt.legend()
plt.xlabel(r'$a (\AA)$')
plt.ylabel('E (eV/Pd atom)')

plt.ylim(-6.32,-5.8)
plt.xlim(2.8,3.15)
plt.title('PdH')
plt.savefig('latticeE.png',bbox_inches='tight')

plt.ylim(-4.45,-3.9)
plt.xlim(2.6,2.95)
plt.title('PdH.2')
plt.savefig('latticeE_f.2.png',bbox_inches='tight')






#################
# Pd hcp
aval = 2.728
cval = 4.563
x0 = (vec0*1/3 + vec1*2/3)*aval + (vec2*.25)*cval
x1 = (vec0*2/3 + vec1*1/3)*aval + (vec2*.75)*cval
x0, x1

xposlist = []
xposlist.append(x0)
xposlist.append(x1)
fname = 'hcp_pd2.xyz'

# PdH_oct.m3
aval = 2.871
cval = 4.7846
x0 = (vec0*1/3 + vec1*2/3)*aval + (vec2*.25)*cval
x1 = (vec0*2/3 + vec1*1/3)*aval + (vec2*.75)*cval
x2 = 0.050152318*vec2*cval
x3 = 0.550152318*vec2*cval
x0,x1
x2,x3

xposlist = []
xposlist.append(x0)
xposlist.append(x1)
fname = 'hcp_pdh_oct.xyz'

# Pd8H2
aval = 5.5274
cval = 4.6279
vec00 = np.asarray([5.523499668,	0.006836391,	0])
vec01 = np.asarray([-2.755829348,	4.786909226,	0])
vec02 = np.asarray([0,0,4.6279])

vecmat=np.asarray(
[[0.163524232,	0.336388265,	0.260063987],
[0.167043094,	0.832956906,	0.259146747],
[0.667003888,	0.332996112,	0.261012546],
[0.663611735,	0.836475768,	0.260063987],
[0.336475768,	0.163611735,	0.760063987],
[0.332956906,	0.667043094,	0.759146747],
[0.832996112,	0.167003888,	0.761012546],
[0.836388265,	0.663524232,	0.760063987]])

xposlist = []
for vecscale in vecmat:
    xposlist.append(vec00*vecscale[0]+vec01*vecscale[1]+vec02*vecscale[2])

fname = 'hcp_pd8h2.xyz'

#Pd8H4
aval = 5.6064
cval = 4.6673
vec00 = np.asarray([5.605975643,	0.000676247,	0.000])
vec01 = np.asarray([-2.802402174,	4.855255441,	0.000])
vec02 = np.asarray([0.000,	0.000,	4.667255892])

vecmat=np.asarray(
[0.167768684,	0.33223138,	0.263586258],
[0.167765232,	0.832234768,	0.269859015],
[0.667765237,	0.332234763,	0.269858961],
[0.66776862,    0.832231316,	0.263586258],
[0.332234653,	0.167765211,	0.769859053],
[0.332231438,	0.667768562,	0.763586275],
[0.832231424,	0.167768576,	0.763586295],
[0.832234789,	0.667765347,	0.769859053]])

xposlist = []
for vecscale in vecmat:
    xposlist.append(vec00*vecscale[0]+vec01*vecscale[1]+vec02*vecscale[2])

fname = 'hcp_pd8h4.xyz'

# Pd fcc
aval = 3.887
fname = 'fcc_pd4.xyz'

# PdH fcc
aval = 4.104
fname = 'fcc_pdh.xyz'

def gen_fcc(aval,ncp):
    xoutlist = []
    mm = -2*ncp
    ll = 2*(ncp+1)
    for i in range(mm,ll):
      for j in range(mm,ll):
        for k in range(mm,ll):
          if (i+j+k)%2==0:
            xval = aval/2*np.asarray((i,j,k))
            xoutlist.append(xval)
    return xoutlist

#################
# run energy calculation for FCC
#################

listaval = [3.887, 4.104]

pdlist = []
hlist = []
#octhedral
for aval in listaval:
    pdpos = []
    hpos = []
    for i in [0,1]:
      for j in [0,1]:
        for k in [0,1]:
          xval = aval/2*np.asarray((i,j,k))
          if (i+j+k)%2==0:
            pdpos.append(xval)
          elif (i+j+k)%2==1:
            hpos.append(xval)
    pdlist.append(pdpos)
    hlist.append(hpos)



for ii in range(len(pdlist)):
    aval = listaval[ii]
    pdpos = pdlist[ii]
    hpos = hlist[ii]
    rhopd = np.zeros(4)
    rhoh = np.zeros(4)
    ftot = fpdval = fhval = 0
    phitot = phi_pdpd = phi_pdh = phi_hh = 0
    for jj in range(len(pdpos)):
        pdval = pdpos[jj]
        hval = hpos[jj]
        for i in range(-2,4):
          for j in range(-2,4):
            for k in range(-2,4):
              xval = aval/2*np.asarray((i,j,k))
              drpd = np.linalg.norm(pdval-xval)
              drh = np.linalg.norm(hval-xval)
              if (i+j+k)%2==0:
                if drpd < 5.35:
                    rhopd[jj]+=rhoa_pd(drpd)
                    phi_pdpd += pi_pd(drpd)
                if drh < 5.35:
                    rhoh[jj]+=rhoa_pd(drh)
                    phi_pdh += pi_pdh(drh)
              elif (i+j+k)%2==1:
                if drpd < 5.35:
                    rhopd[jj]+=rhoa_h(drpd)
                    phi_pdh += pi_pdh(drpd)
                if drh < 5.35:
                    rhoh[jj]+=rhoa_h(drh)
                    phi_hh += pi_hh(drh)
        fpdval += fpd(rhopd[jj])
        fhval += fh(rhoh[jj])
    ftot = fpdval + fhval
    phitot = phi_pdpd + phi_pdh + phi_hh
    Etot = ftot + phitot
    print('aval, Etot, ftot, phitot, fpdval, fhval, phi_pdpd, phi_pdh, phi_hh')
    print(aval, Etot, ftot, phitot, fpdval, fhval, phi_pdpd, phi_pdh, phi_hh)
    print(rhoh)









# get min & max r positions
pdxyz = np.genfromtxt('PdH_xyz_AddDomType.xyz',skip_header=2,skip_footer=1)[:,1:]
pdpos = pdxyz[:,:3]
rmin = np.min(pdpos,axis=0)-3
rmax = np.max(pdpos,axis=0)+3
npd = len(pdpos)

nz=int(np.ceil((rmax-rmin)[2]/cval))
ny=int(np.ceil((rmax-rmin)[1]/aval/vec1[1]))
nx=int(np.ceil((-vec1*aval*ny+(rmax-rmin))[0]/aval))

xoutlist = []
for i in range(nx):
  for j in range(ny):
    for k in range(nz):
      xval =  rmin+ aval*vec0*i + aval*vec1*j + cval*vec2*k-1
      if (xval[0]>=rmin[0] and xval[0] <= rmax[0] and xval[1]<=rmax[1]) :
        for xi in xposlist:
            xpos = xval + xi
            xoutlist.append(' Pd%10.3f000%10.3f000%10.3f000\n' % (xpos[0],xpos[1],xpos[2]))


for i in range(npd):
    xpos = pdpos[i]
    xoutlist.append(' Po%10.3f000%10.3f000%10.3f000\n' % (xpos[0],xpos[1],xpos[2]))


fout = open(fname, "w+")
fout.write(str(len(xoutlist))+'\n')
fout.write('hcp_pd2\n')
for xoutline in xoutlist:
    ival = fout.write(xoutline)


fout.write(' -1')
fout.close()










# Generate Pd position array

for i in range(npd-1):
    for j in range(i+1,npd):
        x0 = pdpos[i]
        x1 = pdpos[j]

rhopd = np.zeros((npd,3))
rhopdgrid = np.zeros((ngrid,3))



#This derivative function is obtaied by
#import sympy as sym
#rhop = sym.Symbol('rhop')
#sym.diff(fpdu(rho))
def fpdup(rho):
    rhop=rho/50
    fpdupval = rhop*(rhop - 0.081228755904399)*(rhop + 0.0529881103461595)*(2*rhop - 2.42426169049628)*(295878.900303866*rhop - 60897.663172466)*(rhop*(rhop - 2.13762746237401) + 1.21692156898226)*(rhop*(rhop - 1.64860079897268) + 0.815982525533977)*(rhop*(rhop - 1.07492041103385) + 0.420074913366884)*(rhop*(rhop - 0.512805604793381) + 0.124686853311675) + rhop*(rhop - 0.081228755904399)*(rhop + 0.0529881103461595)*(2*rhop - 2.13762746237401)*(295878.900303866*rhop - 60897.663172466)*(rhop*(rhop - 2.42426169049628) + 1.47918998862496)*(rhop*(rhop - 1.64860079897268) + 0.815982525533977)*(rhop*(rhop - 1.07492041103385) + 0.420074913366884)*(rhop*(rhop - 0.512805604793381) + 0.124686853311675) + rhop*(rhop - 0.081228755904399)*(rhop + 0.0529881103461595)*(2*rhop - 1.64860079897268)*(295878.900303866*rhop - 60897.663172466)*(rhop*(rhop - 2.42426169049628) + 1.47918998862496)*(rhop*(rhop - 2.13762746237401) + 1.21692156898226)*(rhop*(rhop - 1.07492041103385) + 0.420074913366884)*(rhop*(rhop - 0.512805604793381) + 0.124686853311675) + rhop*(rhop - 0.081228755904399)*(rhop + 0.0529881103461595)*(2*rhop - 1.07492041103385)*(295878.900303866*rhop - 60897.663172466)*(rhop*(rhop - 2.42426169049628) + 1.47918998862496)*(rhop*(rhop - 2.13762746237401) + 1.21692156898226)*(rhop*(rhop - 1.64860079897268) + 0.815982525533977)*(rhop*(rhop - 0.512805604793381) + 0.124686853311675) + rhop*(rhop - 0.081228755904399)*(rhop + 0.0529881103461595)*(2*rhop - 0.512805604793381)*(295878.900303866*rhop - 60897.663172466)*(rhop*(rhop - 2.42426169049628) + 1.47918998862496)*(rhop*(rhop - 2.13762746237401) + 1.21692156898226)*(rhop*(rhop - 1.64860079897268) + 0.815982525533977)*(rhop*(rhop - 1.07492041103385) + 0.420074913366884) + 295878.900303866*rhop*(rhop - 0.081228755904399)*(rhop + 0.0529881103461595)*(rhop*(rhop - 2.42426169049628) + 1.47918998862496)*(rhop*(rhop - 2.13762746237401) + 1.21692156898226)*(rhop*(rhop - 1.64860079897268) + 0.815982525533977)*(rhop*(rhop - 1.07492041103385) + 0.420074913366884)*(rhop*(rhop - 0.512805604793381) + 0.124686853311675) + rhop*(rhop - 0.081228755904399)*(295878.900303866*rhop - 60897.663172466)*(rhop*(rhop - 2.42426169049628) + 1.47918998862496)*(rhop*(rhop - 2.13762746237401) + 1.21692156898226)*(rhop*(rhop - 1.64860079897268) + 0.815982525533977)*(rhop*(rhop - 1.07492041103385) + 0.420074913366884)*(rhop*(rhop - 0.512805604793381) + 0.124686853311675) + rhop*(rhop + 0.0529881103461595)*(295878.900303866*rhop - 60897.663172466)*(rhop*(rhop - 2.42426169049628) + 1.47918998862496)*(rhop*(rhop - 2.13762746237401) + 1.21692156898226)*(rhop*(rhop - 1.64860079897268) + 0.815982525533977)*(rhop*(rhop - 1.07492041103385) + 0.420074913366884)*(rhop*(rhop - 0.512805604793381) + 0.124686853311675) + (rhop - 0.081228755904399)*(rhop + 0.0529881103461595)*(295878.900303866*rhop - 60897.663172466)*(rhop*(rhop - 2.42426169049628) + 1.47918998862496)*(rhop*(rhop - 2.13762746237401) + 1.21692156898226)*(rhop*(rhop - 1.64860079897268) + 0.815982525533977)*(rhop*(rhop - 1.07492041103385) + 0.420074913366884)*(rhop*(rhop - 0.512805604793381) + 0.124686853311675)
    fpdupval /=50
    return fpdupval


#fpdup0 = fpdup(rho0_pd)

def gen_fcc(nx,ny,nz):
    matfcc = np.zeros((nx*2+1,ny*2+1,nz*2+1),dtype=int)
    for i in range(nx*2+1):
        for j in range(ny*2+1):
            for k in range(nz*2+1):
                if (i+j+k)%2==0:
                    matfcc[i,j,k]=1
    return matfcc


####
# hcp oct xpos list
    x0 = (vec0*1/3 + vec1*2/3)*aval + (vec2*.25)*cval
    x1 = (vec0*2/3 + vec1*1/3)*aval + (vec2*.75)*cval
    x2 = 0*vec2*cval
    x3 = (vec2*.5)*cval



# get the number of CPUs and assign it as the number of
# processes to create
num_proc = multiprocessing.cpu_count()
print('You have {0:1d} CPUs'.format(num_proc))

# Create the processes
p_list=[]
for i in range(1,num_proc+1):
    p = Process(target=f, name='Process'+str(i), args=(i,))
    p_list.append(p)
    print('Process:: '+p.name)
    p.start()
    print('Was assigned PID:: '+p.pid)

# Wait for all the processes to finish
for p in p_list:
    p.join()

# Create the processes
p_list=[]

# Create the Queue which will have the partial products
product=Queue()

for i in range(1,num_proc+1):
    # A pair of queues per process for the two arrays
    aq = Queue()
    bq = Queue()

    # push the chunks into the queue
    aq.put(a[(i-1)*chunk:i*chunk])
    bq.put(b[(i-1)*chunk:i*chunk])

    # create the process
    p = Process(target=add, args=(aq,bq,product))
    p.start()
    p.join()

# collect the individual sums
items=[]
for i in range(product.qsize()):
    items.append(product.get())

# final product: sum of individual products
print "Dot product:: ",sum(items)
