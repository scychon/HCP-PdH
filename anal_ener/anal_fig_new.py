from __future__ import print_function
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
from matplotlib.font_manager import FontProperties
import os
import multiprocessing
from multiprocessing import Process, Queue

font = {'weight' : 'bold',
        'size'   : 18}
plt.rc('font',**font)
fontP = FontProperties()
#fontP.set_family('normal')
fontP.set_weight('bold')
fontP.set_size('18')

plt.ion()
plt.show()
plt.clf()
listptl = ['anal_nanoptl/','anal_hcp_pdh/','anal_fcc_pdh/']
listlabel = ['AET', 'HCP', 'FCC']
listcolor = ['k','r','b']

listr = [6,8,10,12,15,17,20,22,25,27,30]
listx = [6,8,10,12,15,17,20,22,25,27,30,39]
listptl = ['anal_hcp_pdh/','anal_fcc_pdh/']
listlabel = ['HCP', 'FCC']
listcolor = ['r','b']
strsymbol = ['--^','-.x']
listh = ['.2','.25','.3','.4','.5','.6','.75','.8','1']
listhval = np.asarray(listh,dtype='float')

########################
# Energy vs particle radius
########################
hval='1'
plt.clf()
 #istep, nPdAtoms,  nHAtoms, Etot, E_pd, E_h, ftot, ftot_pd, ftot_h, rho_pd_avg, rho_h_avg, rho_pdh_avg, rho_hpd_avg, phitot, phitot_pdpd, phitot_pdh, phitot_hpd, phitot_hh

listptl = ['anal_nanoptl/','anal_hcp_pdh/','anal_fcc_pdh/']
listlabel = ['AET', 'HCP', 'FCC']
listcolor = ['k','r','b']
strsymbol = ['-o','-s','-^']

listptl = ['anal_hcp_pdh/','anal_fcc_pdh/']
listlabel = ['HCP', 'FCC']
listcolor = ['r','b']
strsymbol = ['-s','-^']
rscale = [18.84359132,18.69217742]

enersum = np.zeros((len(listptl),len(listr)+1,len(listh),17))       # nH

for i in range(len(listptl)):
    for j in range(len(listr)):
      for k in range(len(listh)):
        hval = listh[k]
        strener = '../htraj/'+listptl[i]+'r'+str(listr[j])+'_h'+hval+'/ener.log'
        #enersum[i,j] = np.average(np.genfromtxt(strener,skip_header=1)[-100:],axis=0)[1:]
        enersum[i,j,k] = np.genfromtxt(strener,skip_header=1)[-1][1:]
    strener = '../htraj/'+listptl[i]+'h1/ener.log'
    enersum[i,-1,-1] = np.genfromtxt(strener,skip_header=1)[-1][1:]

plt.figure()
for k in range(len(listh)):
  plt.clf()
  enerdata = np.zeros((len(listptl)*2,11))
  for i in range(len(listptl)):
    #plt.plot(listr, (enersum[i,:11,k,2])/enersum[i,:11,k,0],strsymbol[i],label=listlabel[i])
    plt.plot(enersum[i,:11,k,0]**(1./3)*30/rscale[i], (enersum[i,:11,k,2])/enersum[i,:11,k,0],strsymbol[i],label=listlabel[i],color=listcolor[i])
    #plt.plot(listr, (enersum[i,:11,k,2])/enersum[i,:11,k,0],strsymbol[i],color=listcolor[i],label=listlabel[i])
    #plt.xscale('log')
    #plt.xlabel(r'$N_{Pd}$')
    plt.xlabel(r'Particle radius ($\AA$)')
    plt.ylabel('energy per PdH (eV)')
    plt.legend( prop=fontP,handletextpad=0.1)
    enerdata[i*2,:] = enersum[i,:11,k,0]**(1./3)*30/rscale[i]
    enerdata[i*2+1:,:] = enersum[i,:11,k,2]/enersum[i,:11,k,0]
  np.savetxt('energy_r_h'+listh[k]+'.dat',enerdata.T)
  plt.savefig('energy_n_vs_r_h'+listh[k]+'.png',bbox_inches='tight')




plt.clf()
for i in range(11):
  for k in range(len(listh)):
    if (k==8 and i==0):
        plt.scatter(enersum[0,i,k,0]**(1./3)*30/rscale[0], float(listh[k]),marker='s',color='r',label='HCP')
    elif (k>5 and i==0):
        plt.scatter(enersum[0,i,k,0]**(1./3)*30/rscale[0], float(listh[k]),marker='s',color='r')
    elif (i==10 and k==0):
        plt.scatter(enersum[0,i,k,0]**(1./3)*30/rscale[0], float(listh[k]),marker='^',color='b',label='FCC')
    elif (i>5 or k>5):
        plt.scatter(enersum[0,i,k,0]**(1./3)*30/rscale[0], float(listh[k]),marker='^',color='b')
    elif (i>3 and k==5):
        plt.scatter(enersum[0,i,k,0]**(1./3)*30/rscale[0], float(listh[k]),marker='^',color='b')


plt.legend( prop=fontP,handletextpad=0.1)
plt.xlabel(r'Particle radius ($\AA$)')
plt.ylabel('H content (x)')
plt.savefig('phase.png',bbox_inches='tight')
    elif (i<2):
        plt.scatter(enersum[0,i,k,0]**(1./3)*30/rscale[0], float(listh[k]),marker='^',color='b')

    else:
        plt.scatter(enersum[0,i,k,0]**(1./3)*30/rscale[0], float(listh[k]),marker='^',color='b')

k=8
plt.clf()
for i in range(len(listptl)):
    #plt.plot(listr, (enersum[i,:11,k,2])/enersum[i,:11,k,0],strsymbol[i],label=listlabel[i])
    plt.plot(enersum[i,:11,k,0]**(1./3)*30/rscale[i], (enersum[i,:11,k,2])/enersum[i,:11,k,0],strsymbol[i],color=listcolor[i],label=listlabel[i])
    #plt.xscale('log')
    plt.xlabel(r'Particle radius ($\AA$)')
    plt.ylabel('Energy per PdH (eV)')
    plt.legend( prop=fontP,handletextpad=0.1)

plt.title(r'$\rm{PdH}_{1}$')
plt.savefig('energy_n_vs_r_pdh.png',bbox_inches='tight')


for k in range(len(listh)):
  plt.figure()

k=0
for i in range(len(listptl)):
    plt.plot(enersum[i,:11,0,0], (enersum[i,:11,k,2])/enersum[i,:11,k,0],strsymbol[i],label=listlabel[i],color=listcolor[i])
    plt.xscale('log')
    plt.xlabel(r'$N_{Pd}$')
    plt.ylabel('energy per PdH (eV)')
    plt.legend( prop=fontP,handletextpad=0.1)

plt.title(r'$\rm{PdH}_{.2}$')
plt.savefig('energy_n_vs_r_pdh.2.png',bbox_inches='tight')

k=8
plt.clf()
plt.title(r'$\rm{PdH}_{1}$')
plt.savefig('energy_n_vs_r_pdh.png',bbox_inches='tight')

    #plt.plot(listr, (enersum[i,:11,k,2])/enersum[i,:11,k,0],strsymbol[i],color=listcolor[i],label=listlabel[i])


plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel('energy per PdH (eV)')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('energy_r.png',bbox_inches='tight')

enerdata = np.zeros((len(listptl)+1,11))
enerdata[0,:] = listr
enerdata[1:,:] = enersum[:,:11,2]/enersum[:,:11,0]
np.savetxt('energy_r.dat',enerdata.T)

plt.clf()
j=0
for i in range(len(listptl)):
    plt.plot(listhval, (enersum[i,j,:,2])/enersum[i,j,:,0],strsymbol[i],label=listlabel[i])

plt.clf()
for k in range(len(listh)):
  for i in range(len(listptl)):
    plt.plot(listr, (enersum[i,:11,k,2])/enersum[i,:11,k,0],strsymbol[i],label=listlabel[i])

plt.clf()
for j in range(len(listr)):
    plt.plot(listhval, (enersum[1,j,:,2]/enersum[1,j,:,0]-enersum[0,j,:,2]/enersum[0,j,:,0]),label=listr[j])

plt.legend()


plt.clf()
strsymbol = ['-o','-s','-^']
for i in range(len(listptl)):
    plt.plot(listr, (enersum[i,:11,3]),strsymbol[i],color=listcolor[i],label=listlabel[i])


plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, enersum[i,:11,3],'--x',color=listcolor[i],label=listlabel[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel('energy per PdH (eV)')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('energy_r_pd.png',bbox_inches='tight')

plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, enersum[i,:11,4],'-.^',color=listcolor[i],label=listlabel[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel('energy per PdH (eV)')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('energy_r_h.png',bbox_inches='tight')

########################
# Energy difference for smallest nanoparticle
########################
plt.clf()
 #istep, nPdAtoms,  nHAtoms, Etot, E_pd, E_h, ftot, ftot_pd, ftot_h, rho_pd_avg, rho_h_avg, rho_pdh_avg, rho_hpd_avg, phitot, phitot_pdpd, phitot_pdh, phitot_hpd, phitot_hh

listcolor=['k','r','b','y','m']
listlabel=listh
strsymbol = ['-o','-s','-^', '--','-.']
enersum = np.zeros((len(listh),len(listptl),len(listr)+1,17))       # nH

for i in range(len(listh)):
  for iptl in range(len(listptl)):
    for j in range(len(listr)):
        strener = '../htraj/'+listptl[iptl]+'r'+str(listr[j])+'_h'+listh[i]+'/ener.log'
        enersum[i,iptl,j] = np.genfromtxt(strener,skip_header=1)[-1][1:]
    strener = '../htraj/'+listptl[iptl]+'h1/ener.log'
    enersum[i,iptl,-1] = np.genfromtxt(strener,skip_header=1)[-1][1:]
  plt.plot(listr, enersum[i,1,:11,4]-enersum[i,0,:11,4],strsymbol[i],color=listcolor[i],label=listlabel[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel('energy per PdH (eV)')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('energy_r_hdiff.png',bbox_inches='tight')

plt.figure()
plt.clf()
listhval = [float(s) for s in listh]
for i in range(len(listr)):
    plt.plot(listhval, enersum[:,1,i,4]-enersum[:,0,i,4],label=listr[i])

plt.legend( prop=fontP,handletextpad=0.1)


########################
# Hsite type vs particle radius
########################

 #istep, nPdAtoms, nHAtoms, nHType(1:6), nHPdGroupNN(1:3), nHPdGroupCoord(1:3)
listptl = ['anal_hcp_pdh/','anal_fcc_pdh/']
listlabel = ['HCP', 'FCC']
listcolor = ['r','b']
#strsymbol = ['--^','-.x']
strsymbol = ['-^','-s']
htypesum = np.zeros((len(listptl),len(listr),14))       # nH

for i in range(len(listptl)):
    for j in range(len(listr)):
        strhtype = '../htraj/'+listptl[i]+'r'+str(listr[j])+'_h1/htype.dat'
        htypesum[i,j] = np.genfromtxt(strhtype,skip_header=1)[-1][1:15]

plt.clf()
for i in range(len(listptl)):
    plt.plot(enersum[i,:11,k,0]**(1./3)*30/rscale[i], (htypesum[i,:,2]+htypesum[i,:,4]+htypesum[i,:,6])/htypesum[i,:,1],strsymbol[i],color=listcolor[i],label=listlabel[i])
#    plt.plot(listr, (htypesum[i,:,2]+htypesum[i,:,6])/htypesum[i,:,0],'-o',color=listcolor[i],label=listlabel[i])
#    plt.plot(enersum[i,:11,k,0]**(1./3)*30/rscale[i], (htypesum[i,:,2])/htypesum[i,:,0],'-.o',color=listcolor[i],label=listlabel[i])
#    plt.plot(listr, (htypesum[i,:,4])/htypesum[i,:,0],'-.x',color=listcolor[i])
#    plt.plot(listr, htypesum[i,:,3]/htypesum[i,:,0],'--^',color=listcolor[i])

plt.xlabel(r'Particle radius ($\AA$)')
plt.ylim(0,0.85)
plt.legend( prop=fontP,handletextpad=0.1)
plt.ylabel(r'Highly coordinated H fraction')
plt.title(r'$\rm{PdH}_{1}$')
plt.savefig('prob_ocall_new.png',bbox_inches='tight')
plt.savefig('prob_ocall.png',bbox_inches='tight')

plt.clf()
for i in range(len(listptl)):
    plt.plot(enersum[i,:11,k,0]**(1./3)*30/rscale[i], (htypesum[i,:,2])/htypesum[i,:,1],strsymbol[i],color=listcolor[i],label=listlabel[i])

plt.ylabel(r'Octahedral H fraction')
plt.savefig('prob_oc_new.png',bbox_inches='tight')

plt.clf()
#strsymbol = ['--o','--s','--^']
strsymbol = ['--^','--s']
for i in range(len(listptl)):
    plt.plot(enersum[i,:11,k,0]**(1./3)*30/rscale[i], (htypesum[i,:,5])/htypesum[i,:,1],strsymbol[i],color=listcolor[i],label=listlabel[i])

plt.xlabel(r'Particle radius ($\AA$)')
plt.ylabel(r'Surface H fraction')
plt.ylim(0,0.5)
plt.legend( prop=fontP,handletextpad=0.1)
plt.title(r'$\rm{PdH}_{1}$')
plt.savefig('prob_surf_new.png',bbox_inches='tight')
plt.savefig('prob_surf.png',bbox_inches='tight')

plt.clf()
#strsymbol = ['-.o','-.s','-.^']
strsymbol = ['-.^','-.s']
for i in range(len(listptl)):
    plt.plot(enersum[i,:11,k,0]**(1./3)*30/rscale[i], (htypesum[i,:,3])/htypesum[i,:,1],strsymbol[i],color=listcolor[i],label=listlabel[i])

plt.xlabel(r'Particle radius ($\AA$)')
plt.ylabel(r'Tetrahedral H fraction')
plt.ylim(0,0.5)
plt.legend( prop=fontP,handletextpad=0.1)
plt.title(r'$\rm{PdH}_{1}$')
plt.savefig('prob_te_new.png',bbox_inches='tight')
plt.savefig('prob_te.png',bbox_inches='tight')

htypedata = np.zeros((len(listptl)+1,len(listr)))
htypedata[0,:] = listr
htypedata[1:,:] = (htypesum[:,:,2]+htypesum[:,:,4]+htypesum[:,:,6])/htypesum[:,:,0]
np.savetxt('htype_ocall.txt',htypedata.T)
htypedata[1:,:] = (htypesum[:,:,3])/htypesum[:,:,0]
np.savetxt('htype_te.txt',htypedata.T)
htypedata[1:,:] = (htypesum[:,:,5])/htypesum[:,:,0]
np.savetxt('htype_surf.txt',htypedata.T)

## h.2
htypesum = np.zeros((len(listptl),len(listr),14))       # nH

for i in range(len(listptl)):
    for j in range(len(listr)):
        strhtype = '../htraj/'+listptl[i]+'r'+str(listr[j])+'_h.2/htype.dat'
        htypesum[i,j] = np.genfromtxt(strhtype,skip_header=1)[-1][1:15]

plt.clf()
strsymbol = ['-^','-s']
for i in range(len(listptl)):
    plt.plot(enersum[i,:11,k,0]**(1./3)*30/rscale[i], (htypesum[i,:,2]+htypesum[i,:,4]+htypesum[i,:,6])/htypesum[i,:,1],strsymbol[i],color=listcolor[i],label=listlabel[i])

plt.xlabel(r'Particle radius ($\AA$)')
plt.ylim(0,0.85)
plt.legend( prop=fontP,handletextpad=0.1)
plt.ylabel(r'Highly coordinated H fraction')
plt.title(r'$\rm{PdH}_{.2}$')
plt.savefig('prob_ocall_h.2_new.png',bbox_inches='tight')

strsymbol = ['--^','--s']
plt.clf()
for i in range(len(listptl)):
    plt.plot(enersum[i,:11,k,0]**(1./3)*30/rscale[i], (htypesum[i,:,5])/htypesum[i,:,1],strsymbol[i],color=listcolor[i],label=listlabel[i])

plt.xlabel(r'Particle radius ($\AA$)')
plt.ylabel(r'Surface H fraction')
plt.ylim(0,0.5)
plt.legend( prop=fontP,handletextpad=0.1)
plt.title(r'$\rm{PdH}_{.2}$')
plt.savefig('prob_surf_h.2_new.png',bbox_inches='tight')

plt.clf()
strsymbol = ['-.^','-.s']
for i in range(len(listptl)):
    plt.plot(enersum[i,:11,k,0]**(1./3)*30/rscale[i], (htypesum[i,:,3])/htypesum[i,:,1],strsymbol[i],color=listcolor[i],label=listlabel[i])

plt.xlabel(r'Particle radius ($\AA$)')
plt.ylabel(r'Tetrahedral H fraction')
plt.ylim(0,0.8)
plt.legend( prop=fontP,handletextpad=0.1)
plt.title(r'$\rm{PdH}_{.2}$')
plt.savefig('prob_te_h.2_new.png',bbox_inches='tight')

htypedata = np.zeros((len(listptl)*2,len(listr)))
htypedata[0,:] = enersum[0,:11,k,0]**(1./3)*30/rscale[0]
htypedata[1,:] = enersum[1,:11,k,0]**(1./3)*30/rscale[1]
htypedata[2:,:] = (htypesum[:,:,2]+htypesum[:,:,4]+htypesum[:,:,6])/htypesum[:,:,0]
np.savetxt('htype_ocall_h.2.txt',htypedata.T)
htypedata[2:,:] = (htypesum[:,:,3])/htypesum[:,:,0]
np.savetxt('htype_te_h.2.txt',htypedata.T)
htypedata[2:,:] = (htypesum[:,:,5])/htypesum[:,:,0]
np.savetxt('htype_surf_h.2.txt',htypedata.T)

#####################
# Htype bar graph for r=6
#####################
# ==> anal_fcc_pdh/r6_h1/htype.dat <==
#       49991          55          55           7          20           2          26           0           0           0           0           0           0           0           0

#==> anal_hcp_pdh/r6_h1/htype.dat <==
#       49991          57          57          10           9          13          25           0           0           0           0           0           0           0           0

#==> anal_nanoptl/r6_h1/htype.dat <==
#       49991          57          57           5          14          11          26           1           0           0           0           0           0           0           0

#nhtype (1:Octahedral, 2:Tetrahedral, 3:InBetween, 4:Surface, 5:Crowded, 6:etc)
#   surf, te, f=5, oc, crowded
#   1-3, 4, 5, 6, more
listptl = ['anal_nanoptl/','anal_hcp_pdh/','anal_fcc_pdh/']
listlabel = ['AET', 'HCP', 'FCC']
listcolor = ['k','r','b']

hlist = np.zeros((3,3))
hlist[0] = [26,14,11+6]     # nanoptl
hlist[1] = [25,9,13+10]     # HCP
hlist[2] = [26,20,2+7]      # FCC

#xlabels = ['0-3\nsurface','4\n\ntetrahedral','5-8\nhighly coordinated']
xlabels = ['0-3\nsurface','4\n\ntetrahedral','5-6\nhighly coordinated']
xlabels = ['0-3\nsurface','4\ntetrahedral','5-6\nbipyramidal\noctahedral']

xlist = np.arange(3)

fig = plt.figure()
ax = fig.add_axes([0,0,1,1])

listbar = [-.25,0,.25]
for i in range(1,len(listptl)):
    ax.bar(xlist+listbar[i], hlist[i],color=listcolor[i],width=0.25,label=listlabel[i])


plt.xlabel(r'Number of coordinating Pd ($N_{Pd}^{coord}$)')
plt.xticks(xlist,xlabels)
plt.ylabel(r'Number of H atoms ($N_H$)')
plt.legend( prop=fontP,handletextpad=0.1)
plt.legend(bbox_to_anchor=(0.25,0.95), loc="lower left",ncol=3,prop=fontP,handletextpad=0.1)
plt.savefig('htype_bar_new.png',bbox_inches='tight')
plt.savefig('htype_bar.png',bbox_inches='tight')




################
# lattice energies
################
rc_pdpd = 5.35
rc_pdh = 4.95
rc_pdh = 5.35
rc_hh = 5.35
rs = 0.3
rho0_pd = 10.261

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

#################
# Run fraction & lattice constant vs E calculation
plt.figure()
xlist = np.linspace(3.887,4.078,5)
xlist1 = np.linspace(2.733,2.889,5)
flist = [0,.25,.5,.75,1]
#xlist = [3.887  , 3.93475, 3.9825 , 4.03025, 4.078, 4.185]
#xlist1 = [2.733, 2.772, 2.811, 2.85 , 2.889, 2.959 ]
covera_hcp = np.sqrt(8/3)
listfccoc = np.zeros((len(xlist),8))
listfccte = np.zeros((len(xlist),8))
listhcpoc = np.zeros((len(xlist),8))
listhcpte = np.zeros((len(xlist),8))
for i in range(len(xlist)):
    fhval = flist[i]
    listfccoc[i] = ener_fcc_pdhoct(xlist[i],fhval)[0:8]
    listfccte[i] = ener_fcc_znblend(xlist[i],fhval)[0:8]
    listhcpoc[i] = ener_hcp_oct(xlist1[i],covera_hcp*xlist1[i],fhval)[0:8]
    listhcpte[i] = ener_hcp_tetra(xlist1[i],covera_hcp*xlist1[i],fhval)[0:8]

plt.clf()
plt.plot(flist,(listfccoc[:,0]/4-listhcpoc[:,0]/2),'-^',color='b',label='FCC-O')
plt.plot(flist,(listfccte[:,0]/4-listhcpoc[:,0]/2),'--^',color='b',markerfacecolor="None",label='FCC-T')
plt.plot(flist,(listhcpoc[:,0]/2-listhcpoc[:,0]/2),'-s',color='r',label='HCP-O')
plt.plot(flist,(listhcpte[:,0]/2-listhcpoc[:,0]/2),'--s',color='r',markerfacecolor="None",label='HCP-T')

plt.scatter([0,1], [-0.0347,-0.0511],marker='^',color='k')
plt.scatter(1, 0.073,marker='s',facecolor="None",color='k')
#plt.legend( prop=fontP,handletextpad=0.1)
#plt.legend(bbox_to_anchor=(-0.15,1), loc="lower left",ncol=2,handleheight=2.4,labelspacing=-10,prop=fontP,handletextpad=0.1)
plt.legend(bbox_to_anchor=(0.1,0.85), loc="lower left",ncol=2,prop=fontP,handletextpad=0.1)

plt.xlabel(r'$x$')
plt.ylabel('relative E (eV/f.u.)')
plt.savefig('latticeEcomp.png',bbox_inches='tight')

latticedata = np.zeros((5,len(flist)))
latticedata[0,:] = flist
latticedata[1,:] = listfccoc[:,0]/4-listhcpoc[:,0]/2
latticedata[2,:] = listfccte[:,0]/4-listhcpoc[:,0]/2
latticedata[3,:] = listhcpoc[:,0]/2-listhcpoc[:,0]/2
latticedata[4,:] = listhcpte[:,0]/2-listhcpoc[:,0]/2
np.savetxt('lattice_f.txt',latticedata.T)

#################
# trajectory H map
#################

strdir = '/data/cson/PdH/init_coord/'
# Read nanoparticle postions
pdpos = np.genfromtxt(strdir+'nanoptl_all.xyz',skip_header=2,skip_footer=1)[:,1:4]
rcom = np.average(pdpos,axis=0)
npd = len(pdpos)

strdir = '/data/cson/PdH/htraj/nanoptl/h1/'
hpos0 = np.genfromtxt(strdir+'first.xyz',skip_header=2)[:,1:4]
hpos1 = np.genfromtxt(strdir+'last.xyz',skip_header=2)[:,1:4]
nh = len(hpos0)
nh1 = len(hpos1)

if (npd!=nh or nh != nh1):
    print('ERROR ! numper pf pd particles '+str(npd)+' does not match number of H particles '+str(nh)+' '+str(nh1)+' !')
    exit()



Res = 0.3479
SlInd=[212,206,199,193,186,179,173,166,160,153,147,140,133,127,120,114,107,101,94,88,81,75,68,61,55,48,42,36,30,23,16]
ndim = 230
ngrid = np.array((len(SlInd),ndim,ndim))

listhidxs = []
for j in range(len(SlInd)):
    listhidxs.append([])

hposidx = hpos1/Res
idxcut = rc_hh/Res
for i in range(nh):
    for j in range(len(SlInd)):
        idx = SlInd[j]
        if abs(hposidx[i,2]-idx) < idxcut:
            listhidxs[j].append(i)



pos0 = np.array((0.5,0.5,0.5))*Res
xpos = np.zeros((len(SlInd),230,230,3))
for i in range(len(SlInd)):
  for j in range(230):
    for k in range(230):
        xpos[i,j,k] = np.array((SlInd[i],j,k))*Res

rhoval = np.zeros((len(SlInd),230,230))
hpos = hpos1
for i in range(len(SlInd)):
  for j in range(230):
    for k in range(230):
      for idx in listhidxs[i]:
        dr = np.linalg.norm(hpos[idx]-xpos[i,j,k])
        rhoval[i,j,k] += rhoa_h(dr)
    print(str(j)+' th column')
  print(str(i)+' th line')

######################
# trajectory H map
#####################
strdir = '/data/cson/PdH/htraj/anal_nanoptl/h1_hmap/'
hmapdata = np.genfromtxt(strdir+'hmap.dat')
hmapdata = hmapdata[:,2:].reshape((230,230,31))
Res=0.3479
sig=10/Res

