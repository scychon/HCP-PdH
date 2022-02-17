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
listr = [6,8,10,12,15,17,20,22,25,27,30]
listx = [6,8,10,12,15,17,20,22,25,27,30,39]
listlabel = ['AET', 'HCP', 'FCC']
listcolor = ['k','r','b']
strsymbol = ['-o','--^','-.x']

plt.clf()
 # irbin, pdgroup, drPdNNAtoms(1:12)
drnnsum = []
strsymbol = ['-o','-','-']
for i in range(len(listptl)):
    strdrnn = '../htraj/'+listptl[i]+'h1_new/dr_pdnn.dat'
    drnn = np.genfromtxt(strdrnn,skip_header=1)
    drnnsum.append(drnn)
    idxs = np.where(np.max(drnn[:,2:],axis=1)<4)
    hist,bins = np.histogram(np.average(drnn[idxs[0],2:],axis=1),bins=30,normed=True)
    centers = (bins[1:]+bins[:-1])/2
    plt.plot(centers,hist,strsymbol[i],label=listlabel[i])
    if i==0:
        np.average(np.average(drnn[idxs[0],2:],axis=1)), np.std(np.average(drnn[idxs[0],2:],axis=1))

plt.ylim(0,7)
plt.xlabel(r'local lattice constant ($\AA$)')
plt.ylabel('probability density')
plt.legend( prop=fontP,handletextpad=0.1)

drnn = drnnsum[0]
listgrp = ['grp1','grp2','nogroup']
colorgrp = ['b','r','k']
listgrp = [1,2,-1]
listrbin = [0,20,30,40]
strsymgrp = ['-o','--x','-.^']
for i in range(3):
    idxs0 = np.where(drnn[:,1]==listgrp[i])[0]
    idxs1 = np.where(np.max(drnn[:,2:],axis=1)<4)[0]
    idxs = np.intersect1d(idxs0,idxs1)
    hist,bins = np.histogram(np.average(drnn[idxs,2:],axis=1),bins=30,normed=True)
    centers = (bins[1:]+bins[:-1])/2
    plt.plot(centers,hist,'--x',color=colorgrp[i])

plt.savefig('lattice_dist.png',bbox_inches='tight')

for i in [2]:
    idxs0 = np.where(drnn[:,1]==listgrp[i])[0]
    for j in range(3):
        idxs1 = np.where(drnn[:,0]>=listrbin[j])[0]
        idxs2 = np.where(drnn[:,0]<listrbin[j+1])[0]
        idxs1 = np.intersect1d(idxs1,idxs2)
        idxs2 = np.where(np.max(drnn[:,2:],axis=1)<4)[0]
        idxs1 = np.intersect1d(idxs1,idxs2)
        idxs = np.intersect1d(idxs0,idxs1)
        hist,bins = np.histogram(np.average(drnn[idxs,2:],axis=1),bins,normed=True)
        plt.plot(centers,hist,strsymgrp[j],color=colorgrp[i])

########################
# Energy vs particle radius
########################

plt.clf()
 #istep, nPdAtoms,  nHAtoms, Etot, E_pd, E_h, ftot, ftot_pd, ftot_h, rho_pd_avg, rho_h_avg, rho_pdh_avg, rho_hpd_avg, phitot, phitot_pdpd, phitot_pdh, phitot_hpd, phitot_hh

enersum = np.zeros((len(listptl),len(listr)+1,17))       # nH

for i in range(len(listptl)):
    for j in range(len(listr)):
        strener = '../htraj/'+listptl[i]+'r'+str(listr[j])+'_h1/ener.log'
        #enersum[i,j] = np.average(np.genfromtxt(strener,skip_header=1)[-100:],axis=0)[1:]
        enersum[i,j] = np.genfromtxt(strener,skip_header=1)[-1][1:]
    strener = '../htraj/'+listptl[i]+'h1/ener.log'
    enersum[i,-1] = np.genfromtxt(strener,skip_header=1)[-1][1:]

plt.clf()
strsymbol = ['-o','-s','-^']
for i in range(len(listptl)):
    plt.plot(listr, enersum[i,:11,2]/enersum[i,:11,0],strsymbol[i],color=listcolor[i],label=listlabel[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel('energy per PdH (eV)')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('energy_r.png',bbox_inches='tight')

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
# Hsite type vs particle radius
########################

 #istep, nPdAtoms, nHAtoms, nHType(1:6), nHPdGroupNN(1:3), nHPdGroupCoord(1:3)
htypesum = np.zeros((len(listptl),len(listr),14))       # nH

for i in range(len(listptl)):
    for j in range(len(listr)):
        strhtype = '../htraj/'+listptl[i]+'r'+str(listr[j])+'_h1/htype.dat'
        htypesum[i,j] = np.genfromtxt(strhtype,skip_header=1)[-1][1:]

plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, htypesum[i,:,2]/htypesum[i,:,0],'-o',color=listcolor[i],label=listlabel[i])
    plt.plot(listr, htypesum[i,:,3]/htypesum[i,:,0],'--^',color=listcolor[i])
    plt.plot(listr, (htypesum[i,:,2]+htypesum[i,:,4])/htypesum[i,:,0],'-.x',color=listcolor[i])

for i in range(len(listptl)):
    plt.plot(listr, htypesum[i,:,2]/htypesum[i,:,0],strsymbol[i],color=listcolor[1],label=listlabel[i])
    plt.plot(listr, htypesum[i,:,3]/htypesum[i,:,0],strsymbol[i],color=listcolor[2])

strsymbol = ['-o','-s','-^']
plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, (htypesum[i,:,2]+htypesum[i,:,4]+htypesum[i,:,6])/htypesum[i,:,0],strsymbol[i],color=listcolor[i],label=listlabel[i])
#    plt.plot(listr, (htypesum[i,:,2]+htypesum[i,:,6])/htypesum[i,:,0],'-o',color=listcolor[i],label=listlabel[i])
#    plt.plot(listr, (htypesum[i,:,2])/htypesum[i,:,0],'-.o',color=listcolor[i],label=listlabel[i])
#    plt.plot(listr, (htypesum[i,:,4])/htypesum[i,:,0],'-.x',color=listcolor[i])
#    plt.plot(listr, htypesum[i,:,3]/htypesum[i,:,0],'--^',color=listcolor[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel(r'highly coordinated H fraction')
plt.ylim(0,0.85)
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('prob_ocall.png',bbox_inches='tight')

plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, (htypesum[i,:,2]+htypesum[i,:,4])/htypesum[i,:,0],'-o',color=listcolor[i],label=listlabel[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel(r'$P_{OC}$')
plt.legend( prop=fontP,handletextpad=0.1)
plt.ylim(0,0.85)
plt.savefig('prob_ocall.png',bbox_inches='tight')


plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, htypesum[i,:,2]/htypesum[i,:,0],'-o',color=listcolor[i],label=listlabel[i])
    if i==0:
        plt.plot(listr, (htypesum[0,:,2]+htypesum[0,:,4]+htypesum[0,:,6])/htypesum[0,:,0],'-.x',color=listcolor[0],label=r'AET$^{*}$')

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel(r'$P_{OC}$')
plt.legend( prop=fontP,handletextpad=0.1)
plt.ylim(0,0.8)
plt.savefig('prob_oconly.png',bbox_inches='tight')

plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, htypesum[i,:,3]/htypesum[i,:,0],'-.^',color=listcolor[i],label=listlabel[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel(r'$P_{TE}$')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('prob_te.png',bbox_inches='tight')

plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, htypesum[i,:,5]/htypesum[i,:,0],'--x',color=listcolor[i],label=listlabel[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel(r'$P_{surf}$')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('prob_surf.png',bbox_inches='tight')

#####################
# Htype bar graph for r=6
#####################
# ==> anal_fcc_pdh/r6_h1/htype.dat <==
#       49991          55          55           7          20           2          26           0           0           0           0           0           0           0           0

#==> anal_hcp_pdh/r6_h1/htype.dat <==
#       49991          57          57          10           9          13          25           0           0           0           0           0           0           0           0

#==> anal_nanoptl/r6_h1/htype.dat <==
#       49991          57          57           5          14          11          26           1           0           0           0           0           0           0           0

#   surf, te, f=5, oc, crowded
#   1-3, 4, 5, 6, more
hlist = np.zeros((3,4))
hlist[0] = [26,14,11,6]
hlist[1] = [25,9,13,10]
hlist[2] = [26,20,2,7]

xlabels = ['0-3\nsurf','4','5','6-8']
xlabels = ['0-3\nsurface','4\ntetrahedral','5\n\nbipyramidal','6-8\noctahedral']

xlist = np.arange(4)

fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.bar(xlist + 0.00, data[0], color = 'k', width = 0.25)
ax.bar(xlist + 0.25, data[1], color = 'r', width = 0.25)
ax.bar(xlist + 0.50, data[2], color = 'b', width = 0.25)

listbar = [-.25,0,.25]
for i in range(len(listptl)):
    ax.bar(xlist+listbar[i], hlist[i],color=listcolor[i],width=0.25,label=listlabel[i])

plt.xlabel(r'number of coordinating Pd ($N_{Pd}^{coord}$)')
plt.xticks(xlist,xlabels)
plt.ylabel(r'number of H atoms ($N_H$)')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('htype_bar.png',bbox_inches='tight')


## htype 5-8 
hlist = np.zeros((3,3))
hlist[0] = [26,14,11+6]
hlist[1] = [25,9,13+10]
hlist[2] = [26,20,2+7]

xlabels = ['0-3\nsurface','4\ntetrahedral','5-8\nhighly coordinated']

xlist = np.arange(3)

fig = plt.figure()
ax = fig.add_axes([0,0,1,1])

listbar = [-.25,0,.25]
for i in range(len(listptl)):
    ax.bar(xlist+listbar[i], hlist[i],color=listcolor[i],width=0.25,label=listlabel[i])

plt.xlabel(r'number of coordinating Pd ($N_{Pd}^{coord}$)')
plt.xticks(xlist,xlabels)
plt.ylabel(r'number of H atoms ($N_H$)')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('htype_bar.png',bbox_inches='tight')


enersum = np.zeros((len(listptl),len(listr)+1,17))       # nH

for i in range(len(listptl)):
  for j in range(len(listr)):
    strener = '../htraj/'+listptl[i]+'r'+str(listr[j])+'_h1_surf/ener.log'
    enersum[i,j] = np.average(np.genfromtxt(strener,skip_header=1)[-100:],axis=0)[1:]


 #istep, nPdAtoms, nHAtoms, nHType(1:6), nHPdGroupNN(1:3), nHPdGroupCoord(1:3)
htypesum = np.zeros((len(listptl),len(listr),14))       # nH

for i in range(len(listptl)):
    for j in range(len(listr)):
        strhtype = '../htraj/'+listptl[i]+'r'+str(listr[j])+'_h1_surf/htype.dat'
        #enersum[i,j] = np.average(np.genfromtxt(strener,skip_header=1)[-100:],axis=0)[1:]
        htypesum[i,j] = np.genfromtxt(strhtype,skip_header=1)[-1][1:15]


plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, (htypesum[i,:,2]+htypesum[i,:,4]+htypesum[i,:,6])/htypesum[i,:,0],'-o',color=listcolor[i],label=listlabel[i])
    plt.plot(listr, (htypesum[i,:,2])/htypesum[i,:,0],'--x',color=listcolor[i],label=listlabel[i])
    plt.plot(listr, (htypesum[i,:,4])/htypesum[i,:,0],'-.^',color=listcolor[i],label=listlabel[i])
    plt.plot(listr, (htypesum[i,:,5])/htypesum[i,:,0],'-x',color=listcolor[i],label=listlabel[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel(r'$P_{OC}$')
plt.legend( prop=fontP,handletextpad=0.1)
plt.ylim(0,0.85)
plt.savefig('prob_ocall.png',bbox_inches='tight')


plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, htypesum[i,:,2]/htypesum[i,:,0],'-o',color=listcolor[i],label=listlabel[i])
    if i==0:
        plt.plot(listr, (htypesum[0,:,2]+htypesum[0,:,4]+htypesum[0,:,6])/htypesum[0,:,0],'-.x',color=listcolor[0],label=r'AET$^{*}$')

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel(r'$P_{OC}$')
plt.legend( prop=fontP,handletextpad=0.1)
plt.ylim(0,0.8)
plt.savefig('prob_oconly.png',bbox_inches='tight')

plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, htypesum[i,:,3]/htypesum[i,:,0],'-.^',color=listcolor[i],label=listlabel[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel(r'$P_{TE}$')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('prob_te.png',bbox_inches='tight')

plt.clf()
for i in range(len(listptl)):
    plt.plot(listr, htypesum[i,:,3]/htypesum[i,:,0],'--x',color=listcolor[i],label=listlabel[i])

plt.xlabel(r'particle radius ($\AA$)')
plt.ylabel(r'$P_{surf}$')
plt.legend( prop=fontP,handletextpad=0.1)
plt.savefig('prob_surf.png',bbox_inches='tight')


for i in range(len(listptl)):
    plt.plot(listr, enersum[i,5:,2]/enersum[i,:,0],'-o',color=listcolor[i],label=listlabel[i])


plt.legend()

    enersum[i,j,0]=enerdata[1]      #nHAtoms
    enersum[i,j,1]=enerdata[7]      #totEner

listptl = ['fcc_pd4','fcc_pdh','hcp_pd2','hcp_pdhoct','nanoptl']
listsim = ['mun3','mup0','mup3','r6_mup3','r10_mup3','r15_mup3','r20_mup3']
listsim = ['mup3','r6_mup3','r10_mup3','r15_mup3','r20_mup3']
listrval = [30,6,10,15,20]

npdlist = np.asarray([[12949,55,321,959,2243,7803],
[11022,55,225,791,1961,6531],
[12983,57,293,955,2283,7699],
[11137,57,257,829,1965,6691],
[9954,59,265,863,2031,6868],])


enersum = np.zeros((len(listptl),len(listsim),5))       # nH

#istep,   nHAtoms,     nAccpt,      nReject,      rAccpt,  ftot, phitot, ftot+phitot, ftot_pd, ftot_h, phitot_pdh, phitot_hh

for i in range(len(listptl)):
  for j in range(len(listsim)):
    strener = '../'+listptl[i]+'/'+listsim[j]+'/ener.log'
    enerdata = np.genfromtxt(strener,skip_header=1)[-1]
    npd = npdlist[i,j]
    enersum[i,j,0]=enerdata[1]      #nHAtoms
    enersum[i,j,1]=npd              #nPdAtoms
    enersum[i,j,2]=enerdata[7]      #totEner
    enersum[i,j,3]=enerdata[7]/enerdata[1]      #totEner/nHAtoms
    enersum[i,j,4]=enerdata[7]/npd      #totEner/nPd



xlist=[0,1,2,3,4]
plt.scatter(xlist,enersum[xlist,0,3])
xlist=[1,3,4]
plt.scatter(xlist,enersum[xlist,0,3])

for i in range(len(listptl)):
    xlist

listptl = ['fcc_pdh','hcp_pdhoct','nanoptl']
listsim = ['r6_h1','r10_h1','r15_h1','r20_h1','r30_h1','h1']
listr = [6,10,15,20,30,39]

enersum = np.zeros((len(listptl),len(listsim),3))       # nH

for i in range(len(listptl)):
  for j in range(len(listsim)):
    strener = '../htraj/'+listptl[i]+'/'+listsim[j]+'/ener.log'
    enerdata = np.genfromtxt(strener,skip_header=1)[-1]
    enersum[i,j,0]=enerdata[1]      #nHAtoms
    enersum[i,j,1]=enerdata[7]      #totEner
    enersum[i,j,2]=enerdata[7]/enerdata[1]      #totEner/nHAtoms


xlist = [0,1,2,3,4,5]

