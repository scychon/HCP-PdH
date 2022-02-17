from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits import mplot3d
from matplotlib.ticker import FormatStrFormatter
from matplotlib.font_manager import FontProperties
import os

def writexyz(pdpos,outfile):
    fout = open(outfile, "w+")
    npd = len(pdpos)
    fout.write(str(npd)+'\n')
    fout.write('generated from rescale_cell.py\n')
    for i in range(len(pdpos)):
        xpos = pdpos[i]
        fout.write(' Pd%16.6f%16.6f%16.6f\n' % (xpos[0],xpos[1],xpos[2]))
    return


listr = [6,8,10,12,15,17,20,22,25,27,30]
listfrac = ['.2', '.25','.3','.4', '.5','.6', '.75', '.8']


##########
# for hcp
##########
for i in range(len(listr)):
    strfile = 'hcp_pdh_r'+str(listr[i])+'.xyz'
    for j in range(len(listfrac)):
        pdxyz = np.genfromtxt(strfile,skip_header=2)[:,1:]
        strout = 'hcp_pdh'+listfrac[j]+'_r'+str(listr[i])+'.xyz'
        ascale = ((2.871-2.728)*float(listfrac[j])+2.728)/2.871
        cscale = ((4.7846-4.563)*float(listfrac[j])+4.563)/4.7846
        pdxyz[:,0:2] = pdxyz[:,0:2]*ascale
        pdxyz[:,2] = pdxyz[:,2]*cscale
        writexyz(pdxyz,strout)

##########
# for fcc
##########
for i in range(len(listr)):
    strfile = 'fcc_pdh_r'+str(listr[i])+'.xyz'
    for j in range(len(listfrac)):
        pdxyz = np.genfromtxt(strfile,skip_header=2)[:,1:]
        strout = 'fcc_pdh'+listfrac[j]+'_r'+str(listr[i])+'.xyz'
        ascale = ((4.104-3.8905)*float(listfrac[j])+3.8905)/4.104
        listfrac[j], ascale
        pdxyz_frac = pdxyz[:,0:3]*ascale
        writexyz(pdxyz_frac,strout)


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




def gen_fcc(nx,ny,nz):
    matfcc = np.zeros((nx*2+1,ny*2+1,nz*2+1),dtype=int)
    for i in range(nx*2+1):
        for j in range(ny*2+1):
            for k in range(nz*2+1):
                if (i+j+k)%2==0:
                    matfcc[i,j,k]=1
    return matfcc

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





# Pd fcc
aval = 3.887
fname = 'fcc_pd4.xyz'

# PdH fcc
aval = 4.104
fname = 'fcc_pdh.xyz'

# get min & max r positions
pdxyz = np.genfromtxt('PdH_xyz_AddDomType.xyz',skip_header=2,skip_footer=1)[:,1:]
pdpos = pdxyz[:,:3]
rmin = np.min(pdpos,axis=0)-3
rmax = np.max(pdpos,axis=0)+3
npd = len(pdpos)

nz=int(np.ceil((rmax-rmin)[2]/aval))
ny=int(np.ceil((rmax-rmin)[1]/aval))
nx=int(np.ceil((rmax-rmin)[0]/aval))


xoutlist = []
for i in range(nx*2+1):
  for j in range(ny*2+1):
    for k in range(nz*2+1):
      if (i+j+k)%2==0:
        xval =  rmin+ aval/2*np.asarray((i,j,k))
        if (xval[2]<=rmax[2] and xval[0] <= rmax[0] and xval[1]<=rmax[1]) :
            xoutlist.append(' Pd%10.3f000%10.3f000%10.3f000\n' % (xval[0],xval[1],xval[2]))


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


vec0=np.array((1,0,0))
vec1 = np.array((np.cos(np.pi*2/3),np.sin(np.pi*2/3),0))
vec2 = np.array((0,0,1))

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
