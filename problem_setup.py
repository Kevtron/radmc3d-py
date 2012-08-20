import numpy as np
import math
import progressbar as pbar
import operator
progress = pbar.ProgressBar()

def readDensity(densfile):
	dfile=open(densfile,'r')
	test=dfile.readlines()
	nx = int(test[18].split()[0])
	ny = int(test[18].split()[1])
	nz = int(test[18].split()[2])
	temp = np.zeros((nx,ny*nz))
	rhod = np.zeros((nx,ny,nz))
	for j in progress(range(ny*nz)):
		for i in range(nx):
			temp[i,j]=float(test[19+j].split()[i])
			if temp[i,j] < 0.25:
				temp[i,j]=0.0
	rhod=temp.reshape(nx,ny,nz)
	dfile.close()
	return rhod,nx,ny,nz

def writeDensities(density,nx,ny,nz)
	dust=open('dust_density.inp','w')
	stars=open('stellarsrc_density.inp','w')
	dust.write('1\n')
	dust.write('%d\n' % nx*ny*nz)
	dust.write('1\n')
	stars.write('1\n')
	stars.write('%d\n' % nx*ny*nz)
	stars.write('1\n')
	for iz in range(nz):
		for iy in range(ny):
			for ix in range(nx):
				dust.write('%f\n' % density[ix,iy,iz])
				stars.write('%f\n' % density[ix,iy,iz])
	dust.close()
	stars.close()
	return

def makeAMRgrid(nx,ny,nz,sizex,sizey,sizez):
	xi = [-sizex + 2*sizex*(i)/(nx) for i in range(nx+1)] 
	yi = [-sizey + 2*sizey*(i)/(ny) for i in range(ny+1)]
	zi = [-sizey + 2*sizez*(i)/(nz) for i in range(nz+1)]
	amr=open('amr_grid.inp','w')
	amr.write('1\n') 					#iformat
	amr.write('0\n')					#AMR grid style
	amr.write('0\n')					#coordinate system
	amr.write('0\n')					#grid info
	amr.write('1,1,1\n')				#include x,y,z
	amr.write('%i,%i,%i\n' % nx,ny,nz)				#size of grid
	[amr.write('%f\n' % xi[i]) for i in range(nx+1)]
	[amr.write('%f\n' % yi[i]) for i in range(ny+1)]
	[amr.write('%f\n' % zi[i]) for i in range(nz+1)]
	amr.close()

def writeControlFiles(nphot):
	radmc3d=open('radmc3d.inp','w')
	dustopac=open('dustopac.inp','w')
	dustopac.write("""2               Format number of this file
1               Nr of dust species
1               Way in which this dust species is read
0               0=Thermal grain
starstuff       Extension of name of dustkappa_***.inp file""")
	dustopac.close()
	radmc3d.write('nphot = %f\n' % nphot)
	radmc3d.write('scattering_mode_max = 1\n')
	radmc3d.write('iranfeqmode = 1\n')
	radmc3d.write('nphot_spec = 100000\n')
	radmc3d.write('incl_userdef_srcalp = 1\n')
	radmc3d.write('incl_freefree = 1\n')
	radmc3d.write('incl_dust = 1\n')

def readTemperature(tempfile,MH,KB):
	dfile=open(densfile,'r')
	test=dfile.readlines()
	nx = int(test[18].split()[0])
	ny = int(test[18].split()[1])
	nz = int(test[18].split()[2])
	array = np.zeros((nx,ny*nz))
	temp = np.zeros((nx,ny,nz))
	for j in progress(range(ny*nz)):
		for i in range(nx):
			array[i,j]=float(test[19+j].split()[i])
			if array[i,j] < 0.25:
				array[i,j]=0.0
	temp=array.reshape(nx,ny,nz)*(2.0*MH)/(3.0*KB)
	dfile.close()
	return temp,nx,ny,nz

def writeTemperature(temperature,nx,ny,nz)
	temp=open('dust_temperature.dat')
	temp.write('1\n')
	temp.write('%f\n' % nx*ny*nz)
	temp.write('1\n')
	for iz in range(nz):
		for iy in range(ny):
			for ix in range(nx):
				temp.write('%f\n' % temperature[ix,iy,iz])
	temp.close()
	return
	
def makeWavelength():	
	lambda1=0.1
	lambda2=7.0
	lambda3=25.0
	lambda4=1e4
	n12=20
	n23=100
	n34=30
	lam12 = [lambda1*(lambda2/lambda1)**(i/(1.0*n12)) for i in range(n12)]
	lam23 = [lambda2*(lambda3/lambda2)**(i/(1.0*n23)) for i in range(n23)]
	lam34 = [lambda3*(lambda4/lambda3)**(i/(1.0*n34)) for i in range(n34)]
	[lam12.append(lam23[i]) for i in range(n23)]
	[lam12.append(lam34[i]) for i in range(n34)]
	#
	# Print wavelength.pinp
	#
	f=open('wavelength_micron.inp','w')
	f.write('%d\n' % len(lam12))
	for i in range(len(lam12)):
	f.write('%f\n' % lam12[i])
	f.close()
	return lam12,len(lam12)

def stellarTemplate(lam,nlam,nx,ny,nz,tstar,RS,MS):
	template=open('stellarsrc_template.inp')
	template.write('2\n')
	template.write('%d\n' % (nx*ny*nz))
	template.write('%d\n' % nlam)
	[template.write('%f\n'%lam[i]) for i in range(nlam)]
	template.write('-%d\n' % tstar)
	template.write('%f' % RS)
	template.write('%f' % MS)
	template.close()
	return
	
def dustOpacity(filename)
	ifile = open(filename,'r')
	ofile = open('dustkappa_starstuff.inp','w')	
	itemp=ifile.readlines()
	wav=reduce(operator.add, [[float(itemp[2+j].split()[i]) for i in range(9)] for j in range(106)]) 	
	wav.append(float(itemp[109]))
	totalopacity=float(itemp[110].split()[7])
	opacity=reduce(operator.add, [[float(itemp[111+j].split()[i]) for i in range(6)] for j in range(357)])
	[ofile.write('%f %f %f\n' % wav[i], opacity[2*i]*totalopacity, opacity[2*i+1]*totalopacity) for i in range(len(wav))]
	ifile.close()
	ofile.close()
#Physical Constants
#
AU = 1.49598e13 # Astronomical Unit 	[cm]
pc = 3.08572e18 # Parsec		[cm]
MS = 1.98892e33 # Solar Mass		[g]
TS = 5.78e3     # Solar Temp		[K]
LS = 3.8525e33	# Solar Luminosity	[erg/s]
RS = 6.96e10	# Solar radius		[cm]
#
# Monte Carlo parameters
#
nphot = 1000000
#
# Read the dust density stuff
#
densfile = open('/Users/ksooley/Documents/Research/relaxsingle/grids/MASS008_001_density_g_cm3_grid.dat','r')
test=densfile.readlines()
nx = int(test[18].split()[0])
ny = int(test[18].split()[1])
nz = int(test[18].split()[2])
temp = np.zeros((nx,ny*nz))
rhod = np.zeros((nx,ny,nz))
for j in progress(range(ny*nz)):
	for i in range(nx):
		temp[i,j]=float(test[19+j].split()[i])
rhod=temp.reshape(nx,ny,nz)
#
# Model Parameters
#

#
# Star parameters
# 
mstar = MS
rstar = RS
tstar = TS
pstar = (0.0,0.0,0.0)
sizex = 1.17*RS
sizey = 1.17*RS
sizez = 1.09*RS
#
# Make Coordinates
#
xi = [-sizex + 2*sizex*(i)/(nx) for i in range(nx+1)] 
yi = [-sizey + 2*sizey*(i)/(ny) for i in range(ny+1)]
zi = [-sizey + 2*sizez*(i)/(nz) for i in range(nz+1)]
#
# Make Density Model
#rhod = [[[rho0 * math.exp(-((rr[i][j][k])**2/radius**2)/2) for i in range(nx)] for j in range(ny)] for k in range(nz)] 
rhod=[[[0 for i in range(nx)] for j in range(ny)] for k in range(nz)]
for k in range(nz):
	for j in range(ny):
		for i in range(nx):
			rhod[i][j][k]=rho0*math.exp(-((rr[i][j][k]**2)/(radius**2))/2)
#
# Wavelength information
#
lambda1=0.1
lambda2=7.0
lambda3=25.0
lambda4=1e4
n12=20
n23=100
n34=30
lam12 = [lambda1*(lambda2/lambda1)**(i/(1.0*n12)) for i in range(n12)]
lam23 = [lambda2*(lambda3/lambda2)**(i/(1.0*n23)) for i in range(n23)]
lam34 = [lambda3*(lambda4/lambda3)**(i/(1.0*n34)) for i in range(n34)]
[lam12.append(lam23[i]) for i in range(n23)]
[lam12.append(lam34[i]) for i in range(n34)]
#
# Print wavelength.pinp
#
f=open('wavelength_micron.inp','w')
f.write('\t%d\n' % len(lam12))
for i in range(len(lam12)):
	f.write('\t%3.9f\n' % lam12[i])
f.close()
f.close
#
# Write the stars file
#
f=open('stars.inp','w')
f.write('\t%d\n' % 2)
f.write('\t%d\t%d\n' % (1, len(lam12)))
f.write('\n\t %e \t %e \t %3.6f \t %3.6f \t %3.6f \n' % (rstar, mstar, pstar[0], pstar[1], pstar[2]))
f.write('\n')
for i in range(len(lam12)):
	f.write('\t%3.9f\n' % lam12[i])
f.write('\n')
f.write('\t -%f' %tstar )
f.close()
#
# Write the grid file
#
f=open('amr_grid.inp','w')
f.write('\t 1 \n') #iformat
f.write('\t 0 \n') #AMR grid style
f.write('\t 0 \n') #coordinate system
f.write('\t 0 \n') #grid info
f.write('\t 1 \t 1 \t 1 \n') # x y z coords
f.write('\t %d \t %d \t %d \n' % (nx,ny,nz)) # size of grid
for i in range(nx+1):
	f.write('\t %e \n' % xi[i])
for i in range(ny+1):
	f.write('\t %e \n' % yi[i])
for i in range(nz+1):
	f.write('\t %e \n' % zi[i])
f.close()
#
# Write the density file
#
f=open('dust_density.inp','w')
f.write('\t 1 \n') #format number
f.write('\t %d \n' % (nx*ny*nz)) # number of cells
f.write('\t 1 \n') # number of dust species
for k in range(nz):
	for j in range(ny):
		for i in range(nx):
			f.write('\t %e \n' % rhod[i][j][k])
#
# Dust opacity control file
#
f=open('dustopac.inp','w')
f.write("""2               Format number of this file
1               Nr of dust species
============================================================================
1               Way in which this dust species is read
0               0=Thermal grain
silicate        Extension of name of dustkappa_***.inp file
----------------------------------------------------------------------------""")
f.close
#
# radcmc3d.inp control file
#
f=open('radmc3d.inp','w')
f.write('nphot = \t %d \n' % nphot)
f.write('scattering_mode_max = \t  0\n')
f.write('iranfreqmode = 1')
f.close
