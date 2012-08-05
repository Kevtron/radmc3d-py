#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#             PYTHON PACKAGE TO READ AND PLOT RESULTS FROM RADMC-3D
#		     Based on the included readradmc.pro
#			Written by: Kevin Sooley
#			      July 22, 2012
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#           		     ROUTINES FOR IMAGES
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
#                     READ THE RECTANGULAR TELESCOPE IMAGE
#---------------------------------------------------------------------------

def readimage(filename=None, imagefile=None ):
	if filename is None:
		filename='image.out'
	f=open(filename)
	lines=f.readlines()
	iformat=int(lines[0])
	nx,ny=lines[1].split()
	nf=int(lines[2])
	sizepix_x,sizepix_y=lines[3].split()
	wavelength=lines[4]
	temp=np.zeros(100*100)
	for i in range(len(lines[6:-1])): 
		temp[i]=lines[6+i]
	image=temp.reshape(-1,100).transpose()
	xi=(np.arange(int(nx))-int(nx)/2+0.5)*float(sizepix_x)
	yi=(np.arange(int(ny))-int(ny)/2+0.5)*float(sizepix_y)
	# do a contour plot of the image
	fig = plt.figure()
	ax = fig.add_subplot(111,aspect='equal')
	cax = ax.imshow(image,cmap=plt.cm.bone, interpolation='nearest',extent=[-float(sizepix_x)*int(nx)/2,float(sizepix_x)*int(nx)/2,-float(sizepix_y)*int(ny)/2,float(sizepix_y)*int(ny)/2])
	ax.contour(image,cmap=plt.cm.bone)
	# compute the flux
	flux=np.sum(image)
	pc=3.0857200e+18
	flux=flux/pc**2
	cbar = fig.colorbar(cax)
	plt.title(r'flux at 1 pc = %s erg/cm^2/s/Hz' % flux)
	if imagefile is None:
		imagefile="image.png"
	fig.savefig(imagefile)
	#plt.show()
	return 

#---------------------------------------------------------------------------
#               read and plot the spectrum
#---------------------------------------------------------------------------
def B_lambda(T,wavelength):
	h=6.62606957e-27
	c=29979245800
	kb=1.3806488e-16
	return (2.0*h*c**2/wavelength**5)*(1.0/(np.exp(h*c/(wavelength*kb*T))-1.0)) 

def readspectrum(filename = None, imagefile = None):
	if filename is None:
		filename = 'spectrum.out'
	f = open(filename)
	lines = f.readlines()
	nlam = int(lines[1])
	wavelength = np.zeros(nlam)
	flux = np.zeros(nlam)
	for i in range(nlam):
		wavelength[i],flux[i] = lines[i+3].split()
	wave=np.logspace(-1,12,num=300)
	plt.loglog(wavelength,flux,'r-')
	#plt.loglog(wave,B_lambda(4000,wave*10**-4),'b-')
	plt.show()

#---------------------------------------------------------------------------
# 	Let's get fancy and call radmc3d from inside python
#---------------------------------------------------------------------------
#def makeimage():

if __name__ == "__main__":
	#readimage('image_U.out', 'image_U.png')
	#readimage('image_B.out', 'image_B.png')
	#readimage('image_V.out', 'image_V.png')
	#readimage('image_R.out', 'image_R.png')
	#readimage('image_I.out', 'image_I.png')
	readspectrum()
