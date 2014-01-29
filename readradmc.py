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
from optparse import OptionParser
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
	temp=np.zeros(int(nx)*int(ny))
	for i in range(len(lines[6:-1])): 
		temp[i]=lines[6+i]
	image=temp.reshape(-1,int(nx)).transpose()
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
	c=299792458
	kb=1.3806488e-16
	return (2.0*h*c**2/wavelength**5)*(1.0/(np.exp(h*c/(wavelength*kb*T))-1.0)) 

def readspectrum(filename = None, imagefile = None):
	h=6.62606957e-27
	c=299792458
	kb=1.3806488e-16
	if filename is None:
		filename = 'spectrum.out'
	if imagefile is None:
		imagefile = 'spectrum.png'
	f = open(filename)
	lines = f.readlines()
	nlam = int(lines[1])
	wavelength = np.zeros(nlam)
	flux = np.zeros(nlam)
	for i in range(nlam):
		wavelength[i],flux[i] = lines[i+3].split()
	wave=np.arange(1e-6,5,10000)
	#plt.xlim(1e-6,1e-5)
	#plt.semilogy(wavelength/1e-6,flux/flux[0],'r-')
	plt.plot(wave,(2.0*h*c**2/wave**5)*(1.0/(np.exp(h*c/(wave*kb*4000))-1.0)),'r-')
	plt.xlabel(r'$\lambda (m)$')
	plt.ylabel(r'Flux')
	plt.savefig(imagefile)

if __name__ == "__main__":
	usage = """ %prog [options] FILE """
	parser = OptionParser(usage=usage, version = "%prog")
	parser.add_option("-i", "--image", action = "store_true", dest = "im", default = False, help = "Read image")
	parser.add_option("-s", "--spectrum", action = "store_true", dest = "spec", default = False, help = "read spectrum")
	(options, args) = parser.parse_args()
        if options.im == True:
		readimage()
	if options.spec == True:
		readspectrum()
