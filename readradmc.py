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
	plt.show()
if __name__ == "__main__":
	readimage()
