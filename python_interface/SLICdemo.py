import os
import subprocess
import ctypes
from ctypes import cdll
from ctypes import POINTER,c_int,c_double,c_bool
from PIL import Image
# from skimage.io import imread,imshow
import numpy as np
from timeit import default_timer as timer

# Create image with random colors
def createimage(labels,numlabels):
	lo,hi = 80,256
	r = np.random.randint(lo,hi,(numlabels))
	g = np.random.randint(lo,hi,(numlabels))
	b = np.random.randint(lo,hi,(numlabels))
	limg =  np.array([r[labels],g[labels],b[labels]],dtype=np.uint8).transpose(1,2,0)
	return limg

# Draw boundaries
def drawBoundaries(imgname,labels,numlabels):

	img = Image.open(imgname)
	img = np.array(img)

	ht,wd = labels.shape

	for y in range(1,ht-1):
		for x in range(1,wd-1):
			if labels[y,x-1] != labels[y,x+1] or labels[y-1,x] != labels[y+1,x]:
				img[y,x,:] = 0

	return img

def segment(imgname,numsuperpixels,compactness,doRGBtoLAB):
	#--------------------------------------------------------------
	# read image and change image shape from (h,w,c) to (c,h,w)
	#--------------------------------------------------------------
	img = Image.open(imgname)
	img = np.asarray(img)

	dims = img.shape
	h,w,c = dims[0],dims[1],1
	if len(dims) > 1:
		c = dims[2]
		img = img.transpose(2,0,1)
		print(c, "channels")
	#--------------------------------------------------------------
	# Reshape image to a single dimensional vector
	#--------------------------------------------------------------
	img = img.reshape(-1)
	#--------------------------------------------------------------
	# Prepare pointers to pass to the C function
	#--------------------------------------------------------------
	inp = img.astype(ctypes.c_double)
	labels = np.zeros((h,w), dtype = ctypes.c_int) # does not matter if the shape is (1,h*w) or (h,w)
	numlabels = np.zeros(1,dtype = ctypes.c_int)

	pinp 		= inp.ctypes.data_as(POINTER(c_double))
	plabels 	= labels.ctypes.data_as(POINTER(c_int))
	pnumlabels 	= numlabels.ctypes.data_as(POINTER(c_int))
	ww = ctypes.c_int(w)
	hh = ctypes.c_int(h)
	cc = ctypes.c_int(c)
	numsp = ctypes.c_int(numsuperpixels)
	comp = ctypes.c_double(compactness)
	colcon = ctypes.c_bool(doRGBtoLAB)

	#--------------------------------------------------------------
	# Load library and call
	#--------------------------------------------------------------
	lib = cdll.LoadLibrary('libslic.so')
	start = timer()
	lib.SLICmain(pinp,ww,hh,cc,numsp,comp,colcon,plabels,pnumlabels)
	end = timer()
	#--------------------------------------------------------------
	# Collect labels
	#--------------------------------------------------------------
	outlabels = np.array(np.fromiter(plabels, dtype=np.int32, count=labels.size))
	# print(outlabels.size)
	print("number of output superpixels: ", numlabels[0])
	print("segmentation time taken in seconds: ", end-start)
	#--------------------------------------------------------------
	# Display labels
	#--------------------------------------------------------------
	# labelimg = Image.fromarray(outlabels.reshape((h,w)))
	return outlabels.reshape(h,w),numlabels[0]



def slicdemo():
	#--------------------------------------------------------------
	# Create shared library
	#--------------------------------------------------------------
	if os.path.exists("libslic.so"):
		print("Compiled library exists")
	else:
		subprocess.call("gcc -c -fPIC slicpython.c -o slic.o", shell=True)
		subprocess.call("gcc -shared slic.o -o libslic.so", shell = True)
		print("library compiled")

	numsuperpixels = 500
	compactness = 20.0
	doRGBtoLAB = True # only works if it is a three channel image
	# imgname = "/Users/achanta/Pictures/classics/lena.png"
	imgname = "bee.png"
	labels,numlabels = segment(imgname,numsuperpixels,compactness,doRGBtoLAB)

	# labelimg = createimage(labels,numlabels)
	segimg = drawBoundaries(imgname,labels,numlabels)
	Image.fromarray(segimg).show()

	return


slicdemo()



