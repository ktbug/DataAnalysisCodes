#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 14:14:13 2018

#Translated by:Katie Kosak
### Original Authors: Christopher Goodard in IDL
### Original Description of Code from IDL
##
#
# :Description:
#    Calculates a curved slit between minimum three points.
#    The points are connected with a spline.
#    In the case of two points, a linear slit is given 
#  Syntax:
#  IDL>curved_time_distance, data, prof, points=points, w=w
#
# :Params:
#    INPUT:
#    data   --> image array
#    points --> FLTARR[2,*] array which given x (points[0,*]) and y (points[1,*])coordinates
#               of points which define the curved strip. If exists, the coordinates stored in the 
#               variable points will be load. If not, the cursor procedure will start.
#       w   --> width of the slit
#
#    OUTPUT:
#    prof   --> time-distance map array
#  
#   CALLS:
#   make_norm.pro
#   HYSTORY:
#
#  Improvementes contribution from Sergey Anfinogentov
#
#  30 May, 2013: /seq keyword for making a sequence of image from the path
#  03 Jun, 2013: add of the keyword array
#  16 Feb, 2014: add the keyword "poly"
"""
import numpy as np
from scipy.interpolate import interp1d,griddata,interp2d
import math 
import scipy
import numpy as np
import scipy.io
import os
import numpy as np
import matplotlib.pyplot as plt
import os.path
import subprocess
from tempfile import mkdtemp
import matplotlib.animation as animation
from astropy.convolution import convolve, Box1DKernel
import itertools
name='ccube-kumar2013-aia131.sav'
Movie=scipy.io.readsav(name,python_dict=False)
data=Movie['ccube']
frameCount=data.shape[0]
width=5
if (width%2==0):
    width+1
name1='data.sav'
Comparison=scipy.io.readsav(name1,python_dict=False)


#plt.figure()
#plt.imshow(data[100,:,:],vmin=0,vmax=30)
#Points=plt.ginput(7)
#plt.show()
#Points=np.array(Points,dtype=float)
#
#
Points=np.array([[359.3815323 , 184.3942354 ],
       [353.7582283 , 193.65614787],
       [347.14257654, 216.48014644],
       [319.68762173, 246.58136195],
       [265.4392773 , 285.94448993],
       [237.98432249, 285.61370734],
       [225.74536673, 250.55075301]])

x=Points[:,0]
y=Points[:,1]


#### Now Fit a Spline
xnew=np.arange(x.min(),x.max(),0.5)
f=interp1d( x,y)
ynew=f(xnew)
xr=xnew
yr=ynew

def make_norm(r):
    rx=r[1]
    ry=r[0]
    if rx==0:
        values=[1,0]
    elif ry==0:
        values=[0,1]
    else:
        nx=np.sqrt((ry**2)/(rx**2+ry**2))*(ry/abs(ry))
        ny=-rx*(nx/ry)
        values=np.array([ny,nx])
    return values

######## Calculate the Grid ---- Sergey

n=xr.size # number of points along the slit
time_distance=np.zeros((n,frameCount),dtype=float)
sx=np.zeros((n,width),dtype=float)
sy=np.zeros((n,width),dtype=float)
y=(np.arange(width)-width/2) 
y=np.array([y,]*n)  #replicate(1.,n); x,y- coordinates in slit coordinate system
x=np.arange(n)
x=np.array([x,]*width).transpose()
###replicate(1.,width)*0.0 
for i in range(n):
    g1=np.array([yr[i]-yr[i-1],xr[i]-xr[i-1]])

    g1=g1/math.sqrt((g1**2).sum())          # g1- a basis vector along the slit direction
    g2=make_norm(g1) ### Create the g2 basis vector orthogonal to slit direction      
   
    sy[i,:]=g1[0]*x[i,:]+g2[0]*y[i,:]+xr[0] #sx,sy - coordinate in image CS 
    sx[i,:]=g1[1]*x[i,:]+g2[1]*y[i,:]+yr[0]

time_distance=np.zeros((n,frameCount),dtype=float)
#### Let's Iterate over the n points and over the width of the slit 
####### If the width==1, essentially it is just using the values of xnew,ynew from intiial slit

#for i in range(frameCount):
#    for k in range(width):
#        temp=np.empty(width)
#        for j in range(n):
#            upper_x=np.ceil(sx[j,k])+1
#            lower_x=np.floor(sx[j,k])
#            upper_y=np.ceil(sy[j,k])+1
#            lower_y=np.floor(sy[j,k])
#            data_xind=np.array([lower_x,upper_x],dtype=int)
#            data_yind=np.array([lower_y,upper_y],dtype=int)
#            ####### To save on computation time 
#            ############# Create a grid so that you have just the area surrounding the pixels of interest
#            ############## We use linear interpolation--- 4 closest pixels
#            xx,yy=np.meshgrid(data_yind,data_xind)
#            data_temp=data[i,data_yind[0]:data_yind[1],data_xind[0]:data_xind[1]]
#            grid_t=interp2d(xx,yy,data_temp,kind='linear')
#            tmp=grid_t(xr[j],yr[j])
#            temp[k]=tmp
#        value=np.sum(temp)/width
#        time_distance[j,i]=value
box1=0
for i in range(frameCount):
    data2=data[i,:,:]
    for j in range(n):
        upper_x=np.ceil(xr[j])+box1
        lower_x=np.floor(xr[j])-box1
        upper_y=np.ceil(yr[j])+box1
        lower_y=np.floor(yr[j])-box1
        data_xind=np.arange(lower_x,upper_x+1,dtype=np.intp)
        data_yind=np.arange(lower_y,upper_y+1,dtype=np.intp)
        ####### To save on computation time 
        ############# Create a grid so that you have just the area surrounding the pixels of interest
        ############## We use linear interpolation--- 4 closest pixels
        yy,xx=np.meshgrid(data_yind,data_xind)
        data_temp=data2[np.ix_(data_yind,data_xind)]
        grid_t=interp2d(yy,xx,data_temp,kind='linear')
        #grid_t=scipy.interpolate.RectBivariateSpline(data_yind,data_xind,data_temp)#,kind='cubic')
        tmp=grid_t(xr[j],yr[j])

        time_distance[j,i]=tmp

plt.figure()
plt.imshow(time_distance,vmin=0,vmax=30)
plt.show()


#### Let's Do the Fourier Transform on Different Slits
#### Defining three frames
frame1=0
frame2=100
frame3=40

slice1=time_distance[:,frame1]
slice2=time_distance[:,frame2]
slice3=time_distance[:,frame3]
slice3=np.pad(slice3,244,'constant',constant_values=0)
smoothed_signal = convolve(slice3, Box1DKernel(20))
frequency= np.fft.fftfreq(len(smoothed_signal)) 

fourier3=np.fft.fft(smoothed_signal)
ps3 = np.abs(fourier3)**2


plt.figure()
plt.plot(smoothed_signal)
plt.show()

plt.figure()
plt.plot(frequency,ps3)
plt.xlim(0.08,0.5)
plt.ylim(0,10000)
plt.show()
############ Let's Animate the Power Spectrum over the Frames

####### Since we are using ffmpeg, we will plot the figuers in a subdirectory
#########
subfolder='Data'
try:
    os.mkdir(subfolder)
except Exception:
    pass
# Change to the directory
os.chdir(subfolder)

for i in range(frameCount):
    image_slice=time_distance[:,i]
    image_slice=np.pad(image_slice,244,'constant',constant_values=0)
    image_slice= convolve(image_slice, Box1DKernel(20))
    fourier=np.fft.fft(image_slice) # Cycles per unit
    freq=np.fft.fftfreq(n)
    fl=int(np.floor(n/2))
    fourier=fourier[0:fl]
    freq=freq[0:fl]
    ps1=np.abs(fourier)**2
    idx=np.argsort(freq)
    ps1=ps1[idx]
    
    plt.figure()
    ## Plot Original Signal
    plt.subplot(221)
    plt.plot(image_slice)
    plt.title('Original Signal')
    
    
    
    ### Plot Power Spectra 
    plt.subplot(223)
    plt.plot(freq,ps1)
    plt.ylim(0,50000)
    plt.xlim(0.08,np.max(freq))
    plt.title('Power Spectra')
    plt.xlabel('Frequency')
    ### Show the Slit of the Image
    plt.subplot(122)
    plt.imshow(time_distance.T,vmin=0,vmax=100)
    plt.xlabel('FrameCount')
    plt.plot([0,time_distance.shape[0]-1],[i,i], 'm', lw=3)
    #plt.tight_layout()
    plt.savefig(str(i)+'.png')
    
os.chdir('..')
###### Animage the png files in the folder ##########
 #Movie Name
result_movie='Spectrum.mp4'
title='Spectrum'
####### Generate Movie
filepath1=os.path.join( os.path.dirname(os.path.realpath("__file__")),subfolder)
params='ffmpeg -framerate 10 -i ' + '%d.png '+ result_movie
p=subprocess.Popen(params,cwd=filepath1,shell=True)