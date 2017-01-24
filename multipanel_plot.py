"""This script generates a multipanel plot as seen in Fig. 2 of Fox et al. 2014."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,FormatStrFormatter
from astropy.io import fits
from astropy.coordinates import SkyCoord

#Files located in the folder
from rebin_ndarray import bin_ndarray
from readcol import *

#Global variables
c=2.99792e5 #Speed of light in km/s

#User input variables
Sightline = 'CTS487'

#Galactic l and b
l = np.deg2rad(np.asarray(SkyCoord.from_name(Sightline).galactic.l)); 
b = np.deg2rad(np.asarray(SkyCoord.from_name(Sightline).galactic.b))  #In radians

#Velocity correction factors
v_corr_LSR = 9.0*np.cos(l)*np.cos(b)+12.0*np.sin(l)*np.cos(b)+7.0*np.sin(b) #v_LSR of the object

#Read FITS files
g130m = np.genfromtxt('Sightlines/CTS487/CTS487_spec-G130M')
g160m = np.genfromtxt('Sightlines/CTS487/CTS487_spec-G160M')

#Read GASS HI spectra file
gass = np.genfromtxt('Sightlines/CTS487/CTS487_HIsp-GASS.dat')

#Extract Wavelength, Flux and Error
w130 = g130m[:,0]
f130 = g130m[:,1]
e130 = g130m[:,2]

w160 = g160m[:,0]
f160 = g160m[:,1]
e160 = g160m[:,2]

def array_truncater(array):
    """In order to rebin by 5 pixels, count the length of the arrays 
    and delete data from extremeties such that the length is 5.
    
    Input: array - enter a 1-D array.
    Output: array - returns the truncated input array such that it is a multiple of 5."""
    
    count = len(array) % 5
    
    while(count != 0):
        array = np.delete(array, [0])
        count = count - 1
        
    return array

"""--------------------------------------------------------------------------------------------------------------------------------------------------"""

#Truncate the arrays so that they can be rebinned by 5 pixels
w130 = array_truncater(w130); w160 = array_truncater(w160)
f130 = array_truncater(f130); f160 = array_truncater(f160)
e130 = array_truncater(e130); e160 = array_truncater(e160)

#Produces the entire spectra
#Combine both grating
w = np.concatenate((w130,w160))
f = np.concatenate((f130,f160))
e = np.concatenate((e130,e160))

#Two-dimensional array consisting of the wavelength and flux
a130 = np.array([w130,f130])
a160 = np.array([w160,f160])
a = np.array([w,f])

#Rebin data by 5 px
tmp130 = bin_ndarray(a130,(2,(len(w130)/5)),'avg')
tmp160 = bin_ndarray(a160,(2,(len(w160)/5)),'avg')
tmp = bin_ndarray(a,(2,(len(w)/5)),'avg')

"""----------------------------------------------------------------------------------------------"""
"""Produce panel plots"""
#http://matplotlib.org/examples/pylab_examples/subplots_demo.html
#http://stackoverflow.com/questions/6963035/pyplot-axes-labels-for-subplots

plt.clf()
f, axes = plt.subplots(6, 2, sharex='col', sharey='row',figsize = (15,15))

###Create HI 21cm spectra###

#Set tick marks for plot
xmajorLocator = MultipleLocator(200)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(100)
ymajorLocator = MultipleLocator(1)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator = MultipleLocator(0.08)

tmp1 = np.linspace(-400,400,10) #Create dashed line at 0
tmp2 = np.full((1,10),0)[0]

#Only apply GSR correction since tmp_hi[0] is already in LSR
#tmp_hi_GSR = tmp_hi[0] + v_corr_GSR

axes[0][1].plot(gass[:,1],gass[:,2],'k')
axes[0][1].plot(tmp1,tmp2,'k--')

axes[0][1].axis([-400,400,-.12,0.2])
axes[0][1].xaxis.set_major_locator(xmajorLocator)
axes[0][1].xaxis.set_major_formatter(xmajorFormatter)
axes[0][1].xaxis.set_major_locator(xminorLocator)
axes[0][1].yaxis.set_major_locator(ymajorLocator)
axes[0][1].yaxis.set_major_locator(yminorLocator)
axes[0][1].set_ylabel("T$_{B}$ [K]")
axes[0][1].text(-380,.2,"H I 21cm",fontsize=15,fontname="serif")

###Create Si III 1206###

#Set tick marks for plot
xmajorLocator = MultipleLocator(200)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(100)
ymajorLocator = MultipleLocator(.8)
yminorLocator = MultipleLocator(.4)

#Set these global variables according to the type of ions. Read from the qal.lst file. 
w0 = 1206.5 #Rest wavelength
fv = 1.66 #Oscillator strengh/f-value

#Set wavelength, velocity and flux
Lambda = w0
velocity = (tmp[0]-Lambda)/Lambda*c + v_corr_LSR #velocity = del_lambda/lambda*(speed of light in km/s); then apply  the LSR and GSR correction
flux = tmp[1]

#Fit the continuum
low1 = -500
high1 = -300
low2 = 230
high2 = 500
x1 = velocity[(velocity>=low1) & (velocity<=high1)] #lower bound range
x2 = velocity[(velocity>=low2) & (velocity<=high2)] #higher bound range
X = np.append(x1,x2)
y1 = flux[(velocity>=low1) & (velocity<=high1)] #lower bound range
y2 = flux[(velocity>=low2) & (velocity<=high2)] #higher bound range
Y = np.append(y1,y2)
Z = np.polyfit(X,Y,1) #Linear fit

#Generate data to plot continuum
xp = np.linspace(-500,501,len(flux[(flux>=-500) & (flux<=500)]))
p = np.poly1d(Z)

#Normalize flux
norm_flux = flux/p(xp)

tmp1 = np.linspace(-400,400,10)
tmp2 = np.full((1,10),1)[0]

axes[1][1].step(velocity,norm_flux,'k')
axes[1][1].plot(tmp1,tmp2,'k--')

axes[1][1].axis([-400,400,-.12,1.4])
axes[1][1].xaxis.set_major_locator(xmajorLocator)
axes[1][1].xaxis.set_major_formatter(xmajorFormatter)
axes[1][1].xaxis.set_major_locator(xminorLocator)
axes[1][1].yaxis.set_major_locator(ymajorLocator)
axes[1][1].yaxis.set_major_locator(yminorLocator)

axes[1][1].text(-380,.2,"Si III 1206",fontsize=15,fontname="serif")


###Create C II 1334###

#Set tick marks for plot
xmajorLocator = MultipleLocator(200)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(100)
ymajorLocator = MultipleLocator(.8)
yminorLocator = MultipleLocator(.4)

w0 = 1334.5323 #Rest wavelength
fv = 0.1278 #Oscillator strengh/f-value

#Set wavelength, velocity and flux
Lambda = w0
velocity = (tmp[0]-Lambda)/Lambda*c + v_corr_LSR #velocity = del_lambda/lambda*(speed of light in km/s); then apply  the LSR and GSR correction
flux = tmp[1]

#Fit the continuum
low1 = -500
high1 = -200
low2 = 350
high2 = 500
x1 = velocity[(velocity>=low1) & (velocity<=high1)] #lower bound range
x2 = velocity[(velocity>=low2) & (velocity<=high2)] #higher bound range
X = np.append(x1,x2)
y1 = flux[(velocity>=low1) & (velocity<=high1)] #lower bound range
y2 = flux[(velocity>=low2) & (velocity<=high2)] #higher bound range
Y = np.append(y1,y2)
Z = np.polyfit(X,Y,1) #Linear fit

#Generate data to plot continuum
xp = np.linspace(-500,501,len(flux[(flux>=-500) & (flux<=500)]))
p = np.poly1d(Z)

#Normalize flux
norm_flux = flux/p(xp)

tmp1 = np.linspace(-400,400,10)
tmp2 = np.full((1,10),1)[0]

axes[0][0].step(velocity,norm_flux,'k')
axes[0][0].plot(tmp1,tmp2,'k--')

axes[0][0].annotate('CII*',xy=(270,0.3),xytext=(270, 0.1), arrowprops=dict(arrowstyle="->"),color='black')

axes[0][0].axis([-400,400,-.12,1.4])
axes[0][0].xaxis.set_major_locator(xmajorLocator)
axes[0][0].xaxis.set_major_formatter(xmajorFormatter)
axes[0][0].xaxis.set_major_locator(xminorLocator)
axes[0][0].yaxis.set_major_locator(ymajorLocator)
axes[0][0].yaxis.set_major_locator(yminorLocator)

axes[0][0].text(-380,.2,"C II 1334",fontsize=15,fontname="serif")


###Create C IV 1548###

#Set tick marks for plot
xmajorLocator = MultipleLocator(200)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(100)
ymajorLocator = MultipleLocator(.8)
yminorLocator = MultipleLocator(.4)

#Set these global variables according to the type of ions. Read from the qal.lst file. 
w0 = 1548.195 #Rest wavelength
fv = 0.1908 #Oscillator strengh/f-value

#Set wavelength, velocity and flux
Lambda = w0
velocity = (tmp[0]-Lambda)/Lambda*c + v_corr_LSR #velocity = del_lambda/lambda*(speed of light in km/s); then apply  the LSR and GSR correction
flux = tmp[1]

#Fit the continuum
low1 = -500
high1 = -220
low2 = 200
high2 = 400
x1 = velocity[(velocity>=low1) & (velocity<=high1)] #lower bound range
x2 = velocity[(velocity>=low2) & (velocity<=high2)] #higher bound range
X = np.append(x1,x2)
y1 = flux[(velocity>=low1) & (velocity<=high1)] #lower bound range
y2 = flux[(velocity>=low2) & (velocity<=high2)] #higher bound range
Y = np.append(y1,y2)
Z = np.polyfit(X,Y,1) #Linear fit

#Generate data to plot continuum
xp = np.linspace(-500,501,len(flux[(flux>=-500) & (flux<=500)]))
p = np.poly1d(Z)

#Normalize flux
norm_flux = flux/p(xp)

tmp1 = np.linspace(-400,400,10)
tmp2 = np.full((1,10),1)[0]

axes[2][1].step(velocity,norm_flux,'k')
axes[2][1].plot(tmp1,tmp2,'k--')

axes[2][1].axis([-400,400,-.12,1.4])
axes[2][1].xaxis.set_major_locator(xmajorLocator)
axes[2][1].xaxis.set_major_formatter(xmajorFormatter)
axes[2][1].xaxis.set_major_locator(xminorLocator)
axes[2][1].yaxis.set_major_locator(ymajorLocator)
axes[2][1].yaxis.set_major_locator(yminorLocator)

axes[2][1].text(-380,.2,'C IV 1548',fontsize=15,fontname="serif")


###Create C IV 1550###

#Set tick marks for plot
xmajorLocator = MultipleLocator(200)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(100)
ymajorLocator = MultipleLocator(.8)
yminorLocator = MultipleLocator(.4)

#Set these global variables according to the type of ions. Read from the qal.lst file. 
w0 = 1550.770 #Rest wavelength
fv = 0.095220 #Oscillator strengh/f-value

#Set wavelength, velocity and flux
Lambda = w0
velocity = (tmp[0]-Lambda)/Lambda*c + v_corr_LSR #velocity = del_lambda/lambda*(speed of light in km/s); then apply  the LSR and GSR correction
flux = tmp[1]

#Fit the continuum
low1 = -400
high1 = -180
low2 = 100
high2 = 500
x1 = velocity[(velocity>=low1) & (velocity<=high1)] #lower bound range
x2 = velocity[(velocity>=low2) & (velocity<=high2)] #higher bound range
X = np.append(x1,x2)
y1 = flux[(velocity>=low1) & (velocity<=high1)] #lower bound range
y2 = flux[(velocity>=low2) & (velocity<=high2)] #higher bound range
Y = np.append(y1,y2)
Z = np.polyfit(X,Y,1) #Linear fit

#Generate data to plot continuum
xp = np.linspace(-500,501,len(flux[(flux>=-500) & (flux<=500)]))
p = np.poly1d(Z)

#Normalize flux
norm_flux = flux/p(xp)

tmp1 = np.linspace(-400,400,10)
tmp2 = np.full((1,10),1)[0]

axes[3][1].step(velocity,norm_flux,'k')
axes[3][1].plot(tmp1,tmp2,'k--')

axes[3][1].axis([-400,400,-.12,1.4])
axes[3][1].xaxis.set_major_locator(xmajorLocator)
axes[3][1].xaxis.set_major_formatter(xmajorFormatter)
axes[3][1].xaxis.set_major_locator(xminorLocator)
axes[3][1].yaxis.set_major_locator(ymajorLocator)
axes[3][1].yaxis.set_major_locator(yminorLocator)

axes[3][1].text(-380,.2,'C IV 1550',fontsize=15,fontname="serif")


###Create Al II 1670###

#Set tick marks for plot
xmajorLocator = MultipleLocator(200)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(100)
ymajorLocator = MultipleLocator(.8)
yminorLocator = MultipleLocator(.4)

#Set wavelength, velocity and flux
Lambda = w0
velocity = (tmp[0]-Lambda)/Lambda*c + v_corr_LSR #velocity = del_lambda/lambda*(speed of light in km/s); then apply  the LSR and GSR correction
flux = tmp[1]

#Fit the continuum
low1 = -500
high1 = -200
low2 = 180
high2 = 500
x1 = velocity[(velocity>=low1) & (velocity<=high1)] #lower bound range
x2 = velocity[(velocity>=low2) & (velocity<=high2)] #higher bound range
X = np.append(x1,x2)
y1 = flux[(velocity>=low1) & (velocity<=high1)] #lower bound range
y2 = flux[(velocity>=low2) & (velocity<=high2)] #higher bound range
Y = np.append(y1,y2)
Z = np.polyfit(X,Y,1) #Linear fit

#Generate data to plot continuum
xp = np.linspace(-500,501,len(flux[(flux>=-500) & (flux<=500)]))
p = np.poly1d(Z)

#Normalize flux
norm_flux = flux/p(xp)

tmp1 = np.linspace(-400,400,10)
tmp2 = np.full((1,10),1)[0]

axes[1][0].step(velocity,norm_flux,'k')
axes[1][0].plot(tmp1,tmp2,'k--')

axes[1][0].axis([-400,400,-.12,1.4])
axes[1][0].xaxis.set_major_locator(xmajorLocator)
axes[1][0].xaxis.set_major_formatter(xmajorFormatter)
axes[1][0].xaxis.set_major_locator(xminorLocator)
axes[1][0].yaxis.set_major_locator(ymajorLocator)
axes[1][0].yaxis.set_major_locator(yminorLocator)

axes[1][0].text(-380,.2,"Al II 1670",fontsize=15,fontname="serif")


#axes[5][1].xlabel('V$_LSR$', fontsize = 15)

###Create Si II 1190###

#Set tick marks for plot
xmajorLocator = MultipleLocator(200)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(100)
ymajorLocator = MultipleLocator(.8)
yminorLocator = MultipleLocator(.4)

#Set these global variables according to the type of ions. Read from the qal.lst file. 
w0 = 1190.4158 #Rest wavelength
fv = 0.25020 #Oscillator strengh/f-value

#Set wavelength, velocity and flux
Lambda = w0
velocity = (tmp[0]-Lambda)/Lambda*c + v_corr_LSR #velocity = del_lambda/lambda*(speed of light in km/s); then apply  the LSR and GSR correction
flux = tmp[1]

#Fit the continuum
low1 = -500
high1 = -200
low2 = 200
high2 = 500
x1 = velocity[(velocity>=low1) & (velocity<=high1)] #lower bound range
x2 = velocity[(velocity>=low2) & (velocity<=high2)] #higher bound range
X = np.append(x1,x2)
y1 = flux[(velocity>=low1) & (velocity<=high1)] #lower bound range
y2 = flux[(velocity>=low2) & (velocity<=high2)] #higher bound range
Y = np.append(y1,y2)
Z = np.polyfit(X,Y,1) #Linear fit

#Generate data to plot continuum
xp = np.linspace(-500,501,len(flux[(flux>=-500) & (flux<=500)]))
p = np.poly1d(Z)

#Normalize flux
norm_flux = flux/p(xp)

tmp1 = np.linspace(-400,400,10)
tmp2 = np.full((1,10),1)[0]

axes[4][0].step(velocity,norm_flux,'k')
axes[4][0].plot(tmp1,tmp2,'k--')

axes[4][0].axis([-400,400,-.12,1.4])
axes[4][0].xaxis.set_major_locator(xmajorLocator)
axes[4][0].xaxis.set_major_formatter(xmajorFormatter)
axes[4][0].xaxis.set_major_locator(xminorLocator)
axes[4][0].yaxis.set_major_locator(ymajorLocator)
axes[4][0].yaxis.set_major_locator(yminorLocator)

axes[4][0].text(-380,.2,"Si II 1190",fontsize=15,fontname="serif")


###Create Si II 1193###

#Set tick marks for plot
xmajorLocator = MultipleLocator(200)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(100)
ymajorLocator = MultipleLocator(.8)
yminorLocator = MultipleLocator(.4)


#Set wavelength, velocity and flux
Lambda = w0
velocity = (tmp[0]-Lambda)/Lambda*c + v_corr_LSR #velocity = del_lambda/lambda*(speed of light in km/s); then apply  the LSR and GSR correction
flux = tmp[1]

#Fit the continuum
low1 = -500
high1 = -200
low2 = 200
high2 = 500
x1 = velocity[(velocity>=low1) & (velocity<=high1)] #lower bound range
x2 = velocity[(velocity>=low2) & (velocity<=high2)] #higher bound range
X = np.append(x1,x2)
y1 = flux[(velocity>=low1) & (velocity<=high1)] #lower bound range
y2 = flux[(velocity>=low2) & (velocity<=high2)] #higher bound range
Y = np.append(y1,y2)
Z = np.polyfit(X,Y,1) #Linear fit

#Generate data to plot continuum
xp = np.linspace(-500,501,len(flux[(flux>=-500) & (flux<=500)]))
p = np.poly1d(Z)

#Normalize flux
norm_flux = flux/p(xp)

tmp1 = np.linspace(-400,400,10)
tmp2 = np.full((1,10),1)[0]

axes[3][0].step(velocity,norm_flux,'k')
axes[3][0].plot(tmp1,tmp2,'k--')

axes[3][0].axis([-400,400,-.12,1.4])
axes[3][0].xaxis.set_major_locator(xmajorLocator)
axes[3][0].xaxis.set_major_formatter(xmajorFormatter)
axes[3][0].xaxis.set_major_locator(xminorLocator)
axes[3][0].yaxis.set_major_locator(ymajorLocator)
axes[3][0].yaxis.set_major_locator(yminorLocator)

axes[3][0].text(-380,.2,"Si II 1193",fontsize=15,fontname="serif")


###Create Si II 1260###

#Set tick marks for plot
xmajorLocator = MultipleLocator(200)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(100)
ymajorLocator = MultipleLocator(.8)
yminorLocator = MultipleLocator(.4)

#Set these global variables according to the type of ions. Read from the qal.lst file. 
w0 = 1260.4221 #Rest wavelength
fv = 1.007 #Oscillator strengh/f-value

#Set wavelength, velocity and flux
Lambda = w0
velocity = (tmp[0]-Lambda)/Lambda*c + v_corr_LSR #velocity = del_lambda/lambda*(speed of light in km/s); then apply  the LSR and GSR correction
flux = tmp[1]

#Fit the continuum
low1 = -500
high1 = -350
low2 = 200
high2 = 500
x1 = velocity[(velocity>=low1) & (velocity<=high1)] #lower bound range
x2 = velocity[(velocity>=low2) & (velocity<=high2)] #higher bound range
X = np.append(x1,x2)
y1 = flux[(velocity>=low1) & (velocity<=high1)] #lower bound range
y2 = flux[(velocity>=low2) & (velocity<=high2)] #higher bound range
Y = np.append(y1,y2)
Z = np.polyfit(X,Y,1) #Linear fit

#Generate data to plot continuum
xp = np.linspace(-500,501,len(flux[(flux>=-500) & (flux<=500)]))
p = np.poly1d(Z)

#Normalize flux
norm_flux = flux/p(xp)

tmp1 = np.linspace(-400,400,10)
tmp2 = np.full((1,10),1)[0]

axes[2][0].step(velocity,norm_flux,'k')
axes[2][0].plot(tmp1,tmp2,'k--')

axes[2][0].annotate('Si II',xy=(-200,0.3),xytext=(-200, 0.1), arrowprops=dict(arrowstyle="->"),color='black')

axes[2][0].axis([-400,400,-.12,1.4])
axes[2][0].xaxis.set_major_locator(xmajorLocator)
axes[2][0].xaxis.set_major_formatter(xmajorFormatter)
axes[2][0].xaxis.set_major_locator(xminorLocator)
axes[2][0].yaxis.set_major_locator(ymajorLocator)
axes[2][0].yaxis.set_major_locator(yminorLocator)

axes[2][0].text(-380,.15,"Si II 1260",fontsize=15,fontname="serif")



###Create Si II 1526###

#Set tick marks for plot
xmajorLocator = MultipleLocator(200)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator = MultipleLocator(100)
ymajorLocator = MultipleLocator(.8)
yminorLocator = MultipleLocator(.4)

#Set these global variables according to the type of ions. Read from the qal.lst file. 
w0 = 1526.7066 #Rest wavelength
fv = 0.127 #Oscillator strengh/f-value

#Set wavelength, velocity and flux
Lambda = w0
velocity = (tmp[0]-Lambda)/Lambda*c + v_corr_LSR #velocity = del_lambda/lambda*(speed of light in km/s); then apply  the LSR and GSR correction
flux = tmp[1]

#Fit the continuum
low1 = -500
high1 = -200
low2 = 180
high2 = 500
x1 = velocity[(velocity>=low1) & (velocity<=high1)] #lower bound range
x2 = velocity[(velocity>=low2) & (velocity<=high2)] #higher bound range
X = np.append(x1,x2)
y1 = flux[(velocity>=low1) & (velocity<=high1)] #lower bound range
y2 = flux[(velocity>=low2) & (velocity<=high2)] #higher bound range
Y = np.append(y1,y2)
Z = np.polyfit(X,Y,1) #Linear fit

#Generate data to plot continuum
xp = np.linspace(-500,501,len(flux[(flux>=-500) & (flux<=500)]))
p = np.poly1d(Z)

#Normalize flux
norm_flux = flux/p(xp)

tmp1 = np.linspace(-400,400,10)
tmp2 = np.full((1,10),1)[0]

axes[5][0].step(velocity,norm_flux,'k')
axes[5][0].plot(tmp1,tmp2,'k--')

axes[5][0].axis([-400,400,-.12,1.4])
axes[5][0].xaxis.set_major_locator(xmajorLocator)
axes[5][0].xaxis.set_major_formatter(xmajorFormatter)
axes[5][0].xaxis.set_major_locator(xminorLocator)
axes[5][0].yaxis.set_major_locator(ymajorLocator)
axes[5][0].yaxis.set_major_locator(yminorLocator)

axes[5][0].text(-380,.2,"Si II 1526",fontsize=15,fontname="serif")

#axes[5][0].xlabel('V$_LSR$', fontsize = 15)
"""--------------------------------------------------------------------------------------"""
f.subplots_adjust(hspace=0, wspace=0.08) #Set no space between subplots
plt.setp([a.get_xticklabels() for a in f.axes[:-2]], visible=False) #Bottom two panels show axis label

#Common y-axis label
f.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.ylabel('Normalized Flux', fontsize = 20)
plt.xlabel('V$_{LSR}$', fontsize = 20)

filename = "Panel_plot_" + Sightline + "_LSR.pdf"
plt.savefig(filename, format = "pdf", bbox_inches='tight')