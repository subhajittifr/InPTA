import numpy as np
from numpy import pi, cos, arcsin, sqrt, sin, log, sign
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib
import os, sys


A_X=np.loadtxt(sys.argv[1],dtype=object)
A_Y=np.loadtxt(sys.argv[2],dtype=object)




def get_ra_dec(sample):
    RA_DEC=[]
    for i in range(len(sample)):
        p=sample[i]
        a1=int(p[1:3])
        b1=int(p[3:5])
        e1=int(p[5:8])
        d1=int(p[8:10])
        psrra=     15*(a1+b1/60)*pi/180
        psrdec =   (e1+e1/np.abs(e1)*d1/60)*pi/180
        RA_DEC.append((sample[i],psrra,psrdec))
    return np.array(RA_DEC)





def get_angle(sample):
    s,info=0,[]
    z=get_ra_dec(sample)
    for i in range(len(z)):
        for j in range(i+1,len(z)):
            ra_1=float(z[i][1])
            ra_2=float(z[j][1])
            dec_1=float(z[i][2])
            dec_2=float(z[j][2])
            theta_12=2*arcsin(sqrt(sin((dec_2-dec_1)/2)**2+cos(dec_1)*cos(dec_2)*sin((ra_1-ra_2)/2)**2))
            info.append([z[i][0],z[j][0],theta_12*180/pi])
            s+=1
    return np.array(info,dtype='object')
            


inpta_X=get_angle(A_X)[...,2].ravel()
inpta_Y=get_angle(A_Y)[...,2].ravel()


outfile='hist_'+sys.argv[1][0:7]+'_'+sys.argv[2][0:7]+'.pdf'


plt.figure(figsize=(8,6))

binwidth=5
plt.hist(inpta_X,bins=range(0, 180, binwidth),alpha=0.3,label='InPTA '+str(sys.argv[1][0:7]), fill=False, hatch='XX', ec="red", lw=2)
plt.hist(inpta_Y,bins=range(0, 180, binwidth),range=(0,180),alpha=0.4,label='InPTA '+str(sys.argv[2][0:7]), fill=False, hatch='...', ec="green")
plt.xticks(np.arange(0.0, 190.0, 20), fontsize=14)
plt.yticks(np.arange(0, 20, 2), fontsize=14)
plt.xlabel("Angular separation", fontsize=15)
plt.ylabel("Number of pairs", fontsize=15)
plt.legend(loc='upper right', fontsize=15)
plt.savefig(outfile)
plt.show()







