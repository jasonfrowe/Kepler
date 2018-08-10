import numpy as np 
from numpy import zeros
from numpy import ones
import tfit5
import fittransitmodel as ftf
import matplotlib  #ploting
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
import math   #used for floor command

def readtt(*files):
    "reading in TT files"
    nmax=0 #we will first scan through the files to determine what size of array we need.
    for filename in files:
        if filename == 'null':
            i=0
        else:
            f = open(filename,'r')
            i=0
            for line in f:
                i+=1
            f.close()
        
        nmax=max(i,nmax) #save the largest number.
        
    np=len(files) #number of planets to read in
    ntt =zeros(np, dtype="int32") #allocate array for number of TTVs measured
    tobs=zeros(shape=(np,nmax)) #allocate array for timestamps
    omc =zeros(shape=(np,nmax)) #allocate array for O-C
    
    i=-1 #counter for files scanned 
    for filename in files: #scan through all files from input
        i+=1
        if filename == 'null':
            ntt[i]=0
        else:
            f = open(filename,'r')
            j=-1 #counter for valid O-C read 
            for line in f:
                line = line.strip() #get rid of line breaks
                columns = line.split()
                if float(columns[2]) > 0.0 :
                    j+=1
                    tobs[i,j]=float(columns[0])
                    omc[i,j]=float(columns[1])
            ntt[i]=j+1
    return ntt, tobs, omc;
    

def transitplot(time,flux,sol,nplanetplot=1, itime=-1, ntt=0, tobs=0, omc=0, dtype=0):
    "plot the transit model"
    nplanet=int((len(sol)-8)/10) #number of planets

    #deal with input vars and translating to FORTRAN friendly.
    
    if type(itime) is int :
        if itime < 0 :
            itime=ones(len(time))*0.020434
        else:
            itime=ones(len(time))*float(itime)
    
    if type(ntt) is int :
        nttin=  zeros(nplanet, dtype="int32") #number of TTVs measured 
        tobsin= zeros(shape=(nplanet,len(time))) #time stamps of TTV measurements (days)
        omcin=  zeros(shape=(nplanet,len(time))) #TTV measurements (O-C) (days)
    else:
        nttin=ntt
        tobsin=tobs
        omcin=omc
 
    if type(dtype) is int :
        dtypein=ones(len(time), dtype="int32")*int(dtype) #contains data type, 0-photometry,1=RV data
    
    #remove other planets for plotting
    sol2=np.copy(sol)
    for i in range(1,nplanet+1):
        if i!=nplanetplot:
            nc=8+10*(i-1)
            sol2[nc+3]=0.0 #rdrs
    tmodel= zeros(len(time)) #contains the transit model
    tfit5.transitmodel(nplanet,sol2,time,itime,nttin,tobsin,omcin,tmodel,dtypein)
    stdev=np.std(flux-tmodel)
    #print(stdev)
    
    #make a model with only the other transits to subtract
    nc=8+10*(nplanetplot-1)
    sol2=np.copy(sol)
    sol2[nc+3]=0.0 #rdrs
    tmodel2= zeros(len(time)) #contains the transit model
    tfit5.transitmodel(nplanet,sol2,time,itime,nttin,tobsin,omcin,tmodel2,dtypein)
    
    epo=sol[nc+0] #time of center of transit
    per=sol[nc+1] #orbital period
    zpt=sol[7] #photometric zero-point
    
    tdur=tfit5.transitdur(sol,1)/3600.0 #transit duration in hours
    
    ph1=epo/per-math.floor(epo/per) #calculate phases
    phase=[]
    tcor=tfit5.lininterp(tobsin,omcin,nplanetplot,nttin,epo)
    #print(tcor,nttin,tobsin[1,1],omcin[1,1])
    for x in time:
        if nttin[nplanetplot-1] > 0:
            tcor=tfit5.lininterp(tobsin,omcin,nplanetplot,nttin,x)
        else:
            tcor=0.0
        t=x-tcor
        ph=(t/per-math.floor(t/per)-ph1)*per*24.0 #phase in hours offset to zero.
        phase.append(ph)
        
    phase = np.array(phase) #convert from list to array
    
    phasesort=np.copy(phase)
    fluxsort=np.copy(tmodel)
    p=ones(len(phase), dtype="int32") #allocate array for output. FORTRAN needs int32 
    tfit5.rqsort(phase,p)    
    i1=0
    i2=len(phase)-1
    for i in range(0,len(phase)):
        phasesort[i]=phase[p[i]-1]
        fluxsort[i]=tmodel[p[i]-1]
        if phasesort[i] < -tdur:
            i1=i
        if phasesort[i] <  tdur:
            i2=i
    
    fplot=flux-tmodel2+1.0
    
    plt.figure(figsize=(12,10)) #adjust size of figure
    matplotlib.rcParams.update({'font.size': 22}) #adjust font
    plt.scatter(phase, fplot, c="blue", s=100.0, alpha=0.35, edgecolors="none") #scatter plot
    plt.plot(phasesort, fluxsort, c="red", lw=3.0)
    plt.xlabel('Phase (hours)') #x-label
    plt.ylabel('Relative Flux') #y-label
    x1,x2,y1,y2 = plt.axis()    #get range of plot
    ymin=min(fluxsort[i1:i2])
    ymax=max(fluxsort[i1:i2])
    y1=ymin-0.1*(ymax-ymin)-2.0*stdev
    y2=ymax+0.1*(ymax-ymin)+2.0*stdev
    plt.axis((-tdur,tdur,y1,y2)) #readjust range
    plt.show()  #show the plot
        
    return;
        
    return;


def transitmodel (sol,time, itime=-1, ntt=0, tobs=0, omc=0, dtype=0 ): 
    "read in transitmodel solution"  
    nplanet=int((len(sol)-8)/10) #number of planets
    
    if type(itime) is float :
        if itime < 0 :
            itime=ones(len(time))*0.020434
        else:
            itime=ones(len(time))*float(itime)
    
    if type(ntt) is int :
        nttin=  zeros(nplanet, dtype="int32") #number of TTVs measured 
        tobsin= zeros(shape=(nplanet,len(time))) #time stamps of TTV measurements (days)
        omcin=  zeros(shape=(nplanet,len(time))) #TTV measurements (O-C) (days)
    else:
        nttin=ntt
        tobsin=tobs
        omcin=omc
    
    if type(dtype) is int :
        dtypein=ones(len(time), dtype="int32")*int(dtype) #contains data type, 0-photometry,1=RV data

    tmodel= zeros(len(time)) #contains the transit model
    tfit5.transitmodel(nplanet,sol,time,itime,nttin,tobsin,omcin,tmodel,dtypein)
    return tmodel;

def readphotometry (filename):
    "reading in Kepler photometry"
    time=[]  #initialize arrays
    flux=[]
    ferr=[]
    f = open(filename, 'r')
    for line in f:
        line = line.strip() #get rid of the \n at the end of the line
        columns = line.split() #break into columns
        time.append(float(columns[0])-54900.0+0.5) #correct for file zero-points to get BJD-2454900
        flux.append(float(columns[1])+1.0) #photometry
        ferr.append(float(columns[2])) #photometric uncertainty 
    f.close()
    time = np.array(time)
    flux = np.array(flux)
    ferr = np.array(ferr)
    return time, flux, ferr;

def readsol (filename):
    "read in transitmodel solution"    
    nplanetmax=9 #maximum number of planets that an n0.dat file can handle
    nplanet=0 #count number of planets found in the solution
    solin=zeros(nplanetmax*10+8) #allocate array to hold parameters. init to zero.
    serrin=zeros(nplanetmax*10+8)
    f = open(filename, 'r')
    for line in f:
        line = line.strip() #get rid of the \n at the end of the line
        columns = line.split() #break into columns
        if columns[0][0:3]=='RHO':
            solin[0]=columns[1]
            serrin[0]=columns[3]
        elif columns[0][0:3]=='NL1':
            solin[1]=columns[1]
            serrin[1]=columns[3]
        elif columns[0][0:3]=='NL2':
            solin[2]=columns[1]
            serrin[2]=columns[3]
        elif columns[0][0:3]=='NL3':
            solin[3]=columns[1]
            serrin[3]=columns[3]
        elif columns[0][0:3]=='NL4':
            solin[4]=columns[1]
            serrin[4]=columns[3]
        elif columns[0][0:3]=='DIL':
            solin[5]=columns[1]
            serrin[5]=columns[3]
        elif columns[0][0:3]=='VOF':
            solin[6]=columns[1]
            serrin[6]=columns[3]
        elif columns[0][0:3]=='ZPT':
            solin[7]=columns[1]
            serrin[7]=columns[3]
        elif columns[0][0:2]=='EP':
            np=float(columns[0][2])
            if np>nplanet:
                nplanet=np
            j=int(10*(np-1)+8+0)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='PE':
            np=float(columns[0][2])
            if np>nplanet:
                nplanet=np
            j=int(10*(np-1)+8+1)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='BB':
            np=float(columns[0][2])
            if np>nplanet:
                nplanet=np
            j=int(10*(np-1)+8+2)
            solin[j]=columns[1] 
            serrin[j]=columns[3]
        elif columns[0][0:2]=='RD':
            np=float(columns[0][2])
            if np>nplanet:
                nplanet=np
            j=int(10*(np-1)+8+3)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='EC':
            np=float(columns[0][2])
            if np>nplanet:
                nplanet=np
            j=int(10*(np-1)+8+4)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='ES':
            np=float(columns[0][2])
            if np>nplanet:
                nplanet=np
            j=int(10*(np-1)+8+5)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='KR':
            np=float(columns[0][2])
            if np>nplanet:
                nplanet=np
            j=int(10*(np-1)+8+6)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='TE':
            np=float(columns[0][2])
            if np>nplanet:
                nplanet=np
            j=int(10*(np-1)+8+7)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='EL':
            np=float(columns[0][2])
            if np>nplanet:
                nplanet=np
            j=int(10*(np-1)+8+8)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='AL':
            np=float(columns[0][2])
            if np>nplanet:
                nplanet=np
            j=int(10*(np-1)+8+9)
            solin[j]=columns[1]
            serrin[j]=columns[3]
    f.close()
    #print(nplanet)
    sol=solin[0:int(nplanet*10+8)]
    serr=serrin[0:int(nplanet*10+8)]
    return sol, serr;

def fittrmodel(time,flux,ferr,sol,serr,itime=-1,ntt=0,tobs=0,omc=0,dtype=0):
    nfit=108 
    npt=len(time)
    solin=zeros(nfit)
    solin[0:len(sol)]=sol
    serrin=zeros(shape=(nfit,2))
    serrin[0:len(serr),1]=serr
    nplanet=int((len(sol)-8)/10) #number of planets
    
    if type(itime) is int :
        if itime < 0 :
            itime=ones(len(time))*0.020434
        else:
            itime=ones(len(time))*float(itime)

    if type(ntt) is int :
        nttin=  zeros(nplanet, dtype="int32") #number of TTVs measured 
        tobsin= zeros(shape=(nplanet,len(time))) #time stamps of TTV measurements (days)
        omcin=  zeros(shape=(nplanet,len(time))) #TTV measurements (O-C) (days)
    else:
        nttin=ntt
        tobsin=tobs
        omcin=omc
            
    if type(dtype) is int :
        dtypein=ones(len(time), dtype="int32")*int(dtype) #contains data type, 0-photometry,1=RV data
        
    nfrho=1 #disable mean stellar density 
    rhoi=np.float(0.0)
    rhoierr=zeros(9) #contains stellar density prior (currently disabled)
    
    #fittransitmodel3(nfit,sol,serr,nplanet,npt,at,am,ae,ait,dtype,ntt,tobs,omc,nfrho,rhoi,rhoierr)
    ftf.fittransitmodel3(nfit,solin,serrin,nplanet,npt,time,flux,ferr,itime,dtypein,nttin,tobsin,\
        omcin,nfrho,rhoi,rhoierr)
    
    solout=zeros(len(sol))
    solout=solin[0:len(sol)]
    
    return solout;