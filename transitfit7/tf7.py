import numpy as np #We will make extensive use of Numpy arrays 
from scipy import stats #For Kernel Density Estimation
import matplotlib  #ploting
#matplotlib.use("Agag")  #some hack to stop the bouncing python icon when plotting of os/x
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import ScalarFormatter
import lcmodel as lc


def triplot(chain,burnin,label,colour,nbin,ntick=5):
    "Making a Triangle Plot"

    matplotlib.rcParams.update({'font.size': 10}) #adjust font
    nullfmt = NullFormatter()       # removing tick labels
    deffmt = ScalarFormatter()      # adding tick labels
    n=len(chain[1,:])               # determine number of subplots
    plt.figure(1, figsize=(20, 20))   # make a square figure
    wsize=0.9/n                     # size of individual panels

    prange=np.zeros(shape=(n,2))    #store range of parameters in chains
    for i in range(0,n):
        prange[i,0]=min(chain[burnin:,i]) #find minimum values
        prange[i,1]=max(chain[burnin:,i]) #find maximum values

    for j in range(0,n):        #loop over each variable
        for i in range(j,n):    #loop again over each variable
            left, width = 0.1+j*wsize, wsize        #determine panel size: left position and width
            bottom, height = 0.9-i*wsize, wsize     #determine panel size: bottom position and width
            rect_scatter = [left, bottom, width, height]#save panel size position
            axScatter = plt.axes(rect_scatter)          #set panel size

            #put histogram on diagonals, scatter plot otherwise
            if i == j:

                plt.hist(chain[burnin:,j],nbin,histtype='stepfilled', density=True, \
                 facecolor=colour[i], alpha=0.8)

                #make a uniform sample across the parameter range
                x_eval = np.linspace(prange[j,0], prange[j,1], num=100)
                kde1 = stats.gaussian_kde(chain[burnin:,j],0.3) #Kernel Density Estimate
                #overlay the KDE
                plt.plot(x_eval, kde1(x_eval), 'k-', lw=3)

                x1,x2,y1,y2 = plt.axis()
                plt.axis((prange[j,0],prange[j,1],y1,y2))
            else:
                axScatter.scatter(chain[burnin:,j],chain[burnin:,i],c="black", s=0.2, alpha=0.1, \
                 edgecolors="none")
                plt.axis((prange[j,0],prange[j,1],prange[i,0],prange[i,1]))

            #to use a sensible x-tick range, we have to do it manually.
            dpr=(prange[j,1]-prange[j,0])/ntick         #make ntick ticks
            rr=np.power(10.0,np.floor(np.log10(dpr)))
            npr=np.floor(dpr/rr+0.5)*rr
            plt.xticks(np.arange(np.floor(prange[j,0]/rr)*rr+npr, \
             np.floor(prange[j,1]/rr)*rr,npr),rotation=30) #make ticks

            if i != j:
                #use a sensible y-tick range, we have to do it manually
                dpr=(prange[i,1]-prange[i,0])/ntick         #make ntick ticks
                rr=np.power(10.0,np.floor(np.log10(dpr)))
                npr=np.floor(dpr/rr+0.5)*rr
                plt.yticks(np.arange(np.floor(prange[i,0]/rr)*rr+npr,\
                 np.floor(prange[i,1]/rr)*rr,npr),rotation=0) #make ticks

            axScatter.xaxis.set_major_formatter(nullfmt)  #default is to leave off tick mark labels
            axScatter.yaxis.set_major_formatter(nullfmt)

            #if we are on the sides, add tick and axes labels
            if i==n-1:
                axScatter.xaxis.set_major_formatter(deffmt)
                plt.xlabel(label[j])
            if j==0 :
                if i > 0:
                    axScatter.yaxis.set_major_formatter(deffmt)
                    plt.ylabel(label[i])


    plt.show()

    return;

def plotchains(chain,label,colour,burnin,psize=0.1):
    nullfmt = NullFormatter()
    matplotlib.rcParams.update({'font.size': 10}) #adjust font
    plt.figure(figsize=(12,37)) #adjust size of figure

    x=np.arange(burnin+1,len(chain)+1,1)
    npar=len(chain[0,:])
    for i in range(0,npar):
        axScatter=plt.subplot(npar, 1, i+1)
        plt.scatter(x,chain[burnin:,i],c=colour[i],s=psize,alpha=0.1)  #plot parameter a
        plt.ylabel(label[i])                   #y-label

        x1,x2,y1,y2 = plt.axis()
        y1=np.min(chain[burnin:,i])
        y2=np.max(chain[burnin:,i])
        plt.axis((x1,x2,y1,y2))


        if i < npar-1:
            axScatter.xaxis.set_major_formatter(nullfmt)

    plt.xlabel('Iteration')           #x-label

    plt.show()

def histplots(chain,sol,serr,burnin,nbin,label,colour):

    matplotlib.rcParams.update({'font.size': 10}) #adjust font
    plt.figure(1, figsize=(20, 20)) 
    nullfmt = NullFormatter()

    nbodies=int((len(sol)-7)/7)
    npl=7+nbodies*7

    n1=int(nbodies)
    n2=int(7)

    jpar=-1
    for j in range(7):
        if np.abs(serr[j]) > 1.0e-30:
            jpar=jpar+1
            npanel=int(jpar+1)
            axScatter=plt.subplot(n1, n2, npanel)

            minx=np.min(chain[burnin:,jpar]) #range of parameter
            maxx=np.max(chain[burnin:,jpar])

            x_eval = np.linspace(minx, maxx, num=100) #make a uniform sample across the parameter range

            kde1 = stats.gaussian_kde(chain[burnin:,jpar],0.3) #Kernel Density Estimate

            #plot the histogram
            plt.hist(chain[burnin:,jpar],nbin,histtype='stepfilled', density=True, facecolor=colour[jpar],\
                     alpha=0.6)

            #overlay the KDE
            plt.plot(x_eval, kde1(x_eval), 'k-', lw=3)

            plt.xlabel(label[jpar])
            #plt.ylabel('Probability Density')
            axScatter.yaxis.set_major_formatter(nullfmt)

    #Dscale parameter
    npanel=npanel+1
    axScatter=plt.subplot(n1, n2, npanel)
    k=chain.shape[1]-1
    minx=np.min(chain[burnin:,k]) #range of parameter
    maxx=np.max(chain[burnin:,k])
    x_eval = np.linspace(minx, maxx, num=100) #make a uniform sample across the parameter range
    kde1 = stats.gaussian_kde(chain[burnin:,k],0.3) #Kernel Density Estimate
    #plot the histogram
    plt.hist(chain[burnin:,k],nbin,histtype='stepfilled', density=True, facecolor=colour[k], alpha=0.6)
    #overlay the KDE
    plt.plot(x_eval, kde1(x_eval), 'k-', lw=3)
    plt.xlabel(label[k])
    #plt.ylabel('Probability Density')
    axScatter.yaxis.set_major_formatter(nullfmt)

    for i in range(nbodies):
        for j in range(7):
            n=7+i*7+j
            #print(n,j+7)
            if np.abs(serr[n]) > 1.0e-30:
#                print(i*n1+j,j)
                jpar=jpar+1
                npanel=int(i*n2+j+1)
                axScatter=plt.subplot(n1, n2, npanel)

                minx=np.min(chain[burnin:,jpar]) #range of parameter
                maxx=np.max(chain[burnin:,jpar])

                x_eval = np.linspace(minx, maxx, num=100) #make a uniform sample across the parameter range

                kde1 = stats.gaussian_kde(chain[burnin:,jpar],0.3) #Kernel Density Estimate

                #plot the histogram
                plt.hist(chain[burnin:,jpar],nbin,histtype='stepfilled', density=True, facecolor=colour[jpar],\
                     alpha=0.6)

                #overlay the KDE
                plt.plot(x_eval, kde1(x_eval), 'k-', lw=3)

                plt.xlabel(label[jpar])
                #plt.ylabel('Probability Density')
                axScatter.yaxis.set_major_formatter(nullfmt)
    
    plt.show()

def readphotometry_v2(filename):
    "reading in Kepler photometry"
    time=[]  #initialize arrays
    flux=[]
    ferr=[]
    itime=[]
    f = open(filename, 'r')
    for line in f:
        line = line.strip() #get rid of the \n at the end of the line
        columns = line.split() #break into columns
        time.append(float(columns[0])-54900.0+0.5) #correct for file zero-points to get BJD-2454900
        flux.append(float(columns[1])+1.0) #photometry
        ferr.append(float(columns[2])) #photometric uncertainty
        if len(columns)>=4:
            itime.append(float(columns[3]))
        else:
            itime.append(float(0.0204340))
    f.close()
    time = np.array(time)
    flux = np.array(flux)
    ferr = np.array(ferr)
    itime = np.array(itime)
    return time, flux, ferr, itime;

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
    solin=np.zeros(nplanetmax*10+8) #allocate array to hold parameters. init to zero.
    serrin=np.zeros(nplanetmax*10+8)
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
        elif columns[0][0:3]=='ZPT':
            solin[6]=columns[1]
            serrin[6]=columns[3]

        elif columns[0][0:2]=='EP':
            npl=float(columns[0][2])
            if npl>nplanet:
                nplanet=npl
            j=int(7*(npl-1)+7+0)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='PE':
            npl=float(columns[0][2])
            if npl>nplanet:
                nplanet=npl
            j=int(7*(npl-1)+7+1)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='BB':
            npl=float(columns[0][2])
            if npl>nplanet:
                nplanet=npl
            j=int(7*(npl-1)+7+2)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='RD':
            npl=float(columns[0][2])
            if npl>nplanet:
                nplanet=npl
            j=int(7*(npl-1)+7+3)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='MA':
            npl=float(columns[0][2])
            if npl>nplanet:
                nplanet=npl
            j=int(7*(npl-1)+7+4)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='EC':
            npl=float(columns[0][2])
            if npl>nplanet:
                nplanet=npl
            j=int(7*(npl-1)+7+5)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='ES':
            npl=float(columns[0][2])
            if npl>nplanet:
                nplanet=npl
            j=int(7*(npl-1)+7+6)
            solin[j]=columns[1]
            serrin[j]=columns[3]
    f.close()
    #print(nplanet)
    sol=solin[0:int(nplanet*7+7)]
    serr=serrin[0:int(nplanet*7+7)]
    return sol, serr;

def nbodymodel(time,sol,itime=-1):
    
    nbodies=int((len(sol)-7)/7)
    npt=len(time)
    tol=1.0e-8

    #Handle integration time vars
    if type(itime) is int :
        if itime < 0 :
            itime=np.ones(npt)*0.020434
        else:
            itime=np.ones(npt)*float(itime)

            
    #Handle variables 
    ntmidmax=0
    tmin=np.min(time)
    tmax=np.max(time)
    for i in range(2,nbodies+1):
        npl=7+7*(i-1)
        per=np.copy(sol[npl+1])
        ntmidmax=np.max([ntmidmax,int( (tmax-tmin)/per )*2])

    ntmid=np.zeros(nbodies, dtype="int32")
    tmid=np.zeros(shape=(nbodies,ntmidmax), order='F')      
    percor=np.zeros(nbodies)
    ans=np.zeros(npt)
    colflag=0 #collision/close encounter flag
    itprint=0 #do not print timing measurement output
    lc.lcmodel_pc(nbodies,npt,tol,sol,time,ntmidmax,ntmid,tmid,percor,colflag,itprint)
    if colflag == 0:
        lc.percorcalc(nbodies,sol,ntmidmax,ntmid,tmid,percor)
    #print(percor)
    itprint=0 #do not print timing measurement output
    itmodel=1 #create a transit model 
    lc.lcmodel(nbodies,npt,tol,sol,time,itime,ntmidmax,ntmid,tmid,percor,ans,colflag,itprint,itmodel)
    
    return ans;

def readtt(files):
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

    npl=len(files) #number of planets to read in
    ntt =np.zeros(npl, dtype="int32") #allocate array for number of TTVs measured
    tobs=np.zeros(shape=(npl,nmax)) #allocate array for timestamps
    omc =np.zeros(shape=(npl,nmax)) #allocate array for O-C
    omcerr=np.zeros(shape=(npl,nmax)) #allocate array for O-C errors

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
                    omcerr[i,j]=float(columns[2])
            ntt[i]=j+1
    return ntt, tobs, omc, omcerr;

def genttvmodel(time,sol):

    nbodies=int((len(sol)-7)/7)
    tol=1.0e-8

    #Handle variables
    ntmidmax=0
    tmin=np.min(time)
    tmax=np.max(time)
    for i in range(2,nbodies+1):
        npl=7+7*(i-1)
        per=np.copy(sol[npl+1])
        ntmidmax=np.max([ntmidmax,int( (tmax-tmin)/per )*2])

    npt=len(time)
    ntmid=np.zeros(nbodies, dtype="int32")
    tmid=np.zeros(shape=(nbodies,ntmidmax), order='F')
    percor=np.zeros(nbodies)
    colflag=0
    itprint=0 #do not print timing measurement output
    lc.lcmodel_pc(nbodies,npt,tol,sol,time,ntmidmax,ntmid,tmid,percor,colflag,itprint)
    if colflag == 0:
        lc.percorcalc(nbodies,sol,ntmidmax,ntmid,tmid,percor) #get period correction
    itprint=1 #create a transit model
    lc.lcmodel_pc(nbodies,npt,tol,sol,time,ntmidmax,ntmid,tmid,percor,colflag,itprint)

    return ntmid,tmid;

def plotTTVandModel(sol,ntmid,tmid,ntt,tobs,omc,omcerr):

    matplotlib.rcParams.update({'font.size': 22}) #adjust font

    nbodies=int((len(sol)-7)/7)
    for i in range(1,nbodies):
        if len(tobs[i-1,:ntt[i-1]]) > 1:
            plt.figure(figsize=(12,5)) #adjust size of figure

            npl=7+7*(i)
            
            xp=tobs[i-1,:ntt[i-1]]  #predicted transit times
            yp=omc[i-1,:ntt[i-1]]   #measured OMC
            ep=1.0/omcerr[i-1,:ntt[i-1]] #weights for bestfit 
            p=np.polyfit(xp,yp,1,w=ep) #get linear fit to data
            #print(p)
            
            #t0=np.copy(sol[npl])
            #per=np.copy(sol[npl+1])
            t0=np.min(xp)
            dxp=[]
            for j in range(1,ntt[i-1]):
                dxp.append(xp[j]-xp[j-1])
            dxp=np.array(dxp)
            per=np.min(dxp)

            #calculate new period
            pernew=per*(1.0+p[0])
            #print(per,pernew,pernew-per)
            
            #Now we need to reconstruct the OMC using the polynomial fit
            xpnew=[]
            ypnew=[]
            for j in range(ntt[i-1]):
                n=np.floor((xp[j]-t0)/per+0.5)
                tc1=t0+n*pernew
                xpnew.append(tc1)
                tobs1=yp[j]+xp[j]
                omc1=tobs1-tc1
                ypnew.append(omc1)
                
                #print(n,xp[j]-tc1,tobs1-tc1)
            
            xpnew=np.array(xpnew)
            ypnew=np.array(ypnew)

            plt.errorbar(xpnew,1440*(ypnew),1440*(omcerr[i-1,:ntt[i-1]]),fmt='o',c='red')

            #plot synthetic model
            x=tmid[i,0:ntmid[i]]
            y=np.zeros(ntmid[i])
            for j in range(ntmid[i]):
                y[j]=(tmid[i,j]-t0)/pernew-np.floor((tmid[i,j]-t0)/pernew)
                if y[j] > 0.5:
                    y[j]=y[j]-1.0
            y=y*pernew #convert from phase to period
            plt.plot(x,1440*(y),c='blue')

            plt.xlabel('Time (days)')
            plt.ylabel('O-C (mins)')

            plt.show()

    return;

def intperc(x,x_eval,kde1,perc=0.6827):
    'find error bounds'
    idx = (np.abs(x_eval-x)).argmin()
    kdea=np.array(kde1(x_eval))

    n=len(x_eval)

    #print(x,idx)

    i1=1
    i2=1
    intval=0.0

    j1=idx
    j2=idx
    j1old=j1
    j2old=j2
    while intval < perc:
        j1test=np.max((0,idx-i1-1))
        j2test=np.min((n-1,idx+i2+1))
        if kdea[j1test] > kdea[j2test]:
            j1=j1test
            i1=i1+1
        elif j2old == n-1:  #case when stuck at edge.
            j1=j1test
            i1=i1+1
        else:
            j2=j2test
            i2=i2+1

        intval=np.trapz(kdea[j1:j2],x_eval[j1:j2])
        #print(j1,j2,intval)

        #make sure we can break from loop
        if j1 == 0 and j2 == n-1:  #break we reach boundaries of array
            #print('break1')
            intval=1.0
        if j1 == j1old and j2 == j2old: #break if stuck in loop.
            #print('break2')
            intval=1.0

        #Update old values to check we are making progress.
        j1old=j1
        j2old=j2

    #print(x_eval[j1],x_eval[j2])
    return x_eval[j1],x_eval[j2];

def modekdestimate(chain,burnin):
    'Estimate Mode with KDE and return KDE'
    #range of data
    minx=np.min(chain[burnin:])
    maxx=np.max(chain[burnin:])
    x_eval = np.linspace(minx, maxx, num=100)
    kde1 = stats.gaussian_kde(chain[burnin:])#,0.3)
    modeval=[]
    modekde=0
    for x in x_eval:
        if kde1(x) > modekde:
            modekde=kde1(x)
            modeval=x
    return modeval,x_eval,kde1 ;
