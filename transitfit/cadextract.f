      program cadextract
C     08/25/2009 - Jason Rowe
C     This program extracts cadences in and out of transit for 
C     generation of difference images.
      implicit none
      integer nunit,nfit,iargc,ncad1,ncad2,ncadmax,i,ncad3,ncad4,
     .  ncadend
      parameter(nfit=16,ncadmax=2000)
      integer nincad(ncadmax),nin,noutcad(ncadmax),nout,id1,id2,kepid
      double precision sol(nfit),serr(nfit,2),Dpvary(nfit),err(nfit),
     .  doe,toff,stime,dt,epo1,epo,transitdur,sec2day,tdur,epoN,dstart,
     .  tdurd2,etime
      character*80 sfile,pfile,dumc,cline
      
C     Needed parameters
      stime=54964.0109938200 !start of Q1 data.
      etime=54997.4812249200 !end of Q1 data
      dt=0.0204337 !30 minute cadence 
      sec2day=86400.0d0 !number of seconds in a day
      ncadend=int((etime-stime)/dt) !largest cadence
       
C     Make sure we have read in at least two commandline arguments
      if(iargc().lt.2) goto 901
      
C     Get file name for transit solution and open file
      call getarg(1,sfile)
      nunit=10
      open(unit=nunit,file=sfile,status='old',err=902)
C     Read in model solution
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
      close(nunit)
      
      call getarg(2,pfile)
      open(unit=nunit,file=pfile,status='old',err=903)
      do 10 i=1,3
        read(nunit,501) dumc
 501    format(A1)
 10   continue
      read(nunit,503) cline
 503  format(A80)
      read(cline,*,err=11,end=11) id1,id2
      goto 12
 11     read(cline,*) id1
        id2=1
 12   continue
      read(nunit,*) kepid
      
      close(nunit)
      
C     Convert HJD to MJD [sol(7) is the model time of first transit]
      epo=(sol(7)+53.0*4.1d-5)/(1.0d0+4.1d-5)-0.50020+54900.0
C     Get Epoch of first transit in Q1 [sol(5) is the period in days]
      dstart=epo-stime !diff between transit model time and Q1 start
      epo1=epo-sol(5)*floor(dstart/sol(5)) !first visible transit time
      tdur=transitdur(nfit,sol)/sec2day !get transit duration [days]
      tdurd2=tdur/2.0d0
      
      epoN=epo1
      nin=0
      nout=0
      do while(epoN.le.etime)
        ncad1=int((epoN-tdurd2-stime)/dt) !cadence of transit start
        if(ncad1.lt.1) ncad1=1
        ncad2=int((epoN+tdurd2-stime)/dt) !candence of transit stop
        if(ncad2.gt.ncadend) ncad2=ncadend
c        write(6,*) epoN,ncad1,ncad2
        nin=nin+1
        nincad(nin)=ncad1 !store start
        nin=nin+1
        nincad(nin)=ncad2 !store end
C       get cadences out of transit, but nearby 
        ncad3=int((epoN-tdur-stime)/dt) !cadence of out-of-transit start
        ncad4=int((epoN+tdur-stime)/dt) !cadence of out-of-transit end
        if(ncad3.lt.1) ncad3=1
        if(ncad4.gt.ncadend) ncad4=ncadend
        if(ncad1.gt.1)then !handle boundaries
            nout=nout+1
            noutcad(nout)=ncad3
            nout=nout+1
            noutcad(nout)=ncad1-1
        endif
        if(ncad2.lt.ncadend)then !handle boundaries
            nout=nout+1
            noutcad(nout)=ncad2+1
            nout=nout+1
            noutcad(nout)=ncad4
        endif
C       Move to next transit
        epoN=epoN+sol(5) !update epoN
      enddo
      
      write(6,502) id1,".0",id2,kepid,nin/2,nout/2
      write(6,500) (nincad(i),nincad(i+1),i=1,nin,2)
      write(6,500) (noutcad(i),noutcad(i+1),i=1,nout,2)
 502  format(I3,A2,I1,1X,I8,1X,I4,1X,I4)
 500  format(500(I4,1X,I4,1X))     
      
      goto 999
 901  write(0,*) "Usage: cadextract <sfile> <plotfile>"
      write(0,*) "Where <sfile> is the transit solution data file"
      write(0,*) "      <plotfile> is the plotting file with KOI,KepID"
      goto 999
 902  write(0,*) "Cannot open: ",sfile
      goto 999
 903  write(0,*) "Cannot open: ",pfile
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This is for circular orbits only.
      implicit none
      integer nfit
      double precision sol(nfit),Psec,M1,M2,R1,R2,asemi,temp(4),incl
C     Read in physical constants (Pi,G,Msun,..etc.)      
      include "utils/physcons.f"
      
      M1=sol(1)*Msun !kg ; mass of star
      M2=sol(2)*Mjup !kg ; mass of planet
      Psec=sol(5)*24.0*60.0*60.0 !sec ; period of planet
      asemi=(Psec*Psec*G*(M1+M2)/(4.0*Pi*Pi))**(1.0/3.0) !m
      R1=sol(3)*Rsun  !radius of star
      R2=sol(4)*Rjup  !radius of planet
      incl=Pi*sol(6)/180.0d0
      
      temp(1)=Psec/Pi
      temp(2)=R1/asemi
      temp(3)=(1+(R2/R1))**2.0-((asemi/R1)*cos(incl))**2.0
      temp(4)=1-cos(incl)*cos(incl)
      
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))
      
      return
      end