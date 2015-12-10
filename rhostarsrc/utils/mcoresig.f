      double precision function mcoresig(tage,tmcore,nmodel,age,mcore)
      implicit none
      integer nmodel,i
      double precision tage(nmodel),tmcore(nmodel),age,mcore,mavg,
     .  ageavg,tcoresum,modelsum,perc(8),chi(8),percore,chisq
      logical update
      data perc/0.0,0.683,0.90,0.954,0.99,0.9973,0.9999,1.0/
      data chi/0.0,1.0,2.71,4.00,6.63,9.00,15.1,99.9/
      
c      write(0,*) age,mcore
      
      update=.true.
      tcoresum=0.0d0
      do 10 i=1,nmodel-1
        mavg=1.0d0-(tmcore(i)+tmcore(i+1))/2.0d0
        ageavg=(tage(i)+tage(i+1))/2.0d0
        tcoresum=tcoresum+mavg*(tage(i+1)-tage(i-1))
        if((age.le.ageavg).and.(update))then
c            write(0,*) "hello"
            modelsum=tcoresum
            update=.false.
        endif
c        write(0,*) ageavg,mavg,tcoresum
 10   continue
c      write(0,*) "tcoresum",tcoresum,modelsum
c      write(0,*) modelsum/tcoresum
      percore=modelsum/tcoresum
      if(age.gt.(tage(nmodel)+tage(nmodel-1))/2.0d0) percore=1.0d0
      
      
      update=.true.
      do 11 i=2,8
        if((percore.le.perc(i)).and.(update))then
            update=.false.
            chisq=chi(i-1)+(chi(i)-chi(i-1))*(percore-perc(i-1))/
     .          (perc(i)-perc(i-1))
c            write(0,*) i,chi(i-1)
        endif
 11   continue
c      write(0,*) "X^2",chisq,sqrt(chisq)
      mcoresig=sqrt(chisq)
      
c      read(5,*)
      
      return
      end