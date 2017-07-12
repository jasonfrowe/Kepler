      !-----------------------------------------------------------------
      ! function  : get_transformation2 (integer)
      !
      ! package   : F1
      !
      ! Language  : Fortran 90
      !
      ! author    : F. Colavecchia (flavioc@lanl.gov)
      !                            (flavioc@cab.cnea.gov.ar)
      !
      ! date      : 3/26/97      version: 0.1
      ! revision  : 6/25/02      version: 1.0
      !
      ! purpose   : Gets the best transformation expression for
      !             a given pair of real variables to compute 
      !             Appell's F1 function following CPC 138 (2001) 29.
      !             Selection is based on the effective distance
      !             of the u,w variables (see section 4.1 and table 1).
      !
      ! input     :    x  -> real variable
      !                y  -> real variable
      !
      !
      ! output    :    integer that corresponds to the selected
      !                continuation. The get_transformation value 
      !                is the number of the equation of paper CPC 138
      !                (2001) 29
      !                to be used by the f1 code.
      !
      ! notes     :  most of the code is f77, few things come from f90.
      !
      !----------------------------------------------------------------
      integer function get_transformation2(x,y)
      implicit none
      integer flag,nopossible
      real*8 tmax1,tmax2,x,y
      logical ispossible,debug,isthere
      real u,w

!
!       Select the analytic continuation
!
!       tmax1 is slightly greater than one to
!       manage some limit cases
!       flag stores the analytic continuation
!       given by the best minimum value of tmax1
      ispossible = .false.    
      debug = .false.
      flag  = 0
      tmax1= 0.9999d0
      tmax2 = max(abs(x),abs(y))
      if(tmax2.lt.tmax1) then
            ispossible=.true.
            if(tmax2.lt.0.5d0) then
                    flag = 0
            else
                    flag = 1
            end if
            tmax1=tmax2
      end if
      if (debug) write(*,*) flag," ", tmax1,"   ",tmax2
      u = x/(x-1)
      w = y/(y-1)
          tmax2=  sqrt(u**2+w**2)
      !
      ! Eq (15)
      !
      if(tmax2.lt.tmax1) then 
                  flag=15
                  ispossible=.true.
                  tmax1=tmax2
          end if
      !
      ! Eq (16)
      !
      u = x/(x-1)
      w = (x-y)/(x-1)
          tmax2 = sqrt(u**2+w**2) 
          if(tmax2.lt.tmax1) then 
                  flag=16
                  ispossible=.true.
          tmax1=tmax2
          end if
      !
      ! Eq (17)
      !
      u = y/(y-1)
      w = (y-x)/(y-1)
          tmax2 = sqrt(u**2+w**2)  
          if(tmax2.lt.tmax1) then 
                  flag=17
                  ispossible=.true.
                  tmax1=tmax2
          end if
          !
      ! Eq (21)
      !
      u = 1-y
      w = 1-x
      tmax2= sqrt(u**2+w**2)
          if(tmax2.lt.tmax1.and.(dabs(1-y) .le. dabs(1-x))) then 
                  flag=21
                  ispossible=.true.
                  tmax1=tmax2
          end if
      !
      ! Eq (22)
      !
          if(tmax2.lt.tmax1.and.(dabs(1-x) .le. dabs(1-y))) then 
                  flag=22
                  ispossible=.true.
                  tmax1=tmax2
          end if
          if (debug) write(*,*) flag," ", tmax1,"   ",tmax2
      !
      ! Eq (23)
      !
      u = x
      w = 1/y
          tmax2= sqrt(u**2+w**2)
          if(tmax2.lt.tmax1) then 
                  flag=23
                  ispossible=.true.
                  tmax1=tmax2
          end if
      if (DEBUG) write(*,*) flag," ", tmax1,"   ",tmax2
          !
      ! Eq (24)
      !
      u = y
      w = 1/x
      tmax2= sqrt(u**2+w**2)
          if(tmax2.lt.tmax1) then 
                  flag=24
                  ispossible=.true.
                  tmax1=tmax2
          end if
      !
      ! Eq (25)
      !
      u = 1-x
      w = 1/y
          tmax2= sqrt(u**2+w**2) 
          if(tmax2.lt.tmax1) then 
                  flag=25
                  ispossible=.true.
                  tmax1=tmax2
          end if
          if (DEBUG) write(*,*) flag," ", tmax1,"   ",tmax2 
      !
      ! Eq (26)
      !
      u = 1-y
      w = 1/x
          tmax2= sqrt(u**2+w**2)          
          if(tmax2.lt.tmax1) then 
                  flag=26
                  ispossible=.true.
                  tmax1=tmax2
          end if
      if (DEBUG) write(*,*) flag," ", tmax1,"   ",tmax2 
      !
      ! Eq (27)
      !
      u = 1/y
      w = 1/x
          tmax2= sqrt(u**2+w**2)
          if(tmax2.lt.tmax1 .and. dabs(x).le.dabs(y)) then 
                  flag=27
          ispossible=.true.
            tmax1=tmax2
          end if
      !
      ! Eq (28)
      !
          if(tmax2.lt.tmax1 .and. dabs(y).le.dabs(x)) then 
                  flag=28
                  ispossible=.true.
                  tmax1=tmax2
          end if
      if (DEBUG) write(*,*) flag," ", tmax1,"   ",tmax2 
      !
      ! Eq (29)
      !
      u = (x-y)/(y*(x-1))
      w = 1/y 
          tmax2= sqrt(u**2+w**2)
          if(tmax2.lt.tmax1 .and. dabs(x-y).le.dabs(1-x)) then 
                  flag=29
                  ispossible=.true.
                  tmax1=tmax2
          end if
      !
      ! Eq (30)
      !    
      if (DEBUG)    write(*,*) flag," ", tmax1,"   ",tmax2 
      u = (x-y)/(x*(y-1))
      w = 1/x 
          tmax2= sqrt(u**2+w**2)          
          if(tmax2.lt.tmax1 .and. dabs(y-x).le.dabs(1-y)) then 
                  flag=30
                  ispossible=.true.
                  tmax1=tmax2
          end if
      if (DEBUG) write(*,*) flag," ", tmax1,"   ",tmax2 
          tmax2 = sqrt(x**2+y**2)
          if(tmax2.lt.1d0) then 
                  flag=1
                  ispossible=.true.
                  tmax1=tmax2
          end if

          if(debug) then
        inquire(file='flag.dat',exist=isthere)
        if(isthere) then
            open(unit=23,file='flag.dat',position='append')
        else
            open(unit=23,file='flag.dat',status='new')
        end if
        if(ispossible) then
                    write(23,'(2f7.3,i4)') x,y,flag
        else
                    write(23,'(2f7.3,a14)') x,y,'  Not possible' 
        end if
              close(23)
          endif
    
    
    
          if(ispossible) then
        if(debug) write(*,*) "Analytic Continuations"
      end if
      get_transformation2 = flag

      return
      end
