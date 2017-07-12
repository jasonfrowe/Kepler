      subroutine writef1(UNIT,ca,cb,cbp,cc,cx,cy,f1)
      implicit none
      logical debug
      integer UNIT
      complex*16 ca,cb,cc,cbp,cx,cy,f1


      debug = .false.
      if(debug) then 
            write(*,*) 
            write(*,*)  'a  =',ca
            write(*,*)  'b  =',cb
            write(*,*)  'bp =',cbp
            write(*,*)  'c  =',cc
            write(*,*)  'x  =',cx  
            write(*,*)  'y  =',cy  
            write(*,*)  'f1 =',f1  
      end if
        call writecomplex(UNIT,'a = ',ca)
        call writecomplex(UNIT,'b1= ',cb)
        call writecomplex(UNIT,'b2= ',cbp)
        call writecomplex(UNIT,'c = ',cc)
        call writecomplex(UNIT,'x = ',cx)
        call writecomplex(UNIT,'y = ',cy)
      write(UNIT,*) '---------------------'  
      call writecomplex(UNIT,' f1 = ',f1)          
      write(UNIT,*) '---------------------'
      return
      end
  
  
      subroutine writecomplex(UNIT,str,ca)
      implicit none
      character*256 text,text1,text2,vari
      character*(*) str
      character*3 sgn,i
      complex*16 ca
      real*8 ar,ai,zero
      integer UNIT
  
  
      text = ''
      text1= ''
      text2= ''
      zero = 1d-12
      sgn = ' + '
      i = ' i '
      ar  = dreal(ca)
      ai  = dimag(ca)
      vari = str//' '
  
      if(ai.lt.0d0) sgn = ' - '
      if(ar.lt.zero.and.ai.gt.0d0) sgn =''
      if(dabs(ar).lt.zero.and.dabs(ai).lt.zero) then
        write(text,*) 0.00
        write(UNIT,*) vari(:len_trim(vari))//' '//text(:len_trim(text))
        return
      end if 
    
  
      if(dabs(ar).gt.zero) then
          if(ar.gt.1d7) then
              write(text1,'(e17.4)') ar
          else
              write(text1,'(f17.8)') ar
          end if
          text = adjustl(text1)
      end if      
      if(dabs(ai).lt.zero) then
         write(UNIT,*) vari(:len_trim(vari))//' '//text(:len_trim(text))
          return
      end if
      text = text(:len_trim(text))//sgn

      if(dabs(ai).gt.1d7) then
          write(text2,'(e17.4)') dabs(ai)
      else
          write(text2,'(f17.8)') dabs(ai)
      end if  
      text2 = adjustl(text2)  
      text = text(:len_trim(text))//' '//text2(:len_trim(text2))//i
      write(UNIT,*) vari(:len_trim(vari))//' '//text(:len_trim(text))
  
      return
      end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
