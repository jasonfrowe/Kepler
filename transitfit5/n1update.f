      program n1update
C     update n1.dat files for change from b^2 to b.
      implicit none
      integer nfit,nunit,nplanet,nplanetmax,i
      parameter(nfit=108,nplanetmax=10)
      double precision sol(nfit),serr(nfit,2),err(nfit)
      character*3 titles(15)
      character*80 inputsol


      call getarg(1,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=901)
C     We start by reading in solution from input file
      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
      write(0,*) "done reading input solution"
      close(nunit) !release unit number as we are done with file

C     change b^2 to b..
      do 25 i=1,nplanet
        sol(11+10*(i-1))=sqrt(sol(11+10*(i-1)))
 25   continue


      call exportfit(nfit,nplanet,sol,serr,err,titles)

      goto 999
 901  write(0,*) "Error opening ",inputsol
      goto 999
 999  end
