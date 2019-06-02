module bfittingmod
	use precision, only: double
	integer,pointer :: nbodies2,nbc2
	real(double), pointer :: t2
	!import mercury vars
	integer, pointer :: algor2,nbig2,ngflag2,opflag2,colflag2,nbod2
	integer, dimension(:), pointer :: opt2
	integer, dimension(:), pointer :: stat2
	real(double), pointer :: rcen2,rmax2,tstart2,tol2
	real(double), dimension(:), pointer :: jcen2,en2,am2
	real(double), dimension(:), pointer :: rphys2,rce2,rcrit2,m2
	real(double), dimension(:,:), pointer :: s2,x2,v2
	!import mercury save vars
	real(double), pointer :: a2(:,:),hrec2,angf2(:,:),ausr2(:,:)
	!current solution
	real(double), dimension(:), pointer :: sol2
end module bfittingmod