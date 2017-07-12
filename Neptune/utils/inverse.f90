function inv(A) result(Ainv)
  use precision
  implicit none
  real(double), dimension(:,:), intent(in) :: A
  real(double), dimension(size(A,1),size(A,2)) :: Ainv

  real(double), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

 !  write(0,*) "Gettings started.. "
  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)
  info=0

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  !write(0,*) "DGETRF",n
  call DGETRF(n, n, Ainv, n, ipiv, info)
 ! write(0,*) "done.. "

  if (info /= 0) then
     write(0,*) "Info: ",info
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv
