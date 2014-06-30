! Testing program for subroutine "dgeqrf"

! Generate matrix A randomly,
! do decomposition,
! then compare A with QR.

program main

  implicit none

  integer , parameter &
    :: m = 100
  integer , parameter &
    :: n = 100

  double precision , dimension(m*n) , target &
    :: a
  double precision , dimension(m*n) , target &
    :: o
  double precision , dimension(m*n) , target &
    :: r
  double precision , dimension(m) &
    :: v
  double precision , dimension(m,m) &
    :: temp

  double precision , dimension(:,:) , pointer &
    :: p_a
  double precision , dimension(:,:) , pointer &
    :: p_r
  double precision , dimension(:,:) , pointer &
    :: p_o

  double precision , dimension(min(m,n)) &
    :: tau
  double precision , dimension(:) , allocatable &
    :: work
  integer &
    :: lwork
  integer &
    :: info

  integer &
    :: index1
  integer &
    :: index2

  call random_number(a)
  o = a

  p_a ( 1:m , 1:n ) => a
  p_r ( 1:m , 1:n ) => r
  p_o ( 1:m , 1:n ) => o

  allocate ( work(1) )
  lwork = -1
  call dgeqrf ( m , n , a , m , tau , work , lwork , info )
  lwork = work(1)
  deallocate ( work )
  allocate ( work(lwork) )
  call dgeqrf ( m , n , a , m , tau , work , lwork , info )
  deallocate ( work )

  write(*,*) info

  r = 0
  do index1 = 1 , m
    do index2 = index1 , n
      p_r ( index1 , index2 ) = p_a ( index1 , index2 )
    end do
  end do

  do index1 = min(m,n) , 1 , -1
    v = 0
    v ( index1 ) = 1
    v ( (index1+1) : m ) = p_a ( (index1+1) : m , index1 )
    temp = matmul ( reshape(v,(/m,1/)) , reshape(v,(/1,m/)) )
    p_r = p_r - tau(index1) * matmul ( temp , p_r )
  end do

  write(*,*) maxval ( abs ( p_r - p_o ) )

  stop

end program main
