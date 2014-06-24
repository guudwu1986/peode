! Testing program for subroutine "dgeqrf"

! Generate matrix A randomly,
! do decomposition,
! then compare A^T A - R^T R.

program main

  implicit none

  integer , parameter &
    :: m = 100
  integer , parameter &
    :: n = 100

  double precision , dimension(m*n) , target &
    :: a
  double precision , dimension(m,m) &
    :: ata

  double precision , dimension(:,:) , pointer &
    :: p

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

  p ( 1:m , 1:n ) => a

  ata = matmul ( transpose(p) , p )

  if ( allocated(work) ) then
    stop
  end if

  allocate ( work(1) )
  lwork = -1
  call dgeqrf ( m , n , a , m , tau , work , lwork , info )
  lwork = work(1)
  deallocate ( work )
  allocate ( work(lwork) )
  call dgeqrf ( m , n , a , m , tau , work , lwork , info )
  deallocate ( work )

  write(*,*) info

  do index1 = 2 , m
    do index2 = 1 , min(index1-1,n)
      p ( index1 , index2 ) = 0
    end do
  end do

  ata = ata - matmul ( transpose(p) , p )
  write(*,*) maxval(abs(ata))

  stop

end program main
