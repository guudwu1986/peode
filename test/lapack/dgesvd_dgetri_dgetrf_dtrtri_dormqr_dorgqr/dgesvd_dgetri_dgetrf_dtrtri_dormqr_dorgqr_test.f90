! Test program

! Generate a random non-singular matrix via SVD "dgesvd".

! Compute its true inverse by analytical method.

! Compute inverse by "dgetri","dgetrf".

! Compute QR decomposition.
! Compute R_inv by "dtrtri".
! Compute Q_inv by "dorgqr".
! Apply Q to R_inv for inverse by "dormqr".

! Compare all inverses with true one.

program main

implicit none

integer , parameter &
  :: N = 3
integer &
  :: ind
integer &
  :: info

double precision , dimension(N,N) &
  :: a
double precision , dimension(N,N) &
  :: ai
double precision , dimension(N,N) &
  :: ai_est
double precision , dimension(N,N) &
  :: qr
double precision , dimension(N,N) &
  :: ri
integer , dimension(N) &
  :: tau

double precision , dimension(N) &
  :: s
double precision , dimension(N,N) &
  :: u
double precision , dimension(N,N) &
  :: vt

integer &
  :: lwork
double precision , dimension(:) , allocatable &
  :: work
integer , dimension(:) , allocatable &
  :: iwork

! Generate random matrix "a" and compute inverse "ai"

call random_number(a)

allocate ( work(1) )
lwork = -1
call dgesvd ( 'A' , 'A' , N , N , a , N , s , u , N , vt , N &
  , work , lwork , info )
lwork = work(1)
deallocate ( work )
allocate ( work(lwork) )
call dgesvd ( 'A' , 'A' , N , N , a , N , s , u , N , vt , N &
  , work , lwork , info )
deallocate ( work )
write(*,*) 'If not zero, generation failed:'
write(*,*) info

a = 0
ai = 0
do ind = 1 , N
  a(ind,ind) = ind
  ai(ind,ind) = 1.0 / ind
end do

a = matmul ( u , a )
a = matmul ( a , vt )
ai = matmul ( transpose(vt) , ai )
ai = matmul ( ai , transpose(u) )

! Compute inverse by "dgetri"("dgetrf")

ai_est = a

allocate ( iwork(N) )
call dgetrf ( N , N , ai_est , N , iwork , info )

allocate ( work(1) )
lwork = -1
call dgetri ( N , ai_est , N , iwork , work , lwork , info )
lwork = work(1)
deallocate ( work )
allocate ( work(lwork) )
call dgetri ( N , ai_est , N , iwork , work , lwork , info )
deallocate ( iwork )
deallocate ( work )

write(*,*) 'Error of inverse by LU:'
write(*,*) maxval ( abs ( ai_est - ai ) )

! Compute QR decomposition

qr = a

allocate ( work(1) )
lwork = -1
call dgeqrf ( N , N , qr , N , tau , work , lwork , info )
lwork = work(1)
deallocate ( work )
allocate ( work(lwork) )
call dgeqrf ( N , N , qr , N , tau , work , lwork , info )
deallocate ( work )

! Compute R_inv by "dtrtri"

ri = 0
do ind = 1 , N
  ri ( 1 : ind , ind ) = qr ( 1 : ind , ind )
end do

call dtrtri ( 'U' , 'N' , N , ri , N , info )

! Compute Q by "dorgqr"

ai_est = qr

allocate ( work(1) )
lwork = -1
call dorgqr ( N , N , N , ai_est , N , tau , work , lwork , info )
lwork = work(1)
deallocate ( work )
allocate ( work(lwork) )
call dorgqr ( N , N , N , ai_est , N , tau , work , lwork , info )
deallocate ( work )

ai_est = matmul ( ri , transpose(ai_est) )
write(*,*) 'Error of inverse by normal QR:'
write(*,*) maxval ( abs ( ai_est - ai ) )

! Directly apply Q to R by "dormqr"

ai_est = ri

allocate ( work(1) )
lwork = -1
call dormqr ( 'R' , 'T' , N , N , N , qr , N , tau , ai_est , N &
  , work , lwork , info )
lwork = work(1)
deallocate ( work )
allocate ( work(lwork) )
call dormqr ( 'R' , 'T' , N , N , N , qr , N , tau , ai_est , N &
  , work , lwork , info )
deallocate ( work )

write(*,*) 'Error of inverse by directly apply Q to R:'
write(*,*) maxval ( abs ( ai_est - ai ) )

stop

end program main
