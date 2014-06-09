module aa_mod
implicit none
double precision :: aa = 4
end module aa_mod

subroutine calfun &
  ( &
    n , x , f &
  )

  use aa_mod

  implicit none

  integer , intent(in) &
    :: n
  double precision , intent(in) , dimension(n) &
    :: x
  double precision , intent(out) &
    :: f

  integer &
    :: ind

  f = ( x(1)-1 ) ** 2

  do ind = 2 , n
    f = f + (x(ind)-aa*x(ind-1))**2
  end do

end subroutine calfun
