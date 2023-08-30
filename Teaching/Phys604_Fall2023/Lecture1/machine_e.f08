program mac_e
  implicit none

  real(8) :: xx = 2.0
  real(8) :: eps = 1.0

  do while(xx+eps /= xx)
     eps=eps/2.0
     
  end do

  write(*,*) 2*eps
  
end program mac_e
