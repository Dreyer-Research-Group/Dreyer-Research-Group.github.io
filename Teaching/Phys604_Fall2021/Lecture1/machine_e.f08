program mac_e
  implicit none

  real(8) :: xx = 2.0_8
  real(8) :: eps = 1.0_8

  do while(xx+eps /= xx)
     eps=eps/2.0_8
     
  end do

  write(*,*) 2*eps
  
end program mac_e
