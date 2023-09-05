program main
  use compute_module

  implicit none
  
  real :: xx,yy,zz

  xx=1.0
  yy=1.0
  zz=compute(xx,yy) ! Call function

  call print_result(zz) ! Call subroutine
  
end program main
