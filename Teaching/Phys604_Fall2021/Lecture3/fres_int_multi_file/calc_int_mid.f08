! * *******************************************************
! * Subroutine to calculate the integral (Midpoint rule)  *
! *********************************************************
subroutine calc_int_mid(nsub,x_low,x_hi,s_integral)
  implicit none

  ! Arguments
  real(8),intent(in) :: x_low,x_hi ! limits of integration
  real(8),intent(out) :: s_integral ! Final result for integral
  integer,intent(in) :: nsub ! number of subintervals
  
  
  ! Local variables
  real(8),parameter :: PI=4.0_8*atan(1.0_8) ! Pi
  real(8),parameter :: TWO=2.0_8 ! Two, 8 byte prec
  real(8) :: dx ! width of subinterval
  real(8) :: xmid ! x point to evaluate S(x) for each subinterval
  real(8) :: s_xmid ! S(x) for points in each subinterval
  integer :: ii ! for loop

  ! Calculate dx
  dx = (x_hi-x_low)/real(nsub)

  ! Initialize the x points for the first subinterval
  xmid = x_low + 0.5*dx

  ! Initialize the integral
  s_integral = 0.0_8

  ! Loop over subintervals
  do ii = 1,nsub        

     ! Evaluate fresnel function at xi,xmid,xf
     s_xmid=dsin(PI*xmid**2/TWO)

     ! Add a term to integral for this subinterval
     s_integral = s_integral + dx*s_xmid
     
     ! Iterate points to next subinterval
     xmid = xmid + dx
           
  end do
  
  return
end subroutine calc_int_mid
