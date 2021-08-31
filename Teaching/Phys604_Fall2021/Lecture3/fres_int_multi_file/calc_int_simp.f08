! * *******************************************************
! * Subroutine to calculate the integral (Simpson's rule) *
! *********************************************************
subroutine calc_int_simp(nsub,x_low,x_hi,s_integral)
  implicit none

  ! Arguments
  real(8),intent(in) :: x_low,x_hi ! limits of integration
  real(8),intent(out) :: s_integral ! Final result for integral
  integer,intent(in) :: nsub ! number of subintervals
  
  
  ! Local variables
  real(8),parameter :: PI=4.0_8*atan(1.0_8) ! Pi 
  real(8),parameter :: TWO=2.0_8,SIX=6.0_8, FOUR=4.0_8 ! 8 byte prec
  real(8) :: dx ! width of subinterval
  real(8) :: xi,xmid,xf ! x points to evaluate S(x) for each subinterval
  real(8) :: s_xi,s_xmid,s_xf ! S(x) for points in each subinterval
  integer :: ii ! for loop

  ! Calculate dx
  dx = (x_hi-x_low)/real(nsub)

  ! Initialize the x points for the first subinterval
  xi = x_low
  xmid = x_low + 0.5*dx
  xf = x_low + dx

  ! Initialize the integral
  s_integral = 0.0_8

  ! Loop over subintervals
  do ii = 1,nsub        

     ! Evaluate fresnel function at xi,xmid,xf
     s_xi=dsin(PI*xi**2/TWO)
     s_xmid=dsin(PI*xmid**2/TWO)
     s_xf=dsin(PI*xf**2/TWO)

     ! Add a term to integral for this subinterval
     s_integral = s_integral + dx*(s_xi+FOUR*s_xmid+s_xf)/SIX
     
     ! Iterate points to next subinterval
     xi = xi + dx
     xmid = xmid + dx
     xf = xf + dx
           
  end do
  
  return
end subroutine calc_int_simp
