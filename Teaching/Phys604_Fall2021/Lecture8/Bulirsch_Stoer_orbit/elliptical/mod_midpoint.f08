! *************************************************
! * Subroutine for Modified Midpoint Method Steps *
! *************************************************
subroutine mod_midpoint(delta_t,n_time,pos_x,pos_y,vel_x,vel_y)
  implicit none

  ! Arguments
  integer,intent(in) :: n_time ! number of timesteps
  real(16), intent(in) :: delta_t ! Timestep
  real(16), intent(inout) :: pos_x,pos_y,vel_x,vel_y ! Positions and velocities (to be used at integer pts)

  ! Local valriables
  real(16),parameter :: HALF=0.5_16 ! 16 byte version of 0.5
  integer :: tt ! For timestep loop
  real(16) :: x_mid,y_mid,vx_mid,vy_mid ! position and velocity at midpoint
  real(16) :: dvdt ! function to calculate acceleration
  real(16) :: k1_vel_x,k1_vel_y ! Temp variables to hold accelerations
  

  ! Initial Euler half step for position and velocity
  x_mid=pos_x+HALF*vel_x*delta_t
  y_mid=pos_y+HALF*vel_y*delta_t     
  vx_mid=vel_x+HALF*dvdt(pos_x,pos_y)*delta_t
  vy_mid=vel_y+HALF*dvdt(pos_y,pos_x)*delta_t

  
  ! Loop over timestep
  do tt=1,n_time

     ! Position at t+\Delta t
     pos_x=pos_x+delta_t*vx_mid
     pos_y=pos_y+delta_t*vy_mid

     ! Velocity at t+\Delta t
     k1_vel_x=dvdt(pos_x,pos_y)*delta_t
     k1_vel_y=dvdt(pos_y,pos_x)*delta_t
     vel_x=vx_mid+HALF*k1_vel_x
     vel_y=vy_mid+HALF*k1_vel_y

     if (tt<n_time) then
     ! Velocity at t+(3/2)\Delta t
     vx_mid=vx_mid+k1_vel_x
     vy_mid=vy_mid+k1_vel_y

     ! Position at t+(3/2)\Delta t
     x_mid=x_mid+vel_x*delta_t
     y_mid=y_mid+vel_y*delta_t
     end if
  end do

  ! Now we average our two results for the position/velocity at the
  ! endpoint. This is our final result for the MMM step
  pos_x=HALF*(pos_x+x_mid+HALF*delta_t*vel_x)
  pos_y=HALF*(pos_y+y_mid+HALF*delta_t*vel_y)
  vel_x=HALF*(vel_x+vx_mid+HALF*delta_t*dvdt(pos_x,pos_y))
  vel_y=HALF*(vel_y+vy_mid+HALF*delta_t*dvdt(pos_y,pos_x))
  
  
  return
  
end subroutine mod_midpoint

! *****************************************
! * Function for calculating acceleration *
! *****************************************
real(16) function dvdt(x_var,y_var)
  implicit none

  ! Arguments
  real(16), intent(in) :: x_var,y_var ! x and y position in AU
  
  ! Local variables
  real(16), parameter :: GRAV_G = 39.47 ! gravitational constant in AU*yr^2/Msun
  real(16), parameter :: M_SUN = 1 ! Mass of the sun, on if G is in the units above

  dvdt = -GRAV_G*M_SUN*x_var/(x_var**2 + y_var**2)**(3.0/2.0)
  
  return
end function dvdt

! *******************************************
! * Subroutine for Richardson Extrapolation *
! *******************************************
subroutine richardson_extrap(R_n_mp1,R_n_m,R_nm1_m,nn,mm,error)
  implicit none

  ! Arguments
  integer :: nn,mm ! n and m indicies
  real(16),intent(out) :: R_n_mp1 ! Extrapolated R_{n,m+1}
  real(16),intent(out) :: error ! Estimated error
  real(16),intent(in) :: R_n_m,R_nm1_m ! Inputted R_{n,m} and R_{n-1,m}} 

  error = (R_n_m-R_nm1_m)/( ( real(nn)/(real(nn)-1) )**(2*mm) - 1)
  R_n_mp1 = R_n_m + error
  
  return
end subroutine richardson_extrap
