  ! Purpose: Program to use 4th order RK method for a body rotating around
  ! the sun
  ! Author: Cyrus Dreyer
  ! Date: 4/18/21
program rk4_orbit
  implicit none

  ! Variable dictionary 
  integer :: n_time,ii ! number of time steps and loop variable
  integer :: lun_out ! for gnuplot file
  real(16) :: x0,y0,vx0,vy0 ! Initial position in AU and velocity in AU/year
  real(16) :: time_f ! final time in years
  real(16) :: delta_t ! time step in years
  real(16) :: time ! loop variable for each time step in years
  real(16) :: dvdt ! function to calculate acceleration
  real(16), allocatable :: pos_x(:),pos_y(:), vel_x(:), vel_y(:) ! Position and velocity

  
  ! For esitmates of changes over time step RK2
  real(16) :: k1_pos_x,k1_pos_y,k1_vel_x,k1_vel_y
  real(16) :: k2_pos_x,k2_pos_y,k2_vel_x,k2_vel_y
  real(16) :: k3_pos_x,k3_pos_y,k3_vel_x,k3_vel_y
  real(16) :: k4_pos_x,k4_pos_y,k4_vel_x,k4_vel_y

  real(16) :: x_mid,y_mid,vx_mid,vy_mid
  
  ! Set the initial conditions
  x0=1.0
  y0=0.0
  vx0=0.0
  vy0=6.28318530718
  time=0.0
  time_f=2.0
  delta_t=1.0/36.50

  ! Determine how many time steps we will have
  n_time = ceiling(time_f/delta_t)
  

  ! Allocate and initialize arrays for position and velocity 
  allocate(pos_x(n_time))
  allocate(pos_y(n_time))
  allocate(vel_x(n_time))
  allocate(vel_y(n_time))

  pos_x(1) = x0
  pos_y(1) = y0
  vel_x(1) = vx0
  vel_y(1) = vy0
  
  ! Loop to integrate differential equations
  do ii=1,n_time-1

     ! Estimate change over time step
     k1_pos_x=vel_x(ii)*delta_t
     k1_pos_y=vel_y(ii)*delta_t
     k1_vel_x=dvdt(pos_x(ii),pos_y(ii))*delta_t
     k1_vel_y=dvdt(pos_y(ii),pos_x(ii))*delta_t

     x_mid=pos_x(ii)+0.5*k1_pos_x
     y_mid=pos_y(ii)+0.5*k1_pos_y
     vx_mid=vel_x(ii)+0.5*k1_vel_x
     vy_mid=vel_y(ii)+0.5*k1_vel_y

     k2_pos_x= vx_mid*delta_t
     k2_pos_y= vy_mid*delta_t
     k2_vel_x= dvdt(x_mid,y_mid)*delta_t
     k2_vel_y= dvdt(y_mid,x_mid)*delta_t

     x_mid=pos_x(ii)+0.5*k2_pos_x
     y_mid=pos_y(ii)+0.5*k2_pos_y
     vx_mid=vel_x(ii)+0.5*k2_vel_x
     vy_mid=vel_y(ii)+0.5*k2_vel_y

     k3_pos_x= vx_mid*delta_t
     k3_pos_y= vy_mid*delta_t
     k3_vel_x= dvdt(x_mid,y_mid)*delta_t
     k3_vel_y= dvdt(y_mid,x_mid)*delta_t

     x_mid=pos_x(ii)+k3_pos_x
     y_mid=pos_y(ii)+k3_pos_y
     vx_mid=vel_x(ii)+k3_vel_x
     vy_mid=vel_y(ii)+k3_vel_y

     k4_pos_x= vx_mid*delta_t
     k4_pos_y= vy_mid*delta_t
     k4_vel_x= dvdt(x_mid,y_mid)*delta_t
     k4_vel_y= dvdt(y_mid,x_mid)*delta_t
     
     
     
     ! Increment position and velocity in time
     pos_x(ii+1)=pos_x(ii) + (1.0/6.0)*(k1_pos_x+2.0*k2_pos_x+2.0*k3_pos_x+k4_pos_x)
     pos_y(ii+1)=pos_y(ii) + (1.0/6.0)*(k1_pos_y+2.0*k2_pos_y+2.0*k3_pos_y+k4_pos_y)
     vel_x(ii+1)=vel_x(ii) + (1.0/6.0)*(k1_vel_x+2.0*k2_vel_x+2.0*k3_vel_x+k4_vel_x)
     vel_y(ii+1)=vel_y(ii) + (1.0/6.0)*(k1_vel_y+2.0*k2_vel_y+2.0*k3_vel_y+k4_vel_y)
     
  end do

  ! Make gnuplot file
  open(newunit=lun_out,file='rk4_orbit_plot.dat',status='replace')

  write(lun_out,*) "set terminal png size 640,640"
  write(lun_out,*) "set output 'rk4_orbit_plot.png'"
  write(lun_out,*) "plot '-' with lines"

  do ii=1,n_time
     write(lun_out,*) pos_x(ii),pos_y(ii)
  end do


  ! Deallocate arrays
  deallocate(pos_x)
  deallocate(pos_y)
  deallocate(vel_x)
  deallocate(vel_y)

    
  stop 0
end program rk4_orbit


! Function for calculating acceleration
real(16) function dvdt(x_var,y_var)
  implicit none

  ! Arguements
  real(16), intent(in) :: x_var,y_var ! x and y position in AU
  
  ! Local variables
  real(16), parameter :: GRAV_G = 39.47 ! gravitational constant in AU*yr^2/Msun
  real(16), parameter :: M_SUN = 1 ! Mass of the sun, on if G is in the units above

  dvdt = -GRAV_G*M_SUN*x_var/(x_var**2 + y_var**2)**(3.0/2.0)
  
  return
end function dvdt
