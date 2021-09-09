  ! Purpose: Program to use Euler's method for a body rotating around
  ! the sun
  ! Author: Cyrus Dreyer
  ! Date: 4/18/21
program euler_orbit
  implicit none

  ! Variable dictionary 
  integer :: n_time,ii ! number of time steps and loop variable
  integer :: lun_out ! for gnuplot file
  real :: x0,y0,vx0,vy0 ! Initial position in AU and velocity in AU/year
  real :: time_f ! final time in years
  real :: delta_t ! time step in years
  real :: time ! loop variable for each time step in years
  real :: dvdt ! function to calculate acceleration
  real, allocatable :: pos_x(:),pos_y(:), vel_x(:), vel_y(:) ! Position and velocity
  

  
  ! Set the initial conditions
  x0=1.0
  y0=0.0
  vx0=0.0
  vy0=6.283185
  time=0.0
  time_f=2.0
  delta_t=1.0/365.0

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

     ! Increment position and velocity in time
     pos_x(ii+1)=pos_x(ii) + vel_x(ii)*delta_t
     pos_y(ii+1)=pos_y(ii) + vel_y(ii)*delta_t
     vel_x(ii+1)=vel_x(ii) + dvdt(pos_x(ii),pos_y(ii))*delta_t
     vel_y(ii+1)=vel_y(ii) + dvdt(pos_y(ii),pos_x(ii))*delta_t
     
  end do

  ! Make gnuplot file
  open(newunit=lun_out,file='euler_orbit_plot.dat',status='replace')

  write(lun_out,*) "set terminal png size 640,640"
  write(lun_out,*) "set output 'euler_orbit_plot.png'"
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

end program euler_orbit


! Function for calculating acceleration
real function dvdt(x_var,y_var)
  implicit none

  ! Arguements
  real, intent(in) :: x_var,y_var ! x and y position in AU
  
  ! Local variables
  real, parameter :: GRAV_G = 39.47 ! gravitational constant in AU*yr^2/Msun
  real, parameter :: M_SUN = 1 ! Mass of the sun, on if G is in the units above

  dvdt = -GRAV_G*M_SUN*x_var/(x_var**2 + y_var**2)**(3.0/2.0)
  
  return
end function dvdt
