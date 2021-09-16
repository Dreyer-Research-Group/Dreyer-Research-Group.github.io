  ! Purpose: Program to use Bulirsch-Stoer method for solving the
  ! orbital problem
  ! Author: Cyrus Dreyer
  ! Date: 09/13/21
program Bul_Stoer_orbit
  implicit none

  ! Variable dictionary 
  integer :: n_time,ii,nn,mm,var ! number of time steps and loop variables
  integer :: lun_out,lun_log ! for gnuplot file and log file
  integer,parameter :: N_RI = 8 ! Max number of Richardson Extrapolation steps
  integer :: ri_conv_step ! Test if the RI steps converged within N_RI
  integer :: tot_timesteps ! Accumulate total number of timesteps
  real(16) :: x0,y0,vx0,vy0 ! Initial position in AU and velocity in AU/year
  real(16) :: time_f ! final time in years
  real(16) :: delta_t,delta_t_bs ! time step in years, and BS timestep
  real(16) :: time ! loop variable for each time step in years
  real(16) :: pos_x_bs,pos_y_bs,vel_x_bs,vel_y_bs ! position and velocity calculated by BS method
  real(16) :: des_error, error ! desired error and actual error for timestep
  real(16) :: Rij(4,N_RI,N_RI)=0.0d0 ! Stores Richardson extrapolation variables; we will limit to N_RI steps
  real(16) :: Rij_err(4,N_RI,N_RI)=0.0d0 ! Stores Richardson extrapolation errors; we will limit to N_RI steps
  real(16), allocatable :: pos_x(:),pos_y(:), vel_x(:), vel_y(:) ! Position and velocity for trajectory 
  real(16), allocatable :: step_size(:) ! Step size for the BS step to converge
  
  ! open a file to serve as a log file
  open(newunit=lun_log,file='log.dat',status='replace')

  ! Set the initial conditions
  x0=0.3
  y0=0.0
  vx0=0.0
  vy0=14.955378073901486 
  time=0.0
  time_f=2.0
  delta_t=1.0/36.50 ! This is now the initial delta_t

  des_error = 1.0e-6 ! Desired accuracy per step

  ! Put this information in a log file
  write(lun_log,*) '# Initial conditions (position in AU, velocity in AU/year:'
  write(lun_log,'(4(a5,e12.4))') 'x0:',x0,'y0:',y0,'vx0:',vx0,'vy0:',vy0
  write(lun_log,'(2(a20,e12.4))') 'Initial time:',time,', final time:',time_f
  write(lun_log,'(a20,e12.4)') 'Initial timestep:',delta_t
  write(lun_log,'(a30,e12.4)') 'Desired error/step:',des_error
  write(lun_log,*)
  
  
  ! Estimate how many time steps we will have. Note that this may not
  ! be acccurate if our adaptive algorithm changes delta_t. We will
  ! deal with this later. Also, divide by two since each step we use
  ! will increment by 2*delta_t
  n_time = ceiling((time_f-time)/(delta_t))
  

  ! Allocate and initialize arrays for position and velocity 
  allocate(pos_x(n_time))
  allocate(pos_y(n_time))
  allocate(vel_x(n_time))
  allocate(vel_y(n_time))
  allocate(step_size(n_time))

  pos_x(1) = x0
  pos_y(1) = y0
  vel_x(1) = vx0
  vel_y(1) = vy0
  step_size(1) = 0
  tot_timesteps = 0
  
  ! We will perform the BS method for each time step
  do ii=1,n_time-1

     ! Initial position and velocity for BS method
     pos_x_bs=pos_x(ii)
     pos_y_bs=pos_y(ii)
     vel_x_bs=vel_x(ii)
     vel_y_bs=vel_y(ii)

     ! Initial application of MMM step
     call mod_midpoint(delta_t,1,pos_x_bs,pos_y_bs,vel_x_bs,vel_y_bs)
     Rij(1,1,1)=pos_x_bs
     Rij(2,1,1)=pos_y_bs
     Rij(3,1,1)=vel_x_bs
     Rij(4,1,1)=vel_y_bs
     
     ! Loop over modified midpoint method intervals.
     ri_conv_step = -1
     do nn = 2,N_RI

        ! Divide up the inteval into nn points
        delta_t_bs=delta_t/nn

        ! Reinitial position and velocity for BS method
        pos_x_bs=pos_x(ii)
        pos_y_bs=pos_y(ii)
        vel_x_bs=vel_x(ii)
        vel_y_bs=vel_y(ii)

        
        ! Modified midpoint step to get R_{n,1}
        call mod_midpoint(delta_t_bs,nn,pos_x_bs,pos_y_bs,vel_x_bs,vel_y_bs)

        Rij(1,nn,1)=pos_x_bs
        Rij(2,nn,1)=pos_y_bs
        Rij(3,nn,1)=vel_x_bs
        Rij(4,nn,1)=vel_y_bs

        ! Loop for the Richardson Extrapolation
        do mm=1,nn-1

           ! Loop through variables
           do var = 1,4
              call richardson_extrap(Rij(var,nn,mm+1),Rij(var,nn,mm),Rij(var,nn-1,mm), &
                   & nn,mm,Rij_err(var,nn,mm))  
           end do
           
        end do ! mm
        
        ! Test if we are below desired error. We care about Euclidian
        ! norm of position.
        error= sqrt(Rij_err(1,nn,nn-1)**2+Rij_err(2,nn,nn-1)**2)
        if (error < des_error) then
           ri_conv_step = nn
           tot_timesteps = tot_timesteps + nn 
           exit
        end if
        
     end do ! nn Loop over MMM intervals

     ! If we have not converged RI
     if (ri_conv_step < 0) then
        write(*,*) 'Richardson Extrapolation not converged within',N_RI,&
             &'steps. Either increase RI steps or decrease timestep.'
        write(lun_log,*) 'Richardson Extrapolation not converged within',N_RI,&
             &'steps. Either increase RI steps or decrease timestep.'

        stop 1
     end if
     
     ! Store converged position and velocity
     pos_x(ii+1) = Rij(1,ri_conv_step,ri_conv_step)
     pos_y(ii+1) = Rij(2,ri_conv_step,ri_conv_step)
     vel_x(ii+1) = Rij(3,ri_conv_step,ri_conv_step)
     vel_y(ii+1) = Rij(4,ri_conv_step,ri_conv_step)
     step_size(ii+1) = delta_t/ri_conv_step

     ! Output Rij and error to log file     
     write(lun_log,*)
     write(lun_log,*) '# Richarson extrapolation for interval',ii
     do var = 1,4
        write(lun_log,*)
        if (var==1) write(lun_log,*) '# x position:'
        if (var==2) write(lun_log,*) '# y position:'
        if (var==3) write(lun_log,*) '# x velocity:'
        if (var==4) write(lun_log,*) '# y velocity:'

        do nn=1,ri_conv_step
           do mm=1,ri_conv_step

              if (mm==ri_conv_step) then
                 write(lun_log,'(e20.10)') Rij(var,nn,mm)
              else
                 write(lun_log,'(e20.10)',advance="no") Rij(var,nn,mm)
              end if
           end do
        end do
     end do
     
     ! Iterate to next time interval
     time=time+delta_t
     time_f=time_f+delta_t
     
  end do ! Time intervals to n_time
     

  write(lun_log,*)
  write(lun_log,*) 'Calculation complete after',tot_timesteps,'timesteps'
  write(*,*) 'Calculation complete after',tot_timesteps,'timesteps'
  
  ! Make output data
  open(newunit=lun_out,file='BS_orbit_plot.dat',status='replace')

  write(lun_out,*) '# x position, y position, x velocity, y velocity'
  
  do ii=1,n_time
     write(lun_out,'(5e30.20)') pos_x(ii),pos_y(ii),vel_x(ii),vel_y(ii),step_size(ii)
  end do


  ! Deallocate arrays, close files
  deallocate(pos_x)
  deallocate(pos_y)
  deallocate(vel_x)
  deallocate(vel_y)

  close(lun_out)
  close(lun_log)
  
  stop 0
end program Bul_Stoer_orbit

