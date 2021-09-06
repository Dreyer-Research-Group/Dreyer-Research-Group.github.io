  ! Purpose: find the root of x^3 + 6 with the Newton-Raphson method
  ! Author: Cyrus Dreyer
  ! Date: 02/28/2019
program NR_root
  implicit none
  real(16) :: epsilon,x_0 ! Input tolerance and guess root
  real(16) :: x_n, x_n_m_1 ! Value for current and previous guess
  integer,parameter :: NMAX = 1000 ! Max iterations before we admit defeat
  integer :: nn ! loop variable
 

  ! Prompt the user and read in initial guess and tolerance
  write(*,*) "Enter an inital guess for the root of (x^3 + 6):"
  read(*,*) x_0
  write(*,*) "Enter a tolerance for the root finding:"
  read(*,*) epsilon


  ! Loop over NR steps
  x_n_m_1=x_0
  do nn=1,NMAX

     ! make sure we do not divide by zero. Derivative is 3*x^2
     if ( abs((3.0*(x_n_m_1)**2)) < 1.0d-10 ) then
        write(*,*) " A derivative of zero has been obtained,", & 
             & " root finding failed, try a different guess"
        stop 1
     end if

     ! Calculate new guess
     x_n=x_n_m_1-( ((x_n_m_1)**3 + 6.0) / (3.0*(x_n_m_1)**2) )

     write(*,'(a30,i6,a5,f20.10)') 'New guess at step:',nn,'is:',x_n
     
     ! Check for convergence, write out root and stop program if converged
     if (abs(x_n) < 1.0d-20 .and. abs(x_n-x_n_m_1) < epsilon) then
        write(*,*) "Root of (x^3 + 6) is at x=",x_n
        stop 0
        
     else if ( abs(x_n-x_n_m_1) < abs(epsilon*x_n)) then
        write(*,*) "Root of (x^3 + 6) is at x=",x_n
        stop 0

     end if

     ! Set x_(n-1) = x_n for next round of the loop
     x_n_m_1=x_n
     
  end do ! nn loop

  ! If we have exceeded max intervals, warn user and stop program
  write(*,*) "After",NMAX,"steps, the root could not converge below",epsilon
  stop 2
  

end program NR_root
