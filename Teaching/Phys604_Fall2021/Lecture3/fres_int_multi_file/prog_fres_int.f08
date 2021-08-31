! Name: Cyrus Dreyer
! Date: 04/21/2020
!
! Purpose: Calculate the integral of the fresnel series

! * ********************************************
! * Main program to calculate fresnel integral *
! **********************************************
program fresnel_int
  implicit none

  real(8) :: x_low, x_hi ! Limits of the integral
  real(8) :: s_int_simp,s_int_trap,s_int_mid ! integral for a given nsub
  integer, parameter :: NSUBMAX=10000 ! Maximum number of subintervals
  integer :: nsub ! number of subintervals to use for the integration
  integer :: out_lun ! LUN for output file
  
  ! Prompt user for limits of the integral
  write(*,*) "Enter point to evaluate fresnel function:"
  read (*,*) x_hi
  x_low = 0 ! Lower bound for Fresnel integral is zero
  
  ! Open file for output
  open(newunit=out_lun,file='fres_int_output.dat',status='replace')
  write(out_lun,*) '# Point:', x_hi 
  write(out_lun,*) '# Number of subintervals, Midpoint, Trapezoid, Simpsons'
  
  ! Loop for the number of subintervals
  do nsub = 1, NSUBMAX

     ! Calculate the integral for nsub subintervals
     call calc_int_mid(nsub,x_low,x_hi,s_int_mid)
     call calc_int_trap(nsub,x_low,x_hi,s_int_trap)
     call calc_int_simp(nsub,x_low,x_hi,s_int_simp)

     ! Write out integral

     write(out_lun,'(i10,3e40.30)') nsub,s_int_mid,s_int_trap,s_int_simp
          
  end do ! End loop for number of subintervals

  close(out_lun)
  
  stop 0
end program fresnel_int

