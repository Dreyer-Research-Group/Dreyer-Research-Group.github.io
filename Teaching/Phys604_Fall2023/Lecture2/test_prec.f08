! Purpose: Test the precision of reals
! Author: Cyrus Dreyer
! Date: 2/4/2019 
program test_prec_reals
 implicit none ! Turn off implicit typing
 ! Variable dictionary
 real :: factor1 ! Variable for factor 1 
 real :: factor2 ! Variable for factor 2  
 real :: prec_test_lhs ! Variable for result 
 real :: prec_test_rhs ! Variable for result
 factor1 = 1.0 ! Assign a value to factor1 
 factor2 = 1.0d-9 ! Assign a value to factor2 

 prec_test_lhs = (factor1-factor1) + factor2 
 prec_test_rhs = factor1 + (-factor1 + factor2) 
 ! Output
 write(*,'(a20,e20.12e2,a20,e20.12e2)') "Prec_test_lhs:",prec_test_lhs, &
      "Prec_test_rhs:", prec_test_rhs

end program test_prec_reals
