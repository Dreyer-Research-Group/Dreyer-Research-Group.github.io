program exp_series_neg_24

  implicit none
  real(8) :: input ! input number to take factorial
  real(8) :: final_exp
  real(8) :: exp_sum ! S(x)
  real(8) :: term_exp,trunc_err ! Term in the expansion
  real(8) :: fact_exp ! output factorial result 
  integer :: loop1 ! number to iterate for loops
  
  ! Prompt the user for number to take the factorial
 ! write(*,*) "Enter a real number to take the exponential of "
  !read(*,*) input  

  input=-1.0
  
  ! Initialize the variables  
  loop1=0
  exp_sum=0.0

  
  do
     
     ! Iterate factorial
     if (loop1==0) then
        fact_exp=1
     else
        fact_exp=fact_exp*loop1
     end if
     
     ! Iterate exp sum
     term_exp=((input)**loop1)/fact_exp

     ! Make sure we have not overflowed (or NaN)
     if (abs(term_exp)>=huge(term_exp).or. &
          & abs(term_exp)<=tiny(term_exp).or.term_exp /= term_exp) then
        write(*,*) "Overflow!"
        loop1=loop1-1
        exit
     else
        trunc_err=term_exp
     end if
     
     ! Test for convergence to machine precision
     if (exp_sum+term_exp /= exp_sum) then
        exp_sum=exp_sum+term_exp
     else
        exit
     end if

     ! Write out progress to the terminal
     write(*,*) "Term in sum:",term_exp, "Sum^24:",exp_sum**24
     
     ! iterate loop variable
     loop1=loop1+1     
     
  end do

  write(*,*)
  write(*,*) "After ",loop1, "steps:"
  write(*,*) "Truncation error:", trunc_err
  write(*,*) "Result from series:",exp_sum**24
  write(*,*) "Result from built-in function:", exp(-24.0)
  
  
  stop 0
  
end program exp_series_neg_24
