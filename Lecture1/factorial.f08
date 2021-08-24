program factorial

  implicit none
  integer(8) :: input ! input number to take factorial
  integer(8) :: out_fact ! output factorial result 
  integer(8) :: loop_iter ! number to iterate for loops
  
  write(*,*) "Enter an INTEGER to take the factorial of: "
  read(*,*) input

  loop_iter=0
  out_fact=1 

  do loop_iter=1,input
     out_fact=out_fact*loop_iter
  end do

  write(*,*) "The factorial of ",input,"is",out_fact

end program factorial
