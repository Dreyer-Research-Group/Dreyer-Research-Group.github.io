subroutine print_result(res)

  use compute_module, only: min_value

  implicit none

  real :: res

  write(*,*) res

end subroutine print_result
