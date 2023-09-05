module compute_module

  implicit none

  real :: min_value

contains
  real function compute(xx, yy)

    real :: xx,yy
    
    min_value = min(xx,yy)
    compute = xx + yy

    return
  end function compute

end module compute_module
