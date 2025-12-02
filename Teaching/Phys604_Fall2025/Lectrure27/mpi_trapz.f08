  ! Purpose:  Integrate a function with the trapezoid rule
  ! Author:    Cyrus Dreyer
  ! Date:        12/01/2025
program trapezoid

  use mpi

  implicit none                

  integer :: ierr, nprocs, rank ! For MPI
  integer :: nsub ! number of subintervals
  integer :: loop, nsubs ! For the loop over subintervals
  integer, allocatable :: nsub_per_proc(:) ! For the loop over subintervals 
  real(8) :: xhi,xlow,dx,xi,xlow_count ! Limits of integration and other x points
  real(8) :: sum_part,sum_final ! partial and final sum
  real(8) :: f ! function 
  real(8) :: t_start,t_finish ! for timing
  real(8),allocatable :: xlows(:)

  ! Set up the MPI calculation
  call MPI_Init(ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr)

  ! Time it all 
  !t_start=MPI_WTIME()

  ! Get some parameters from the user
  if (rank == 0) then
     write(*,*) "Enter the lower and upper bounds:"
     read(*,*) xlow,xhi
     write(*,*) "Enter the number of subintervals:"
     read(*,*) nsub

     dx = (xhi-xlow)/(1.0*real(nsub))
     nsubs=floor(real(nsub)/real(nprocs))

     ! Determine how many subintervals to do on each processor
     allocate(nsub_per_proc(0:nprocs-1))
     nsub_per_proc=0
     do loop=0,nsub-1
        nsub_per_proc(mod(loop,nprocs))=nsub_per_proc(mod(loop,nprocs))+1
     end do

     ! Send the processes their sections of the range
     xlow_count=xlow
     allocate(xlows(0:nprocs-1))
     do loop = 0,nprocs-1
        xlows(loop)=xlow_count
        xlow_count=xlow_count+dx*nsub_per_proc(loop)
     end do

  end if

  ! Send out xhi, dx and nsub to all processors
  call MPI_Bcast(xhi,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(dx,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  ! Scatter the nsub_per_proc and xlows
  call MPI_Scatter(nsub_per_proc,1,MPI_INT,nsub,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(xlows,1,MPI_REAL8,xlow,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  !TEST: Write out info on each rank
  !write(*,*) 'rank: ',rank,'nsub: ',nsub,'xlow: ', xlow

  ! Time just the loop 
  t_start=MPI_WTIME()

  
  sum_part = 0.0                     
  xi = xlow                     
  do loop = 1,nsub
     sum_part = sum_part+0.5*(f(xi)+f(xi+dx))*dx 
     xi = xi+dx                 
     
     ! In case we are at the end of the range
     if (xi+dx > xhi) then
        exit
     end if

  enddo                         
 
  ! Gather the results:
  call MPI_Reduce(sum_part,sum_final,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)   

  t_finish=MPI_WTIME()
  
  if (rank == 0) then
    write(*,*) ' sum = ',sum_final
    write(*,*) ' time= ',t_finish-t_start

 end if

 call MPI_Finalize(ierr)

end program trapezoid

real(8) function f(x)

  implicit none

  real(8),intent(in) :: x

  f=x*sin(x)

end function f
