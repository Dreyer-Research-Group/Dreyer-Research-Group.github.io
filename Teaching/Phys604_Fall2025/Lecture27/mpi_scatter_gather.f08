program scatter_gather

  use mpi

  implicit none

  integer :: ierr, rank, nprocs,status(MPI_STATUS_SIZE)
  integer :: my_array(10),data,sum_data

  call MPI_Init(ierr)

  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, nprocs, ierr)

  my_array=0

  if (rank == 0) then

     my_array=[0,1,2,3,4,5,6,7,8,9]
     write(*,*) "my_array:",my_array

  end if
 
  call MPI_Scatter(my_array,1,MPI_INT,data,1,MPI_INT,0,MPI_COMM_WORLD,ierr)   
    
  write(*,*) "Process:",rank, "data:",data

  data=data*rank**2

  call MPI_Gather(data,1,MPI_INT,my_array,1,MPI_INT,0,MPI_COMM_WORLD,ierr)  
  if (rank == 0) then
     write(*,*) "my_array:",my_array
  end if

  call MPI_Reduce(data,sum_data,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (rank == 0) then
     write(*,*) "sum_data:",sum_data
  end if


  call MPI_Finalize(ierr)
  
end program scatter_gather
