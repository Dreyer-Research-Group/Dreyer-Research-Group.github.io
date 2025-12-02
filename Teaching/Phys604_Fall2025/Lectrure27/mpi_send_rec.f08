program send_rec

  use mpi

  implicit none

  integer :: ierr, rank, nprocs, count, status(MPI_STATUS_SIZE)
  integer :: my_array(5)

  call MPI_Init(ierr)

  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, nprocs, ierr)

  my_array=0

  if (rank == 0) then

     my_array=[1,2,3,4,5]
 
     call MPI_Send(my_array,5, MPI_INT,3, 1, MPI_COMM_WORLD, ierr)
     
  else if (rank == 3) then
     
     call MPI_Recv(my_array,5, MPI_INT,0, 1, MPI_COMM_WORLD, status,ierr)

     write(*,*) "Source rank = ", status(MPI_SOURCE)
     write(*,*) "Tag         = ", status(MPI_TAG)

     call MPI_Get_count(status, MPI_INTEGER, count, ierr)
     write(*,*) "Number of elements received = ", count

     
  end if
    
  write(*,*) "Process:",rank, "my_array:",my_array

  call MPI_Finalize(ierr)
  
end program send_rec
