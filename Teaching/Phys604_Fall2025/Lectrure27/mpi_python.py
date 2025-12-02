from mpi4py import MPI
import numpy as np


comm = MPI.COMM_WORLD       # Communicator
rank = comm.Get_rank()      # Process rank
size = comm.Get_size()      # Number of processes

print(f"Hello from process {rank} of {size}")

#if rank==0:
#    data=np.array([1,2,3,4],dtype='int64')
#    comm.Send(data, dest=2)
#elif rank==2:
#    data=np.zeros(4,dtype='int64')
#    comm.Recv(data, source=0)

#    print(f"data: {data}")

#if rank==0:
#    data=np.array([1,2,3,4],dtype='int64')
#    comm.send(data, dest=2)
#elif rank==2:
#    data=comm.recv(source=0)
#    print(f"data: {data}")

#data = comm.bcast({'x':2,'y':3},root=0)

#print(f"Process {rank}, dict {data}")

comm.Barrier()

if rank==0:
    data={'x':2,'y':3}
    transfer=np.array([data['x'],data['y']],dtype='int64')
else:
    data={'x':5,'y':6}
    transfer=np.zeros(2,dtype='int64')
    
comm.Bcast(transfer,root=0)

if rank != 0:
    data['x']=transfer[0]
    data['y']=transfer[1]

comm.Barrier()

print(f"Process {rank}, dict {data}")
