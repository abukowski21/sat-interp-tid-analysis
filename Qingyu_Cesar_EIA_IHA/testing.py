# mpi_hello.py
from mpi4py import MPI

# Initialize the MPI communicator
comm = MPI.COMM_WORLD
# Get the total number of processes
size = comm.Get_size()
# Get the rank (unique ID) of this process
rank = comm.Get_rank()
# Get the name of the processor
name = MPI.Get_processor_name()

print(f"Hello from rank {rank} out of {size} processors on {name}")

# Ensure all processes finish before exiting
comm.Barrier()
