# Stochastic Simulation of Malaria Epidemic

## Usage
Build the files by running ```$ make``` in the directory. Simulate with ```$ mpirun -np p ./malaria N```. 
* p is the number of processes.
* N is the number of experiments.
* ```malaria``` is the executable file
The results are printed (stdout). 

> Example: ```$ mpirun -np 4 ./malaria 4000```. This runs 4000 experiments with 4 parallel processes, 1000 experiments per process.

**Note**: OpenMPI has to be installed.


## Scalability
<p float="left">
  <img src="/images/strong.png" width=40% />
  <img src="/images/weak.png" width=40% /> 
</p>

