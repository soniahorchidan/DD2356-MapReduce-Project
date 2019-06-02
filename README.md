# DD2356-MapReduce-Project
The final project of the DD2356 Methods in High Performance Computing at KTH Royal Institute of Technology. 

The goal of the project was to implement a <b>MapReduce</b> framework in a high performance context and apply it on a specific use-case: <b>WordCount</b>. 

The design includes <b>MPI</b> non-blocking  point-to-point  communication  and  collective  MPI I/O.  We  also  explored  methods  to  parallelize  the  tasks  even further, at each process’ level, by integrating <b>OpenMP</b>.

To run our code on Beskow:
- Login with your credentials
- ```module swap PrgEnv-cray/5.2.82 PrgEnv-gnu```
- ```cd to/any/folder```
- ```git clone https://github.com/SoniaHorchidan/DD2356-MapReduce-Project.git```
- ```cd DD2356-MapReduce-Project```
- ```make```
- ```aprun -n 32 -N 32 ./bin/project.out -r 1 /cfs/klemming/scratch/s/sergiorg/DD2356/input/wikipedia_20GB.txt```
