# DD2356-MapReduce-Project
The final project of the DD2356 Methods in High Performance Computing at KTH Royal Institute of Technology. 

The goal of the project was to implement a <b>MapReduce</b> framework in a high performance context and apply it on a specific use-case: <b>WordCount</b>. 

The design includes <b>MPI</b> non-blocking  point-to-point  communication  and  collective  MPI I/O.  We  also  explored  methods  to  parallelize  the  tasks  even further, at each process’ level, by integrating <b>OpenMP</b>.

To run our code on Beskow:
1 Login with your credentials
2 ```module swap PrgEnv-cray/5.2.82 PrgEnv-gnu```
3 ```cd to_any_folder```
4 ```git clone https://github.com/SoniaHorchidan/DD2356-MapReduce-Project.git```
5 ```cd DD2356-MapReduce-Project```
6 ```make```
7 ```aprun ```
