Parallelize n-particle simulation which is expected to run in time T = O(n) on a single processor and in time T/p when using p processors. Write parallel codes that approach these expectations

Parallel programming for performance with shared memory and message passing

C with Pthreads, OpenMP and MPI

Given four particle simulation programs of O(n2) time complexity, you are to improve performance of the programs, i.e. develop and to evaluate the following four particle simulation programs

    A sequential program that runs in time T = O(n), where n is the number of particles.
    A parallel program using Pthreads that runs in time close to T/p when using p processors.
    A parallel program using OpenMP (with/without tasks) that runs in time close to T/p when using p processors.
    A parallel program using MPI that runs in time close to T/p when using p processors.
