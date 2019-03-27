#include <iostream>
#include <stdio.h>
#include <string>
#include <stdio.h>
#include <bits/stdc++.h>
using namespace std;
#include "mpi.h"

main(int argc, char** argv){
    
    int my_rank;    //the index of the process
    int p;          //opt->process
    
    //initialize the MPI environment and configuration
    MPI_Init(&argc, &argv);     
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // Get the rank of the process

    int k = 10;     //mod->K
    float dest = 0;
    int tag = 1;
    int source;     //Process sending the possibility 
    MPI_Status status;

    int max_processor = 0;
    int max_logL = INT_MIN;
    int amimax = 0; //int is 1 if processor is max, it's 0 otherwise

    //process do something and get a result for sending;
    srand(time(0)*my_rank);
    float result = rand()/k;
   
    if(my_rank == 0){
       // for(int countK = 0; countK < k; countK++ ){           
            for(int proc = 1; proc < p; proc++){
                cout<<"Process 0 is recieving "<<endl;
                MPI_Recv(&result, 1, MPI_FLOAT, proc, tag*proc, MPI_COMM_WORLD, &status);
                if(max_logL < result){
                    max_logL = result;
                    max_processor = proc;
                }
            }
            cout<<"the biggest one is"<<max_logL<<endl;
            cout<<"the processor with biggest logL is"<<max_processor<<endl;
            //Head node here knows which processor has the maxL
            //So it has to send that info to every processor, so they can proceed
            for (int proc = 1; proc < p; proc++){
                if (proc == max_processor)
                    amimax = 1;
                MPI_Send(&amimax, 1, MPI_INT, proc, tag*proc*7, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);  
       // }
    }
    else
    {
        cout<<"Process "<<my_rank<<" sending "<<endl;
        MPI_Send(&result, 1, MPI_FLOAT, dest, tag*my_rank, MPI_COMM_WORLD);
        cout<<"my result = " << result;

        MPI_Barrier(MPI_COMM_WORLD);
        
        MPI_Recv(&amimax, 1, MPI_INT, 0, tag*my_rank*7, MPI_COMM_WORLD, &status);
        
        if (amimax == 1)
            cout<< "Hey I am processor "<<my_rank<<" and I have the max logl"<<endl;
        else if (amimax == 0)
            cout<< "I am processor "<<my_rank<<" and I don't have the max logl :("<<endl;
        //MPI_Barrier(MPI_COMM_WORLD);
        /*
        if(my_rank = max_k){
            //do printout job;
            cout << "!!! I'm the biggest" <<endl; 
        }
        */
        
    }

    MPI_Finalize();
    return 0;
}
