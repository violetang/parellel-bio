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
    int tag = 0;
    int source;     //Process sending the possibility 
    MPI_Status status;;

    int max_k = 0;
    int max_logL = INT_MIN;

    //process do something and get a result for sending;
    srand(time(0)*my_rank);
    float result = rand()/2;
   
    if(my_rank == 0){
       // for(int countK = 0; countK < k; countK++ ){
            
            for(int proc = 1; proc < p; proc++){
                cout<<"Process 0 is recieving "<<endl;
                MPI_Recv(&result, 1, MPI_FLOAT, proc, tag, MPI_COMM_WORLD, &status);
                if(max_logL < result){
                    max_logL = result;
                    max_k = proc;
                }
            }
            cout<<"the biggest one is"<<max_logL<<endl;
            cout<<"the biggest p is"<<max_k<<endl;
           // MPI_Send(&max_k,1,MPI_INT,max_k,tag,MPI_COMM_WORLD); //need change to send all
            MPI_Barrier(MPI_COMM_WORLD);  
       // }
    }
    else
    {
        cout<<"Process "<<my_rank<<" sending "<<endl;
        MPI_Send(&result, 1, MPI_FLOAT, dest,tag, MPI_COMM_WORLD);
        cout<<"my result = " << result;
        MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Recv(&max_k,1,MPI_INT,dest,tag, MPI_COMM_WORLD);
    
        if(my_rank = max_k){
            //do printout job;
            cout << "!!! I'm the biggest" <<endl; 
        }
        
    }

    MPI_Finalize();
    return 0;
}
