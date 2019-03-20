int count_Pro = opt->process;  //row == my_rank;
int count_k = mod->K;   //col = K's

int max_k;
int max_logL;

//rank=0, receive and save 
if(my_rank == 0){
    for(int countK = 0; countK < col; countK++){ //for each K
        for(int pro = 1; pro < p; pro++){ //for each process result
            cout<<"process 0 recieving " //(change to c print)
            MPI_Recv(&ProcessNewRes,1,MPI_FLOAT,source,tag,MPI_COMM_WORLD,&status);
            if(max_logL < ProcessNewRes){
                max_logL = ProcessNewRes;
                max_k = source;
            }
       }
       //sort get the hihgest one, tell the process to print out
       MPI_Barrier(MPI_COMM_WORLD);
    } 
}
//rank != 0 send the res;
else{
    count<<"Process "<<my_rank<<" sending"; // chang to c print
    MPI_Send(&ProcessNewRes,1,MPI_FLOAT,des,tag,MPI_COMM_WORLD);
    //if(my_rank == highest, print result)
    MPI_Barrier(MPI_COMM_WORLD);
}

//bubbort sort
if(my_rank == 0){
    for(int i = 0; i< col; i++){ //for each col
        for(int j = 0 ; j<row; j++){ //for each column
            if()
        }
    }
}




