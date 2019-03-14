int row = opt->process;
int col = mod->K;

double resContainer[row][col];  //need make this global, so in the main can sort and get the best result

//initial
for(int i = 0; i < row; i++){
    for(int j = 0 ; i < col; j++){
        resContainer[i][j] = 0;
    }
}

int processCount[row];
for(int i = 0 ; i < row ;i ++){
    processCount[i] = 0;
}

//rank=0, receive and save 
if(my_rank == 0){
    for(source = 1; source < p; source++){
        cout<<"process 0 recieving " //(change to c print)
        currRound = processCount[source];//get the index;
        //need make sure what's the out put of the res, a number or an array
        MPI_Recv(&ProcessNewRes,1,MPI_FLOAT,source,tag,MPI_COMM_WORLD,&status);
        resContainer[source][currRound] = ProcessNewRes; //save the res
        processCount[source]++; //increase the round_count of the process;
    }
}
//rank != 0 send the res;
else{
    count<<"Process "<<my_rank<<" sending"; // chang to c print
    MPI_Send(&ProcessNewRes,1,MPI_FLOAT,des,tag,MPI_COMM_WORLD);
}

//bubbort sort
if(my_rank == 0){
    for(int i = 0; i< col; i++){ //for each col
        for(int j = 0 ; j<row; j++){ //for each column
            if()
        }
    }
}




