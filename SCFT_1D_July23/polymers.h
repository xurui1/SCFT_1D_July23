void diblock(double **qA1,double **qdagA1,double **qB1,double **qdagB1,double **w,double ds,int *Ns,double dr){
    
    int i;
    // Here is the for loop for setting the propagator initial conditions to 1.0
    for(i=0;i<Nr;i++){
        qA1[i][0]=1.0;
    }
    
    // Here we solve the diffusion equation for the forwards propagators
    solvediffyQ(qA1,w[0],ds,Ns[0],dr);
    
    for(i=0;i<Nr;i++){
        qB1[i][0]=qA1[i][Ns[0]];
    }
    solvediffyQ(qB1,w[1],ds,Ns[1],dr);
    
    //Complementary propagator
    for(i=0;i<Nr;i++){
        qdagB1[i][0]=1.0;
    }
    solvediffyQ(qdagB1,w[1],ds,Ns[1],dr);
    
    for(i=0;i<Nr;i++){
        qdagA1[i][0]=qdagB1[i][Ns[1]];
    }
    solvediffyQ(qdagA1,w[0],ds,Ns[0],dr);
    
}

void triblock(double **qA2,double **qB2,double **qA3,double **w,double ds,int *Ns,double dr){
    //Triblock is symmetric so we don't need complementary propagator
    int i;
    
    for (i=0;i<Nr;i++){
        qA2[i][0]=1.0;
    }
    
    solvediffyQ(qA2,w[2],ds,Ns[0],dr);
    
    for (i=0;i<Nr;i++){
        qB2[i][0]=qA2[i][Ns[0]];
    }
    solvediffyQ(qB2,w[3],ds,2*Ns[1],dr);
    
    for (i=0;i<Nr;i++){
        qA3[i][0]=qB2[i][2*Ns[1]];
    }
    solvediffyQ(qA3,w[4],ds,Ns[0],dr);
        
}

void homopolymer(double **qC,double **w,double ds,int *Ns,double dr){
    int i;
    for (i=0;i<Nr;i++){
        qC[i][0]=1.0;
    }
    
    solvediffyQ(qC,w[5],ds,Ns[2],dr);
    
}

