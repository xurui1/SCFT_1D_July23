double Conc(double **phi,double **w,int *Ns,double ds,double dr, double *mu,double volume){
    
    double      Q;
    double      **qA1,**qA2,**qA3,**qdagA1;
    double      **qB1,**qB2,**qdagB1;
    double      **qC;
    
    
    //Forwards propagators
    qA1=create_2d_double_array(Nr,Ns[0]+1,"qA1");
    qB1=create_2d_double_array(Nr,Ns[1]+1,"qB1");
    qA2=create_2d_double_array(Nr,Ns[0]+1,"qA2");
    qB2=create_2d_double_array(Nr,2*Ns[1]+1,"qB2");
    qA3=create_2d_double_array(Nr,Ns[0]+1,"qA3");
    qC=create_2d_double_array(Nr,Ns[2]+1,"qC");
    
    //Complementary propagators
    qdagA1=create_2d_double_array(Nr,Ns[0]+1,"qdagA1");
    qdagB1=create_2d_double_array(Nr,Ns[1]+1,"qdagB1");
    
    diblock(qA1,qdagA1,qB1,qdagB1,w,ds,Ns,dr);
    triblock(qA2,qB2,qA3,w,ds,Ns,dr);
    homopolymer(qC,w,ds,Ns,dr);
    
    // Here we get the single chain partition functions Q_AB+Q_C
    Q=q_partition(qB1,qA3,qC,dr,Ns,mu,volume);
        
    //cout<<"Q: "<< Q<<endl;
    
    // Here we do the concentration calculation by integration over chain
    phi_calc(phi,qA1,qdagA1,qB1,qdagB1,qA2,qB2,qA3,qC,Ns,mu,ds);

    //calculation of average concentrations over entire computation box
    phi_total(phi,dr,volume);
    
    loop(qA2,qA3,Ns,phi,volume,dr);
    
    
    //clearing the memory
    destroy_2d_double_array(qA1);
    destroy_2d_double_array(qB1);
    destroy_2d_double_array(qA2);
    destroy_2d_double_array(qB2);
    destroy_2d_double_array(qA3);
    destroy_2d_double_array(qC);
    destroy_2d_double_array(qdagA1);
    destroy_2d_double_array(qdagB1);
    
    return Q;
    
}