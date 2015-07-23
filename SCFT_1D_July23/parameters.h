void parameters(double *chi,double *f,double *ds,int *Ns,double *dr,double *mu){
    
    int Ds=100;
    r_0=1.0;
    double delr;
    
    initial=2;
    Coord=3; //if 1->Cartesian, if 2->Cylindrical, if 3->Spherical coordinate system
    
    //Length ratio of c homopolymer to diblock copolymer
    kappa=1.0;
        
    //Interaction parameters
    chi[0]=25.0;        //Chi_AB
    chi[1]=25.0;        //Chi_BC
    chi[2]=0.0;         //Chi_AC
    
    //Chemical potential array
    mu[0]=0.0;      //AB
    mu[1]=-20.0;      //ABA
    mu[2]=-4.2612;    //C
    
    ifstream inputfA;
    inputfA.open("fA.dat");
    inputfA >> f[0];
    inputfA.close();
    
    //Chain fraction array
    f[1]=1.0-f[0];  //B
    f[2]=kappa*1.0; //C
    
    //Chain length array
    Ns[0]=(Ds*f[0]);            //A blocks
    Ns[1]=(Ds*f[1]);            //B blocks
    Ns[2]=(Ds*f[2]);            //C blocks
    
    //cout<<Ns[0]<<" "<<Ns[1]<<" "<<Ns[2]<<endl;
    
    //Step size in r,z direction
    *dr=0.12;
    delr=*dr;
    
    
    //Step length along polymer
    *ds=1.0/Ds;
    
}

void Xmatrix(double **chiMatrix, double *chi){
    //Interaction Matrix
    chiMatrix[0][0]=0.0;    //ChiA1,A1
    chiMatrix[0][1]=chi[0]; //ChiA1,B1
    chiMatrix[0][2]=0.0;    //ChiA1,A2
    chiMatrix[0][3]=chi[0]; //ChiA1,B2
    chiMatrix[0][4]=0.0;    //ChiA1,A3
    chiMatrix[0][5]=chi[2]; //ChiA1,C
    
    chiMatrix[1][0]=chi[0]; //ChiB1,A1
    chiMatrix[1][1]=0.0;    //ChiB1,B1
    chiMatrix[1][2]=chi[0]; //ChiB1,A2
    chiMatrix[1][3]=0.0;    //ChiB1,B2
    chiMatrix[1][4]=chi[0]; //ChiB1,A3
    chiMatrix[1][5]=chi[1]; //ChiB1,C
    
    chiMatrix[2][0]=0.0;    //ChiA2,A1
    chiMatrix[2][1]=chi[0]; //ChiA2,B1
    chiMatrix[2][2]=0.0;    //ChiA2,A2
    chiMatrix[2][3]=chi[0]; //ChiA2,B2
    chiMatrix[2][4]=0.0;    //ChiA2,A3
    chiMatrix[2][5]=chi[2]; //ChiA2,C
    
    chiMatrix[3][0]=chi[0]; //ChiB2,A1
    chiMatrix[3][1]=0.0;    //ChiB2,B1
    chiMatrix[3][2]=chi[0]; //ChiB2,A2
    chiMatrix[3][3]=0.0;    //ChiB2,B2
    chiMatrix[3][4]=chi[0]; //ChiB2,A3
    chiMatrix[3][5]=chi[1]; //ChiB2,C
    
    chiMatrix[4][0]=0.0;    //ChiA3,A1
    chiMatrix[4][1]=chi[0]; //ChiA3,B1
    chiMatrix[4][2]=0.0;    //ChiA3,A2
    chiMatrix[4][3]=chi[0]; //ChiA3,B2
    chiMatrix[4][4]=0.0;    //ChiA3,A3
    chiMatrix[4][5]=chi[2]; //ChiA3,C
    
    
    chiMatrix[5][0]=chi[2]; //ChiC,A1
    chiMatrix[5][1]=chi[1]; //ChiC,B1
    chiMatrix[5][2]=chi[2]; //ChiC,A2
    chiMatrix[5][3]=chi[1]; //ChiC,B2
    chiMatrix[5][4]=chi[2]; //ChiC,A3
    chiMatrix[5][2]=0.0;    //ChiC,C
    
}
