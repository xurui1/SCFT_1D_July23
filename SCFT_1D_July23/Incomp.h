void Incomp(double *eta, double **phi, double *delphi){
    
    int     i;
    int     chain;
    double  ptot;
    
    ptot=0.0;
    
    for(i=0;i<Nr;i++){
        
        ptot=0.0;
        delphi[i]=0.0;
                
        for(chain=0;chain<ChainType;chain++){
            ptot+=phi[chain][i];
        }
                            
        delphi[i]=1.0-ptot;
        eta[i]-=delphi[i];
        
        if (fabs(delphi[i])>1e2){
            cout<<i<<" incomp: "<<delphi[i]<<endl;
            exit(EXIT_FAILURE);
        }
    }
}


void Pin(double *sigma,double **phi){
    
    int Ntip;
    int i;
    
    Ntip=2*Nr/5;
    
    for (i=0;i<Nr;i++){
        if (i==Ntip){
            sigma[i]=sigma[i]-10.0*(phi[0][i]+phi[2][i]+phi[4][i]-phi[1][i]-phi[3][i]);
        }
        else{
            sigma[i]=0.0;
        }
    }
}