double homogfE(double *mu, double **chimatrix, double *f){
    
    int i,j;
    double pA1_ave,pA2_ave,pA3_ave;
    double pB1_ave,pB2_ave;
    double pC_ave;
    double wA1_ave,wA2_ave,wA3_ave;
    double wB1_ave,wB2_ave;
    double wC_ave;
    double dwA1_ave,dwA2_ave,dwA3_ave;
    double dwB1_ave,dwB2_ave;
    double dwC_ave,dpp_ave;
    double eta_ave;
    double f_int, f_omeg,fE_hom;
    double *p_vect;
    double *w_vect;
    
    p_vect=create_1d_double_array(6,"p_vect");
    w_vect=create_1d_double_array(6,"w_vect");
    
    f_int=0.0;
    f_omeg=0.0;
    
    dwA1_ave=0.0;
    dwA2_ave=0.0;
    dwA2_ave=0.0;
    dwB1_ave=0.0;
    dwB2_ave=0.0;
    dwC_ave=0.0;
    
    eta_ave=0.0;
    
    pA1_ave=0.002;
    pB1_ave=pA1_ave;
    pA2_ave=0.002;
    pB2_ave=pA2_ave;
    pA3_ave=pA2_ave;
    
    
    pC_ave=1.0-(pA1_ave+pA2_ave+pA3_ave+pB1_ave+pB2_ave);
    
    wA1_ave=chimatrix[0][1]*pB1_ave+chimatrix[0][2]*pA2_ave+chimatrix[0][3]*pB2_ave+chimatrix[0][4]*pA3_ave+chimatrix[0][5]*pC_ave+eta_ave;
    
    wB1_ave=chimatrix[1][0]*pA1_ave+chimatrix[1][2]*pA2_ave+chimatrix[1][3]*pB2_ave+chimatrix[1][4]*pA3_ave+chimatrix[1][5]*pC_ave+eta_ave;
    
    wA2_ave=chimatrix[2][0]*pA1_ave+chimatrix[2][1]*pB1_ave+chimatrix[2][3]*pB2_ave+chimatrix[2][4]*pA3_ave+chimatrix[2][5]*pC_ave+eta_ave;
    
    wB2_ave=chimatrix[3][0]*pA1_ave+chimatrix[3][1]*pB1_ave+chimatrix[3][2]*pA2_ave+chimatrix[3][4]*pA3_ave+chimatrix[3][5]*pC_ave+eta_ave;
    
    wA3_ave=chimatrix[4][0]*pA1_ave+chimatrix[4][1]*pB1_ave+chimatrix[4][2]*pA2_ave+chimatrix[4][3]*pB2_ave+chimatrix[4][5]*pC_ave+eta_ave;
    
    wC_ave=chimatrix[5][0]*pA1_ave+chimatrix[5][1]*pB1_ave+chimatrix[5][2]*pA2_ave+chimatrix[5][3]*pB2_ave+chimatrix[5][4]*pA3_ave+eta_ave;
    
    
    for (i=0;i<10000000;i++){
        
        eta_ave=eta_ave-0.05*(1.0-(pA1_ave+pA2_ave+pA3_ave+pB1_ave+pB2_ave+pC_ave));
        
        pA1_ave=exp(mu[0]-wA1_ave*f[0]-wB1_ave*f[1])*f[0];
        pB1_ave=exp(mu[0]-wA1_ave*f[0]-wB1_ave*f[1])*f[1];
        
        pA2_ave=exp(2.0*mu[1]-wA2_ave*f[0]/2.0-wB2_ave*f[1]-wA3_ave*f[0]/2.0)*f[0]/2.0;
        pB2_ave=exp(2.0*mu[1]-wA2_ave*f[0]/2.0-wB2_ave*f[1]-wA3_ave*f[0]/2.0)*(2.0*f[1])/2.0;
        pA3_ave=exp(2.0*mu[1]-wA2_ave*f[0]/2.0-wB2_ave*f[1]-wA3_ave*f[0]/2.0)*f[0]/2.0;
        
        pC_ave=exp(kappa*(mu[2]-wC_ave));
        
        dwA1_ave=(chimatrix[0][1]*pB1_ave+chimatrix[0][2]*pA2_ave+chimatrix[0][3]*pB2_ave+chimatrix[0][4]*pA3_ave+chimatrix[0][5]*pC_ave+eta_ave)-wA1_ave;
        dwB1_ave=(chimatrix[1][0]*pA1_ave+chimatrix[1][2]*pA2_ave+chimatrix[1][3]*pB2_ave+chimatrix[1][4]*pA3_ave+chimatrix[1][5]*pC_ave+eta_ave)-wB1_ave;
        
        dwA2_ave=(chimatrix[2][0]*pA1_ave+chimatrix[2][1]*pB1_ave+chimatrix[2][3]*pB2_ave+chimatrix[2][4]*pA3_ave+chimatrix[2][5]*pC_ave+eta_ave)-wA2_ave;
        dwB2_ave=(chimatrix[3][0]*pA1_ave+chimatrix[3][1]*pB1_ave+chimatrix[3][2]*pA2_ave+chimatrix[3][4]*pA3_ave+chimatrix[3][5]*pC_ave+eta_ave)-wB2_ave;
        dwA3_ave=(chimatrix[4][0]*pA1_ave+chimatrix[4][1]*pB1_ave+chimatrix[4][2]*pA2_ave+chimatrix[4][3]*pB2_ave+chimatrix[4][5]*pC_ave+eta_ave)-wA3_ave;
        
        dwC_ave=(chimatrix[5][0]*pA1_ave+chimatrix[5][1]*pB1_ave+chimatrix[5][2]*pA2_ave+chimatrix[5][3]*pB2_ave+chimatrix[5][4]*pA3_ave+eta_ave)-wC_ave;
        
        dpp_ave=1.0-(pA1_ave+pA2_ave+pA3_ave+pB1_ave+pB2_ave+pC_ave);
        
        wA1_ave=wA1_ave+0.005*dwA1_ave;
        wB1_ave=wB1_ave+0.005*dwB1_ave;
        wA2_ave=wA2_ave+0.005*dwA2_ave;
        wB2_ave=wB2_ave+0.005*dwB2_ave;
        wA3_ave=wA3_ave+0.005*dwA3_ave;
        wC_ave=wC_ave+0.005*dwC_ave;
        
    }
    //not sure if I will need these
    //phiAB_hom=pA_ave+pB_ave;
    //phiC_hom=pC_ave;
    
    
    p_vect[0]=pA1_ave;
    p_vect[1]=pB1_ave;
    p_vect[2]=pA2_ave;
    p_vect[3]=pB2_ave;
    p_vect[4]=pA3_ave;
    p_vect[5]=pC_ave;

    
    w_vect[0]=wA1_ave;
    w_vect[1]=wB1_ave;
    w_vect[2]=wA2_ave;
    w_vect[3]=wB2_ave;
    w_vect[4]=wA3_ave;
    w_vect[5]=wC_ave;
    
    
    
    for (i=0;i<6;i++){
        for (j=0;j<6;j++){
            f_int+=p_vect[i]*p_vect[j]*chimatrix[i][j];
        }
        f_omeg+=p_vect[i]*w_vect[i];
    }
    
    
    fE_hom=f_int/2.0-f_omeg-(exp(mu[0]-wA1_ave*f[0]-wB1_ave*f[1]));
    fE_hom-=(exp(2.0*mu[1]-wA2_ave*f[0]-wA3_ave*f[0]-2.0*wB2_ave*f[1])/2.0);
    fE_hom-=(exp(kappa*(mu[2]-wC_ave))/kappa);
    
    return fE_hom;
}

//This one is wrong I think
double homofE(double **chiMatrix){
    
    int i,j;
    int NABA_triblock,NAB_diblock,NC_homopolymer;
    int *Ls;
    double fE_homo;
    double *p_vect;
    double *phiAve;
    
    p_vect=create_1d_double_array(6, "p_vect");
    phiAve=create_1d_double_array(3, "phiAve");
    Ls=create_1d_integer_array(6, "Ls");
    
    
    phiAve[0]=0.0001; // ABC Triblock
    phiAve[1]=0.9998; // DE Diblock
    phiAve[2]=1.0-(phiAve[0]+phiAve[1]); // F Homopolymer
    
    // Setting the chain lengths________________
    Ls[0]=50;//NA1
    Ls[1]=50;//NB1
    Ls[2]=50;//NA2
    Ls[3]=100;//NB2
    Ls[4]=50;//NA3
    Ls[5]=100;//NC
    NABA_triblock=Ls[2]+Ls[3]+Ls[4];
    NAB_diblock=Ls[0]+Ls[1];
    NC_homopolymer=Ls[5];
    
    // Chains measured with respect to diblock

    
    p_vect[0]=phiAve[1]*((double)Ls[0]/(double)NAB_diblock); // A1
    p_vect[1]=phiAve[1]*((double)Ls[1]/(double)NAB_diblock); // B1
    p_vect[2]=phiAve[0]*((double)Ls[2]/(double)NABA_triblock); // A2
    p_vect[3]=phiAve[0]*((double)Ls[3]/(double)NABA_triblock); // B2
    p_vect[4]=phiAve[0]*((double)Ls[4]/(double)NABA_triblock); // A3
    p_vect[5]=phiAve[2]; // C
    
    fE_homo=0.0;
    for(i=0;i<ChainType;i++){
        for(j=0;j<ChainType;j++){
            fE_homo+=p_vect[i]*p_vect[j]*chiMatrix[i][j];
        }
    }
    fE_homo/=2.0;
    
    destroy_1d_double_array(phiAve);
    destroy_1d_double_array(p_vect);
    destroy_1d_integer_array(Ls);
    
    return fE_homo;
    
    
}
