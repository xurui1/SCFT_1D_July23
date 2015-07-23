double fE(double **newW, double **phi, double **chiMatrix, double dr, double volume){
    
    int i,ii,jj;
    double fEW,fEchi;
    double fE_int;
    fEW=0.0;
    fEchi=0.0;
    fE_int=0.0;
    //ends
    i=0;
    for(ii=0;ii<ChainType;ii++){
        for(jj=0;jj<ChainType;jj++){
            fEchi+=0.5*phi[ii][i]*chiMatrix[ii][jj]*phi[jj][i]*dV(i,dr);
        }
        fEW+=0.5*(newW[ii][i]*phi[ii][i]*dV(i,dr));
    }
    
    
    i=(int)Nr-1;
    for(ii=0;ii<ChainType;ii++){
        for(jj=0;jj<ChainType;jj++){
            fEchi+=0.5*phi[ii][i]*chiMatrix[ii][jj]*phi[jj][i]*dV(i,dr);
        }
        fEW+=0.5*(newW[ii][i]*phi[ii][i]*dV(i,dr));
    }

    
    //middle
    for (i=1;i<(int)Nr-1;i++){
        for(ii=0;ii<ChainType;ii++){
            for(jj=0;jj<ChainType;jj++){
                fEchi+=phi[ii][i]*chiMatrix[ii][jj]*phi[jj][i]*dV(i,dr);
            }
            fEW+=(newW[ii][i]*phi[ii][i]*dV(i,dr));
        }
    }

    //normalize by box size
    fEchi/=2.0*volume;
    fEW/=volume;
    
    fE_int=fEchi-fEW;

    return fE_int;
    
}