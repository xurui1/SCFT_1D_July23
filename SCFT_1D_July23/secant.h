
void secant(double **w, double **phi, double *eta, int *Ns, double ds, double *chi, double dr, double **chiMatrix, double *mu, double *f){
    
    
    double  currentfE, oldfE, deltafE;
    int     maxIter=100000;
    double precision=1.0e-5;          //convergence condition
    double mu1,mu2,mu3;
    double fE1,fE2,fE3;
    int     i,xx,iter,ii,jj;
    int     mmb;
    int coordinate;
    double  Q;
    double  fE_int, fES;            //interaction free energy and chain partition function fE
    double  epsilon, gamma;
    double  *delphi;
    double *sigma;
    double  **delW;
    double  **newW;
    double  deltaW;
    double fE_hom=0.0;
    double volume;
    
    //Arrays for updating the omega fields
    delW=create_2d_double_array(ChainType,Nr,"delW");
    delphi=create_1d_double_array(Nr,"delphi");
    sigma=create_1d_double_array(Nr,"sigma");
    newW=create_2d_double_array(ChainType,Nr,"newW");
    
    currentfE=0.0;
    deltafE=0.0;
    
    mu1=mu[2];
    mu2=mu[2]+0.01;
    
    epsilon=0.05;
    gamma=0.05;
    
    iter=0;
    mmb=1;
    coordinate=Coord;
    Coord=1;
    
    for (xx=0;xx<3;xx++){
        fE_hom=homogfE(mu,chiMatrix,f);
        volume=vol(dr);

    for (iter=0;iter<maxIter;iter++){
        
        fE_int=0.0;
        fES=0.0;
        deltaW=0.0;
        
        
        Q=Conc(phi,w,Ns,ds,dr,mu,volume);          //Calculate Chain partition function for both AB and C
        
        
        Incomp(eta,phi,delphi);              //Enforce incompressibility condition
        output(dr,phi,w);                   //Output some data to file
        
        if (mmb==1){Pin(sigma, phi);}
        
        
        
        //Calculate components for new field and interaction free energies
        for(i=0;i<Nr;i++){
            for(ii=0;ii<ChainType;ii++){
                newW[ii][i]=0.0;            //set field update to zero
                for(jj=0;jj<ChainType;jj++){
                    newW[ii][i]+=(chiMatrix[ii][jj]*phi[jj][i]);
                }
                newW[ii][i]+=eta[i];
                
                if (mmb==1){
                    if (ii==0 || ii == 2 || ii==4){
                        newW[ii][i]-=sigma[i];
                    }
                    else if (ii==1 || ii==3){
                        newW[ii][i]+=sigma[i];
                    }
                }
                delW[ii][i]=newW[ii][i]-w[ii][i];
                w[ii][i]+=(gamma*delW[ii][i]-epsilon*delphi[i]);     //update omega field
                deltaW+=fabs(delW[ii][i])*dV(i,dr);
            }
        }
        fE_int=fE(newW,phi,chiMatrix,dr,volume);
        
        //Normalize by box size
        deltaW/=volume;
        
        
        //Update free energy
        fES=Q;
        oldfE=currentfE;
        currentfE=-fES+fE_int;
        deltafE=fabs(currentfE-oldfE);
        
        //Print free energy, difference in free energy, change in omega field to screen
        cout<<xx<<" "<<iter<<" fE:"<<currentfE<< " dfE:"<<currentfE-fE_hom<<" " << deltaW<<" "<<fE_hom<<endl;
        
        
        if (deltafE<precision && deltaW<precision){break;} //Convergence condition
        
    }
    
        
        if (xx==0){
            fE1=currentfE-fE_hom;
            mu[2]=mu2;
        }
        else if(xx==1){
            fE2=currentfE-fE_hom;
            mu3=mu2-(fE2*((mu2-mu1)/(fE2-fE1)));
            mu[2]=mu3;
        }
        else if (xx==2){
            fE3=currentfE-fE_hom;
            if (abs(fE3)>1.0e-6){
                mu1=mu2;
                mu2=mu3;
                mu[2]=mu1;
                xx=-1;
            }
        }
        cout<<"muC: "<<mu[2]<<endl;
        if (xx>2){cout<<"You messed up secant mod"<<endl;}
    }

 //   mu[2]=mu3;
    cout<<"muC: "<<mu[2]<<" done secant mod"<<endl;
    Coord=coordinate;
    
    ofstream outfile;
    outfile.open("./results/mu.dat");
    outfile << mu[2];
    outfile.close();
    
    destroy_1d_double_array(delphi);
    destroy_2d_double_array(delW);
    destroy_2d_double_array(newW);
    
}
