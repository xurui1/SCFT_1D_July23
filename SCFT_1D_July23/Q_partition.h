
double q_partition(double **qB1,double **qA3,double **qC, double dr, int *Ns, double *mu, double volume){
    
    double Q,Q_AB,Q_C,Q_ABA;
    int i;
 
    Q=0.0;
    Q_AB=0.0;
    Q_C=0.0;
    Q_ABA=0.0;
    
    //ends
    Q_AB+=0.5*qB1[0][Ns[1]]*dV(0,dr);
    Q_ABA+=0.5*qA3[0][Ns[0]]*dV(0,dr);
    Q_C+=0.5*qC[0][Ns[2]]*dV(0,dr);
    
    Q_AB+=0.5*qB1[(int)Nr-1][Ns[1]]*dV(Nr-1,dr);
    Q_ABA+=0.5*qA3[(int)Nr-1][Ns[0]]*dV(Nr-1,dr);
    Q_C+=0.5*qC[(int)Nr-1][Ns[2]]*dV(Nr-1,dr);
    
    //Middle
    for(i=1;i<(int)Nr-1;i++){
        Q_AB+=qB1[i][Ns[1]]*dV(i,dr);
        Q_ABA+=qA3[i][Ns[0]]*dV(i,dr);
        Q_C+=qC[i][Ns[2]]*dV(i,dr);
    }
    
    Q_AB=exp(mu[0])*Q_AB;
    Q_ABA=exp(mu[1]*2.0)*Q_ABA/2.0;
    Q_C=(exp(mu[2]*kappa)*Q_C)/kappa;
    //I'm adding the three single chain partition functions together for the return function
    Q=Q_AB+Q_C+Q_ABA;
    // Normalizing with respect to box volume
    Q/=volume;
    
    
    return Q;
}