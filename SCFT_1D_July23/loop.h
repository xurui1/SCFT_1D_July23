void loop(double **qA2, double **qA3,int *Ns,double **phi,double volume,double dr){
    
    int i;
    double vL; //percentage of looping chains
    double Q_ABA;
    double *v_L;
    
    v_L=create_1d_double_array(Nr,"v_L");
    
    for (i=0;i<Nr;i++){
        Q_ABA+=qA3[i][Ns[0]];
    }
    Q_ABA/=2.0;
    
    vL=0;
    v_L[0]=qA2[0][Ns[0]]*(qA3[0][0])/Q_ABA;
    v_L[Nr-1]=qA2[Nr-1][Ns[0]]*(qA3[Nr-1][0])/Q_ABA;
    
    for (i=1;i<Nr-1;i++){
        v_L[i]+=qA2[i][Ns[0]]*(qA3[i][0])/Q_ABA;
    }

   /* ofstream file1("./results/Loop.dat");
    for (i=0;i<Nr;i++){
        file1 << i*dr <<" "<<v_L[i]<<endl;
    }
    file1.close();
    cout <<"loop fraction: "<<vL<<endl;
    */
    
}