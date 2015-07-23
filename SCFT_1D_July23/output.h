void output(double dr, double **phi,double **w){
    
    int i;
    
    ofstream outputFile1("./results/phi.dat");
    
    for(i=0;i<Nr;i++){
        outputFile1 <<i*dr<<" "<<phi[0][i]<<" "<<phi[1][i]<<" "<<phi[2][i]<<" "<<phi[3][i]<<" "<<phi[4][i]<<" "<<phi[5][i]<<std::endl;
    }
    
    outputFile1.close();
    
    ofstream outputFile2("./results/omega.dat");
    
    for(i=0;i<Nr;i++){
        outputFile2 <<i*dr<<" "<<w[0][i]<<" "<<w[1][i]<<" "<<w[2][i]<<" "<<w[3][i]<<" "<<w[4][i]<<" "<<w[5][i]<<std::endl;
    }
    outputFile2.close();

    
}