void omega(double **w){
    
    int i,x;
    
    //This is for a sinusoidal initial omega field. Adjust 'initial' in the parameters header.
    //Currently will give a lamellar structure in z-direction
    if (initial==0){
        for(i=0;i<Nr;i++){
            w[0][i]=-5.0*cos(2.0*Pi*(double)i/(double)Nr)+(double)((rand() %100)/50.0)-1.0; //A1
            w[1][i]=5.0*cos(2.0*Pi*(double)i/(double)Nr)+(double)((rand() %100)/50.0)-1.0;  //B1
            w[2][i]=-5.0*cos(2.0*Pi*(double)i/(double)Nr)+(double)((rand() %100)/50.0)-1.0; //A2
            w[3][i]=5.0*cos(2.0*Pi*(double)i/(double)Nr)+(double)((rand() %100)/50.0)-1.0;  //B2
            w[4][i]=-5.0*cos(2.0*Pi*(double)i/(double)Nr)+(double)((rand() %100)/50.0)-1.0; //A3
            w[5][i]=10.0+(double)((rand() %100)/50.0-1.0);                                  //C
        }
    }
    
    
    //This is for a bilayer conformation with 50 points
    else if (initial==1){
        ifstream Init1;
        Init1.open("bilayer_N50.dat");
        for (i=0;i<Nr;i++){
                Init1 >> x >> w[0][i]>>w[1][i]>>w[5][i];
        }
        
        for (i=0;i<Nr;i++){
            w[2][i]=w[0][i];
            w[3][i]=w[1][i];
            w[4][i]=w[4][i];
        }
        
        Init1.close();
    }
    //This is for a bilayer conformation with 100 points
    else if (initial==2){
        ifstream Init2;
        Init2.open("bilayer_N100.dat");
        for (i=0;i<Nr;i++){
            Init2 >> x >> w[0][i]>>w[1][i]>>w[5][i];
        }
        
        for (i=0;i<Nr;i++){
            w[2][i]=w[0][i];
            w[3][i]=w[1][i];
            w[4][i]=w[4][i];
        }
        
        Init2.close();
    }
    
    for (i=0;i<Nr;i++){
                std::cout << w[0][i]<<" "<<w[1][i]<<" "<<w[2][i]<<endl;
    }
    
}