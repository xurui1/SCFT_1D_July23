#include <cstring>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

void *smalloc(int n, const char *name){
    
    if (n == 0) return NULL;
  		void *ptr = malloc(n);
    if (ptr == NULL) {
        char str[128];
        std::cout<<"Failed to allocate %d bytes for array %s"<<n<<name<<std::endl;
    }
    return ptr;
}



void sfree(void *ptr){
    
    if (ptr == NULL) return;
    free(ptr);
}

double *create_1d_double_array(int n1, const char *name){
    
    double *array = (double *) smalloc(n1*sizeof(double),name);
    return array;
}




void destroy_1d_double_array(double *array){
    
    if (array == NULL) return;
    sfree(array);
}

int *create_1d_integer_array(int n1, const char *name){
    
    int *array = (int *) smalloc(n1*sizeof(int),name);
    return array;
}

void destroy_1d_integer_array(int *array){
    
    if (array == NULL) return;
    sfree(array);
}


double **create_2d_double_array(int n1, int n2, const char *name){
    
    double *data = (double *) smalloc(n1*n2*sizeof(double),name);
    double **array = (double **) smalloc(n1*sizeof(double *),name);
    
    int n = 0;
    for (int i = 0; i < n1; i++){
        array[i] = &data[n];
        n += n2;
    }
    
    return array;
}



void destroy_2d_double_array(double **array){
    
    if (array == NULL) return;
    sfree(array[0]);
    sfree(array);
}



double ***create_3d_double_array(int n1, int n2, int n3, const char *name){
    
    int i,j;
    
    double *data = (double *) smalloc(n1*n2*n3*sizeof(double),name);
    double **plane = (double **) smalloc(n1*n2*sizeof(double *),name);
    double ***array = (double ***) smalloc(n1*sizeof(double **),name);
    
    int n = 0;
    for (i = 0; i < n1; i++) {
        array[i] = &plane[i*n2];
        for (j = 0; j < n2; j++) {
            plane[i*n2+j] = &data[n];
            n += n3;
        }
    }
    
    return array;
}



void destroy_3d_double_array(double ***array){
    
    if (array == NULL) return;
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
}



double ****create_4d_double_array(int n1, int n2, int n3, int n4, const char *name){
    
    int i,j,k;
    
    double *data = (double *) smalloc(n1*n2*n3*n4*sizeof(double),name);
    double **cube = (double **) smalloc(n1*n2*n3*sizeof(double *),name);
    double ***plane = (double ***) smalloc(n1*n2*sizeof(double **),name);
    double ****array = (double ****) smalloc(n1*sizeof(double ***),name);
    
    int n = 0;
    for (i = 0; i < n1; i++) {
        array[i] = &plane[i*n2];
        for (j = 0; j < n2; j++) {
            plane[i*n2+j] = &cube[i*n2*n3+j*n3];
            for (k = 0; k < n3; k++) {
                cube[i*n2*n3+j*n3+k] = &data[n];
                n += n4;
            }
        }
    }
    return array;
}



void destroy_4d_double_array(double ****array)
{
    if (array == NULL) return;
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
}
