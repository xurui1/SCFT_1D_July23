#include <stdio.h>     //Include the standard input/output libraries
#include <iostream>  //Cout and Cin etc.
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>    //Include standard fucntion libraries
#include <math.h>      //Use the math function libraries
#include <time.h>      //Call system time libraries to define the integer seed for random numbers
#include "smemory.h"  //Use my custom memory handling class
//#include "mpi.h"     //Use this for MPI parallel implimentation later

using namespace std;

#define Nr 100
#define ChainType 6
#define Pi 3.14159


double kappa;
double r_0;
int initial;
int Coord;
