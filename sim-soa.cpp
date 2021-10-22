using namespace std;

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <cassert>
#include <fstream>
#include <array>


class Force{
    public:
        Force(double fx, double fy, double fz){
            x = fx;
            y = fy;
            z = fz;
        }
        double x;
        double y;
        double z;
};


class Universe{
    public:
        Universe(int num_objects, int size_enclosure){
            objects = num_objects;
            size = size_enclosure;
        }
        int objects;
        int size;
        double * px = (double *)malloc(sizeof(double)*objects);
        double * py = (double *)malloc(sizeof(double)*objects);
        double * pz = (double *)malloc(sizeof(double)*objects);
        double * vx = (double *)malloc(sizeof(double)*objects);
        double * vy = (double *)malloc(sizeof(double)*objects);
        double * vz = (double *)malloc(sizeof(double)*objects);
        double * m = (double *)malloc(sizeof(double)*objects);
        double * ax = (double *)malloc(sizeof(double)*objects);
        double * ay = (double *)malloc(sizeof(double)*objects);
        double * az = (double *)malloc(sizeof(double)*objects);
};


// constants
const double g = 6.674 * pow(10, -11);
const double COL_DISTANCE = 1;  // minimum colision distance

// Global variables
int num_objects;
int num_iterations;
int random_seed;
int size_enclosure;
int time_step;


int main(int argc, const char ** argcv){

    // check arguments
    if (argc != 6 || num_objects < 0 || num_iterations < 0 
        || random_seed < 0 || size_enclosure < 0 || time_step < 0){
        cout << "sim-aos invoked with " << argc << "parameters." 
        << endl << "Arguments: "<< endl << " num_objects: " << argcv[1] 
        << endl << " num_iterations: " << argcv[2] << endl << " random_seed: " 
        << argcv[3] << endl << " size_enclosure: " << argcv[4] << endl 
        << " time_step: " << argcv[5] << endl ;
    }

    // parameters init & casting
    num_objects = atoi(argcv[1]);
    num_iterations = atoi(argcv[2]);
    random_seed = atoi(argcv[3]);
    size_enclosure = atoi(argcv[4]);
    time_step = atoi(argcv[5]);

    // distribution generation
    std::random_device rd;
    mt19937_64 gen64;  // generate object
    uniform_real_distribution<> dis(0.0, size_enclosure);
    normal_distribution<double> d{pow(10, 21),pow(10, 15)};

    gen64.seed(random_seed);  // introduce seed

    // big bang
    Universe universe(num_objects, size_enclosure);
    
    ofstream MyFile("init_config.txt", ofstream::out);  // open file
    
    MyFile << argcv[4] << " " << argcv[5] << " " << argcv[1];
    MyFile << endl;
    
    // populate
    for (int i = 0; i < num_objects; i++){
        universe.px[i] = dis(gen64);
        universe.py[i] = dis(gen64);
        universe.pz[i] = dis(gen64);
        universe.vx[i] = 0;
        universe.vy[i] = 0;
        universe.vz[i] = 0;
        universe.m[i] = d(gen64);
        universe.ax[i] = 0;
        universe.ay[i] = 0;
        universe.az[i] = 0;

        // write to file
        MyFile << universe.px[i] << " " << universe.py[i] << " " << universe.pz[i]
        << " " << universe.vx[i] << " " << universe.vy[i] << " " << universe.vz[i] 
        << " " << universe.m[i];
        MyFile << endl;
    }
    
    MyFile.close();
    
    int curr_objects = num_objects;

    bool deleted[num_objects];  // bytemap of objects -> if true, object is deleted

    /* ---
    KERNEL
    --- */
    
    for(int iteration = 0; iteration < num_iterations; iteration++){
        if(curr_objects == 0) break;
        for(int i = 0; i < num_objects; i++){
            if(deleted[i]) continue;
            for(int j = i + 1; j < num_objects; j++){
                if(deleted[j]) continue;
                
                /* ---
                FORCE COMPUTATION
                --- */
                Force fa(0, 0, 0);
                
                // distance
                double dx = universe.px[j] - universe.px[i];
                double dy = universe.py[j] - universe.py[i];
                double dz = universe.pz[j] - universe.pz[i];
                double distance = sqrt(dx*dx + dy*dy + dz*dz);

                if(distance < COL_DISTANCE){
                    /* ---
                    OBJECT COLLISION
                    --- */

                    // merge objects into a (i)
                    universe.m[i] = universe.m[i] + universe.m[j];
                    universe.vx[i] = universe.vx[i] + universe.vx[j];
                    universe.vy[i] = universe.vy[i] + universe.vy[j];
                    universe.vz[i] = universe.vz[i] + universe.vz[j];

                    // delete b (j)
                    curr_objects--;
                    deleted[j] = true;

                    // force between a & b is 0
                } else{
                
                fa.x = (g * universe.m[i] * universe.m[j] * dy) / abs(pow(dy, 3));
                fa.y = (g * universe.m[i] * universe.m[j] * dy) / abs(pow(dy, 3));
                fa.z = (g * universe.m[i] * universe.m[j] * dz) / abs(pow(dz, 3));

                Force fb(- fa.x, -fa.y, -fa.z);

                // b acceleration
                universe.ax[j] -= fb.x/universe.m[j];
                universe.ay[j] -= fb.y/universe.m[j];
                universe.az[j] -= fb.z/universe.m[j];
                }

                // a acceleration
                universe.ax[i] += fa.x/universe.m[i];
                universe.ay[i] += fa.y/universe.m[i];
                universe.az[i] += fa.z/universe.m[i];

            }
            
            /* ---
            UPDATE POSITION
            --- */
            // velocity calculation
            double vx = universe.vx[i] + universe.ax[i] * time_step;
            double vy = universe.vy[i] + universe.ay[i] * time_step;
            double vz = universe.vz[i] + universe.az[i] * time_step;

            universe.vx[i] = vx;
            universe.vy[i] = vy;
            universe.vz[i] = vz;

            // position calculation
            universe.px[i] += vx * time_step;
            universe.py[i] += vz * time_step;
            universe.py[i] += vz * time_step;

            /* ---
            REBOUND EFFECT
            --- */

            if(universe.px[i] <= 0){
                universe.px[i] = 0;
                universe.vx[i] = - universe.vx[i];
            } else if(universe.px[i] >= size_enclosure){
                universe.px[i] = size_enclosure;
                universe.vx[i] = - universe.vx[i];
            }

            if(universe.py[i] <= 0){
                universe.py[i] = 0;
                universe.vy[i] = - universe.vy[i];
            } else if(universe.py[i] >= size_enclosure){
                universe.py[i] = size_enclosure;
                universe.vy[i] = - universe.vy[i];
            }

            if(universe.pz[i] <= 0){
                universe.pz[i] = 0;
                universe.vz[i] = - universe.vz[i];
            } else if(universe.pz[i] >= size_enclosure){
                universe.pz[i] = size_enclosure;
                universe.vz[i] = - universe.vz[i];
            }

            //cout << "iteration " << iteration << ", object " << i << " | " << universe.px[i] << " " << universe.py[i] << " " << universe.pz[i] << " | " << universe.vx[i] << " " << universe.vy[i] << " " << universe.vz[i] << " | " << universe.m[i] << endl;
        }
    }

    /*
    OUTPUT
    */
    ofstream OutFile("final_config.txt");

    OutFile << argcv[4] << " " << argcv[5] << " " << argcv[1];
    OutFile << endl;

    for (int i = 0; i < num_objects; i++){
        if(deleted[i]) continue;
        // write to file
        OutFile << universe.px[i] << " " << universe.py[i] << " " << universe.pz[i] 
        << " " << universe.vx[i] << " " << universe.vy[i] << " " << universe.vz[i] 
        << " " << universe.m[i];
        OutFile << endl;
    }

    OutFile.close();
}