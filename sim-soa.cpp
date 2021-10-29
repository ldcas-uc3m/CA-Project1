using namespace std;

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <cassert>
#include <fstream>



class Universe{
    public:
        Universe(int num_objects, int size_enclosure){
            objects = num_objects;
            size = size_enclosure;
            px = (double *)malloc(sizeof(double) * num_objects);
            py = (double *)malloc(sizeof(double) * num_objects);
            pz = (double *)malloc(sizeof(double) * num_objects);
            vx = (double *)malloc(sizeof(double) * num_objects);
            vy = (double *)malloc(sizeof(double) * num_objects);
            vz = (double *)malloc(sizeof(double) * num_objects);
            m = (double *)malloc(sizeof(double) * num_objects);
            fx = (double *)malloc(sizeof(double) * num_objects);
            fy = (double *)malloc(sizeof(double) * num_objects);
            fz = (double *)malloc(sizeof(double) * num_objects);
        }
        int objects;
        int size;
        double * px;
        double * py;
        double * pz;
        double * vx;
        double * vy;
        double * vz;
        double * m;
        double * fx;
        double * fy;
        double * fz;
};


// constants
const double g = 6.674e-11;
const double COL_DISTANCE = 1;  // minimum colision distance

// Global variables
int num_objects;
int num_iterations;
int random_seed;
double size_enclosure;
double time_step;


int main(int argc, const char ** argcv){
    // check arguments
        if(argc<6){
            switch (argc) {
                case 5:
                    cerr  << "sim-soa invoked with " << argc << " parameters."
                         << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                         << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                         << argcv[3] << endl << " size_enclosure: " << argcv[4] << endl
                         << " time_step: ?"<< endl ;
                    break;
                case 4:
                    cerr << "sim-soa invoked with " << argc << " parameters."
                         << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                         << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                         << argcv[3] << endl << " size_enclosure: ?"<< endl
                         << " time_step: ?"<< endl ;
                    break;
                case 3:
                    cerr << "sim-soa invoked with " << argc << " parameters."
                         << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                         << endl << " num_iterations: " << argcv[2] << endl << " random_seed: ?"
                         << endl << " size_enclosure: ?"<< endl
                         << " time_step: ?"<< endl ;
                    break;
                case 2:
                    cerr << "sim-soa invoked with " << argc << " parameters."
                         << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                         << endl << " num_iterations: ?" << endl << " random_seed: ?"
                         << endl << " size_enclosure: ?"<< endl
                         << " time_step: ?"<< endl ;
                    break;
                case 1:
                    cerr << "sim-soa invoked with " << argc << " parameters."
                         << endl << "Arguments: "<< endl << " num_objects: ?"
                         << endl << " num_iterations: ?" << endl << " random_seed: ?"
                         << endl << " size_enclosure: ?"<< endl
                         << " time_step: ?"<< endl ;
                    break;
            }
            return -1;
    }
        else if(argc > 6){
            cerr  << "sim-soa invoked with " << argc << " parameters."
                  << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                  << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                  << argcv[3] << endl << " size_enclosure: " << argcv[4] << endl
                  << " time_step: "<< argcv[5] << endl ;
            return -1;
        }
    // parameters init & casting
    num_objects = atoi(argcv[1]);
    num_iterations = atoi(argcv[2]);
    random_seed = atoi(argcv[3]);
    size_enclosure = atof(argcv[4]);
    time_step = atof(argcv[5]);

    if(num_objects <= 0) {
        cerr << "Invalid number of object "<<endl << "sim-soa invoked with " << argc << " parameters."
              << endl << "Arguments: "<< endl << " num_objects: " << num_objects
              << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
              << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
              << " time_step: "<< time_step << endl ;
        return -2;
    }
    if(num_iterations <= 0){
        cerr << "Invalid number of iterations "<<endl << "sim-soa invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl ;
        return -2;
    }
    if( random_seed<= 0){
        cerr << "Invalid seed "<<endl << "sim-soa invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl ;
        return -2;
    }
    if (size_enclosure<= 0 || size_enclosure < num_objects){
        cerr << "Invalid box size "<< endl << "sim-soa invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl ;
        return -2;
    }

    // distribution generation
    std::random_device rd;
    mt19937_64 gen64;  // generate object
    uniform_real_distribution<> dis(0.0, size_enclosure);
    normal_distribution<> d{1e21, 1e15};

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
        universe.fx[i] = 0;
        universe.fy[i] = 0;
        universe.fz[i] = 0;

        // write to file
        MyFile << universe.px[i] << " " << universe.py[i] << " " << universe.pz[i]
        << " " << universe.vx[i] << " " << universe.vy[i] << " " << universe.vz[i] 
        << " " << universe.m[i];
        MyFile << endl;
    }
    
    MyFile.close();
    
    int curr_objects = num_objects;

    bool *deleted = (bool *)calloc(num_objects, sizeof(bool));  // bytemap of objects -> if true, object is deleted

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
                // distance
                double dx = universe.px[j] - universe.px[i];
                double dy = universe.py[j] - universe.py[i];
                double dz = universe.pz[j] - universe.pz[i];
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if(distance <= COL_DISTANCE){
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
                
               		double dfx = (g * universe.m[i] * universe.m[j] * dx) / (distance*distance*distance);
                    double dfy = (g * universe.m[i] * universe.m[j] * dy) / (distance*distance*distance);
                    double dfz = (g * universe.m[i] * universe.m[j] * dz) / (distance*distance*distance);

                    // a forces
                    universe.fx[i] += dfx;
                    universe.fy[i] += dfy;
                    universe.fz[i] += dfz;

                    // b forces
                    universe.fx[j] -= dfx;
                    universe.fy[j] -= dfy;
                    universe.fz[j] -= dfz;
                }
            }

            /* ---
            UPDATE POSITION
            --- */
            // acceleration calculation
            double ax = universe.fx[i]/universe.m[i];
            double ay = universe.fy[i]/universe.m[i];
            double az = universe.fz[i]/universe.m[i];

            //reset force 
            universe.fx[i] = 0;
            universe.fy[i] = 0;
            universe.fz[i] = 0;

            // velocity calculation
            universe.vx[i] += ax * time_step;
            universe.vy[i] += ay * time_step;
            universe.vz[i] += az * time_step;

            // position calculation
            universe.px[i] += universe.vx[i] * time_step;
            universe.py[i] += universe.vy[i] * time_step;
            universe.pz[i] += universe.vz[i] * time_step;  

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
