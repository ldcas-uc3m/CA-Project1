using namespace std;

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <cassert>
#include <fstream>


class Object{
    public:
        Object(double x, double y, double z, double mass){
            px = x;
            py = y;
            pz = z;
            m = mass;
            ax = 0;
            ay = 0;
            az = 0;
        }
        double px = 0;
        double py = 0;
        double pz = 0;
        double vx = 0;
        double vy = 0;
        double vz = 0;
        double ax = 0;
        double ay = 0;
        double az = 0;
        double m;
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

    /* ---
    PARAMETERS
    --- */

    // check arguments
    if(argc < 6){
        switch(argc){
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
    }else if(argc > 6){
        cerr << "sim-soa invoked with " << argc << " parameters."
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
    size_enclosure = atoi(argcv[4]);
    time_step = atoi(argcv[5]);

    // chech correct parameters
    if(num_objects <= 0) {
        cerr << "Invalid number of object "<<endl << "sim-aos invoked with " << argc << " parameters."
              << endl << "Arguments: "<< endl << " num_objects: " << num_objects
              << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
              << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
              << " time_step: "<< time_step << endl ;
        return -2;
    }
    if(num_iterations <= 0){
        cerr << "Invalid number of iterations "<<endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl ;
        return -2;
    }
    if( random_seed<= 0){
        cerr << "Invalid seed "<<endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl ;
        return -2;
    }
    if (size_enclosure<= 0){
        cerr << "Invalid box size "<< endl << "sim-aos invoked with " << argc << " parameters."
             << endl << "Arguments: "<< endl << " num_objects: " << num_objects
             << endl << " num_iterations: " << num_iterations << endl << " random_seed: "
             << random_seed << endl << " size_enclosure: " <<size_enclosure << endl
             << " time_step: "<< time_step << endl ;
        return -2;
    }

    // distribution generation
    random_device rd;
    mt19937_64 gen64;  // generate object
    uniform_real_distribution<> dis(0.0, size_enclosure);
    normal_distribution<> d{10e21, 10e15};
    
    gen64.seed(random_seed);  // introduce seed

    // memory alloc
    Object * universe = (Object*)malloc(sizeof(Object) * num_objects);
    
    // init file
    ofstream inFile("init_config.txt", ofstream::out);  // open file
    
    inFile << argcv[4] << " " << argcv[5] << " " << argcv[1];
    inFile << endl;
    
    // populate
    for (int i = 0; i < num_objects; i++){
        universe[i] = Object(dis(gen64), dis(gen64), dis(gen64), d(gen64));
        // write to file
        inFile << universe[i].px << " " << universe[i].py << " " << universe[i].pz 
        << " " << universe[i].vx << " " << universe[i].vy << " " << universe[i].vz 
        << " " << universe[i].m << endl;
    }
    
    inFile.close();

    /* ---
    OUTPUT
    --- */
    ofstream OutFile("final_config.txt");

    OutFile << argcv[4] << " " << argcv[5] << " " << argcv[1] << endl;


    // extra vars
    int curr_objects = num_objects;
    bool *deleted = (bool *)calloc(num_objects, sizeof(bool)); // bytemap of objects -> if true, object is deleted

    /* ---
    KERNEL
    --- */
 
    for(int iteration = 0; iteration < num_iterations; iteration++){
        for(int i = 0; i < num_objects; i++){
            if(deleted[i]) continue;
            Object *a = &universe[i];

            for(int j = i + 1; j < num_objects; j++){
                if(deleted[j]) continue;
                
                Object *b = &universe[j];

                /* ---
                FORCE COMPUTATION
                --- */

                double fax = 0;
                double fay = 0;
                double faz = 0;

                // distance
                double dx = b->px - a->px;
                double dy = b->py - a->py;
                double dz = b->pz - a->pz;
                double distance = sqrt(dx*dx + dy*dy + dz*dz);

                if(distance <= COL_DISTANCE){
                    /* ---
                    OBJECT COLLISION
                    --- */

                    // merge objects into a
                    a->vx = a->vx + b->vx;
                    a->vy = a->vy + b->vy;
                    a->vz = a->vz + b->vz;
                    a->m = a->m + b->m;

                    // del b
                    curr_objects--;
                    deleted[j] = true;

                    // force between a & b is 0
                } else{
                
                    fax = (g * a->m * b->m * dx) / abs(distance*distance*distance);
                    fay = (g * a->m * b->m * dy) / abs(distance*distance*distance);
                    faz = (g * a->m * b->m * dz) / abs(distance*distance*distance);

                    double fbx = -fax;
                    double fby = -fay;
                    double fbz = -faz;

                    // b acceleration
                    b->ax -= fbx/b->m;
                    b->ay -= fby/b->m;
                    b->az -= fbz/b->m;
                }

                // a acceleration
                a->ax += fax/a->m;
                a->ay += fay/a->m;
                a->az += faz/a->m;
            }

            /* ---
            UPDATE POSITION
            --- */
            // velocity calculation
            double vx = a->vx + a->ax * time_step;
            double vy = a->vy + a->ay * time_step;
            double vz = a->vz + a->az * time_step;

            a->vx = vx;
            a->vy = vy;
            a->vz = vz;
            
            // position calculation
            a->px += vx * time_step;
            a->py += vy * time_step;
            a->py += vz * time_step;            
            

            /* ---
            REBOUND EFFECT
            --- */

            if(a->px <= 0){
                a->px = 0;
                a->vx = - a->vx;
            } else if(a->px >= size_enclosure){
                a->px = size_enclosure;
                a->vx = - a->vx;
            }

            if(a->py <= 0){
                a->py = 0;
                a->vy = - a->vy;
            } else if(a->py >= size_enclosure){
                a->py = size_enclosure;
                a->vy = - a->vy;
            }

            if(a->pz <= 0){
                a->pz = 0;
                a->vz = - a->vz;
            } else if(a->pz >= size_enclosure){
                a->pz = size_enclosure;
                a->vz = - a->vz;
            }

            // print to output
            if((iteration = num_iterations - 1) || curr_objects == 1){  // final positions

            OutFile << universe[i].px << " " << universe[i].py << " " << universe[i].pz 
            << " " << universe[i].vx << " " << universe[i].vy << " " << universe[i].vz 
            << " " << universe[i].m << endl;
            }
        }
    }
    OutFile.close();
    return 0;
}
