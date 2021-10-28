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
            ax = (double *)malloc(sizeof(double) * num_objects);
            ay = (double *)malloc(sizeof(double) * num_objects);
            az = (double *)malloc(sizeof(double) * num_objects);
        }
        double * px;
        double * py;
        double * pz;
        double * vx;
        double * vy;
        double * vz;
        double * m;
        double * ax;
        double * ay;
        double * az;
        int objects;
        int size;

};


// constants
const double g = 6.674e-11;
const double COL_DISTANCE = 1;  // minimum colision distance

// Global variables
int num_objects;
int num_iterations;
int random_seed;
int size_enclosure;
int time_step;


int main(int argc, const char ** argcv){

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

    // big bang
    Universe universe(num_objects, size_enclosure);
    
    ofstream inFile("init_config.txt", ofstream::out);  // open file
    inFile << argcv[4] << " " << argcv[5] << " " << argcv[1] << endl;
    
    for(int i = 0; i < num_objects; i++){
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

        inFile << universe.px[i] << " " << universe.py[i] << " " << universe.pz[i]
        << " " << universe.vx[i] << " " << universe.vy[i] << " " << universe.vz[i] 
        << " " << universe.m[i] << endl;
    }

    inFile.close();
    
    /* ---
    OUTPUT
    --- */
    ofstream outFile("final_config.txt");

    outFile << argcv[4] << " " << argcv[5] << " " << argcv[1] << endl;

    // extra vars
    int curr_objects = num_objects;
    bool *deleted = (bool *)calloc(num_objects, sizeof(bool));  // bytemap of objects -> if true, object is deleted


    /* ---
    KERNEL
    --- */
    
    for(int iteration = 0; iteration < num_iterations; iteration++){
        for(int i = 0; i < num_objects; i++){
            if(deleted[i]) continue;

            for(int j = i + 1; j < num_objects; j++){
                if(deleted[j]) continue;
                
                /* ---
                FORCE COMPUTATION
                --- */
                double fax = 0;
                double fay = 0;
                double faz = 0;
                
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
                    universe.vx[i] = universe.vx[i] + universe.vx[j];
                    universe.vy[i] = universe.vy[i] + universe.vy[j];
                    universe.vz[i] = universe.vz[i] + universe.vz[j];
                    universe.m[i] = universe.m[i] + universe.m[j];

                    // delete b (j)
                    curr_objects--;
                    deleted[j] = true;

                    // force between a & b is 0
                } else{
                
                    fax = (g * universe.m[i] * universe.m[j] * dx) / abs(distance*distance*distance);
                    fay = (g * universe.m[i] * universe.m[j] * dy) / abs(distance*distance*distance);
                    faz = (g * universe.m[i] * universe.m[j] * dz) / abs(distance*distance*distance);

                    double fbx = -fax;
                    double fby = -fay;
                    double fbz = -faz;

                    // b acceleration
                    universe.ax[j] -= fbx/universe.m[j];
                    universe.ay[j] -= fby/universe.m[j];
                    universe.az[j] -= fbz/universe.m[j];
                }

                // a acceleration
                universe.ax[i] += fax/universe.m[i];
                universe.ay[i] += fay/universe.m[i];
                universe.az[i] += faz/universe.m[i];

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
            universe.py[i] += vy * time_step;
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

            // print to output

            if((iteration == num_iterations - 1) ||  curr_objects == 1) {  // final positions

            outFile << universe.px[i] << " " << universe.py[i] << " " << universe.pz[i] 
            << " " << universe.vx[i] << " " << universe.vy[i] << " " << universe.vz[i] 
            << " " << universe.m[i] << endl;
            }


            //cout << "iteration: " << iteration << ", element: " << i << " | " << universe.px[i] << " " << universe.py[i] << " " << universe.pz[i] << " | " << universe.vx[i] << " " << universe.vy[i] << " " << universe.vz[i] << endl;
        }
    }

    outFile.close();
    return 0;
}
