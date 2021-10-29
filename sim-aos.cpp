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
            vx = 0;
            vy = 0;
            vz = 0;
            m = mass;
            fx = 0;
            fy = 0;
            fz = 0;
        }
		double m;
        double px;
        double vx;
        double fx;
        double py;
        double vy;
        double fy;
        double pz;
        double vz;
        double fz;
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
    if(argc<6){
        switch (argc) {
            case 5:
                cerr  << "sim-aos invoked with " << argc << " parameters."
                      << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                      << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                      << argcv[3] << endl << " size_enclosure: " << argcv[4] << endl
                      << " time_step: ?"<< endl ;
                break;
            case 4:
                cerr << "sim-aos invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: " << argcv[2] << endl << " random_seed: "
                     << argcv[3] << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl ;
                break;
            case 3:
                cerr << "sim-aos invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: " << argcv[2] << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl ;
                break;
            case 2:
                cerr << "sim-aos invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: " << argcv[1]
                     << endl << " num_iterations: ?" << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl ;
                break;
            case 1:
                cerr << "sim-aos invoked with " << argc << " parameters."
                     << endl << "Arguments: "<< endl << " num_objects: ?"
                     << endl << " num_iterations: ?" << endl << " random_seed: ?"
                     << endl << " size_enclosure: ?"<< endl
                     << " time_step: ?"<< endl ;
                break;
        }
        return -1;
    }else if(argc > 6){
        cerr  << "sim-aos invoked with " << argc << " parameters."
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
    if (size_enclosure<= 0 || size_enclosure < num_objects ){
        cerr << "Invalid box size "<< endl << "sim-aos invoked with " << argc << " parameters."
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

    // memory alloc
    Object * universe = (Object*)malloc(sizeof(Object) * num_objects);
    
    ofstream MyFile("init_config.txt", ofstream::out);  // open file
    
    MyFile << argcv[4] << " " << argcv[5] << " " << argcv[1];
    MyFile << endl;
    
    // populate
    for (int i = 0; i < num_objects; i++){
        double x = dis(gen64);
        double y = dis(gen64);
        double z = dis(gen64);
        double m = d(gen64);

        universe[i] = Object(x, y, z, m);
        // write to file
        MyFile << universe[i].px << " " << universe[i].py << " " << universe[i].pz 
        << " " << universe[i].vx << " " << universe[i].vy << " " << universe[i].vz 
        << " " << universe[i].m << endl;
    }
    
    MyFile.close();

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

                if(distance <= COL_DISTANCE){
                    /* ---
                    OBJECT COLLISION
                    --- */

                    // merge objects into a
                    a->m = a->m + b->m;
                    a->vx = a->vx + b->vx;
                    a->vy = a->vy + b->vy;
                    a->vz = a->vz + b->vz;

                    // del b
                    deleted[j] = true;

                    // force between a & b is 0
                } else{

                    // a forces
                    a->fx += (g * a->m * b->m * (b->px - a->px)) / ((std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz))));
                    a->fy += (g * a->m * b->m * (b->py - a->py)) / ((std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz))));
                    a->fz += (g * a->m * b->m * (b->pz - a->pz)) / ((std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz))));

                    // b forces
                    b->fx -= (g * a->m * b->m * (b->px - a->px)) / ((std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz))));
                    b->fy -= (g * a->m * b->m * (b->py - a->py)) / ((std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz))));
                    b->fz -= (g * a->m * b->m * (b->pz - a->pz)) / ((std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz)))*(std::sqrt((b->px - a->px)*(b->px - a->px) + (b->py - a->py)*(b->py - a->py) + (b->pz - a->pz)*(b->pz - a->pz))));
                }
            }

            /* ---
            UPDATE POSITION
            --- */

            // velocity calculation
            a->vx += (a->fx/a->m) * time_step;
            a->vy += (a->fy/a->m) * time_step;
            a->vz += (a->fz/a->m) * time_step;
			
			// reset force 
            a->fx = 0;
            a->fy = 0;
            a->fz = 0;
            
            // position calculation
            a->px += a->vx * time_step;
            a->py += a->vy * time_step;
            a->pz += a->vz * time_step;  

            /* ---
            REBOUND EFFECT
            --- */

            if(a->px <= 0){
                a->px = 0;
                a->vx = - a->vx;
            }
	    	if(a->px >= size_enclosure){
                a->px = size_enclosure;
                a->vx = - a->vx;
            }

            if(a->py <= 0){
                a->py = 0;
                a->vy = - a->vy;
            }
			if(a->py >= size_enclosure){
                a->py = size_enclosure;
                a->vy = - a->vy;
            }

            if(a->pz <= 0){
                a->pz = 0;
                a->vz = - a->vz;
            }
			if(a->pz >= size_enclosure){
                a->pz = size_enclosure;
                a->vz = - a->vz;
            }     
        }
    }
    /*
    OUTPUT
    */
    ofstream OutFile("final_config.txt");

    OutFile << argcv[4] << " " << argcv[5] << " " << argcv[1] << endl;

    for (int i = 0; i < num_objects; i++){
        if(deleted[i]) continue;
        // write to file
        OutFile << universe[i].px << " " << universe[i].py << " " << universe[i].pz 
        << " " << universe[i].vx << " " << universe[i].vy << " " << universe[i].vz 
        << " " << universe[i].m << endl;
    }

    OutFile.close();
}
