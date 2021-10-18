//
// Created by ivan on 9/10/21.
//
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
const double g = 6.674 * pow(10, -11);
const double COL_DISTANCE = 0;  // minimum colision distance

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
    random_device rd;
    mt19937_64 gen64;  // generate object
    uniform_real_distribution<> dis(0.0, size_enclosure);
    normal_distribution<double> d{pow(10, 21),pow(10, 15)};
    
    gen64.seed(random_seed);  // introduce seed

    // memory alloc
    Object * universe = (Object*)malloc(sizeof(Object) * num_objects);
    
    ofstream MyFile("init_config.txt", ofstream::out);  // open file
    
    MyFile << argcv[4] << " " << argcv[5] << " " << argcv[1];
    MyFile << endl;
    
    // populate
    for (int i = 0; i < num_objects; i++){
        universe[i] = Object(dis(gen64), dis(gen64), dis(gen64), d(gen64));
        // write to file
        MyFile << universe[i].px << " " << universe[i].py << " " << universe[i].pz 
        << " " << universe[i].vx << " " << universe[i].vy << " " << universe[i].vz 
        << " " << universe[i].m;
        MyFile << endl;
    }
    
    MyFile.close();

    int curr_objects = num_objects;

    bool Bmap[num_objects] = {true};  // bytemap of objects

    /* ---
    KERNEL
    --- */

    for(int iteration; iteration < num_iterations; iteration++){
        if(curr_objects == 0) break;
        Object a(0,0,0,0);
        for(int i = 0; i < num_objects; i++){
            if(not Bmap[i]) continue;
            a = universe[i];
            Object b(0,0,0,0);
            for(int j = i + 1; j < num_objects; j++){
                if(not Bmap[i]) continue;
                b = universe[j];

                /* ---
                FORCE COMPUTATION
                --- */
                Force fa(0, 0, 0);

                // distance
                double dx = b.px - a.px;
                double dy = b.py - a.py;
                double dz = b.pz - a.pz;
                double distance = sqrt(dx*dx + dy*dy + dz*dz);

                if(distance < COL_DISTANCE){
                    /* ---
                    OBJECT COLLISION
                    --- */

                    // merge objects into a
                    a.m = a.m + b.m;
                    a.vx = a.vx + b.vx;
                    a.vy = a.vy + b.vy;
                    a.vz = a.vz + b.vz;

                    // del object
                    delete &(universe[j]);
                    curr_objects--;
                    Bmap[j] = false;

                    // force between a & b is 0
                } else{
                
                fa.x = (g * a.m * b.m * dy) / abs(pow(dy, 3));
                fa.y = (g * a.m * b.m * dy) / abs(pow(dy, 3));
                fa.z = (g * a.m * b.m * dz) / abs(pow(dz, 3));

                Force fb(- fa.x, -fa.y, -fa.z);

                // b acceleration
                b.ax -= fb.x/a.m;
                b.ay -= fb.y/a.m;
                b.az -= fb.z/a.m;
                }

                // a acceleration
                a.ax += fa.x/a.m;
                a.ay += fa.y/a.m;
                a.az += fa.z/a.m;

            }
            /* ---
            UPDATE POSITION
            --- */
            // velocity calculation
            double vx = a.vx + a.ax * time_step;
            double vy = a.vy + a.ay * time_step;
            double vz = a.vz + a.az * time_step;

            a.vx = vx;
            a.vy = vy;
            a.vz = vz;

            // position calculation
            a.px += vx * time_step;
            a.py += vz * time_step;
            a.py += vz * time_step;

            /* ---
            REBOUND EFFECT
            --- */

            if(a.px <= 0){
                a.px = 0;
                a.vx = - a.vx;
            } else if(a.px >= size_enclosure){
                a.px = size_enclosure;
                a.vx = - a.vx;
            }

            if(a.py <= 0){
                a.py = 0;
                a.vy = - a.vy;
            } else if(a.py >= size_enclosure){
                a.py = size_enclosure;
                a.vy = - a.vy;
            }

            if(a.pz <= 0){
                a.pz = 0;
                a.vz = - a.vz;
            } else if(a.pz >= size_enclosure){
                a.pz = size_enclosure;
                a.vz = - a.vz;
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
        // write to file
        OutFile << universe[i].px << " " << universe[i].py << " " << universe[i].pz 
        << " " << universe[i].vx << " " << universe[i].vy << " " << universe[i].vz 
        << " " << universe[i].m;
        OutFile << endl;
    }

    OutFile.close();
}