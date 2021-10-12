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
        Object(double x, double y, double z, int mass){
            px = x;
            py = y;
            pz = z;
            m = mass;
        }
        double px, py, pz, vx, vy, vz;
        int m;
};

// constants
const double g = 6.674 * pow(10, -11);

// Global variables
int num_objects;
int num_iterations;
int random_seed;
int size_enclosure;
int time_step;
Force ** forces;
//Force forces[1][1];  // TODO: quitar este apaño necesario para acceder en updatePosition()

auto parametersGeneration(const char * argcv[]){

    // parameters init & casting
    num_objects = (int) argcv[0];
    num_iterations = (int) argcv[1];
    random_seed = (int) argcv[2];
    size_enclosure = (int) argcv[3];
    time_step = (int) argcv[4];

    // distribution generation
    mt19937_64 gen64;  // generate object
    uniform_real_distribution<> dis(0.0, size_enclosure);
    normal_distribution<> d{pow(10, 21),pow(10, 15)};

    // introduce seed
    gen64.seed(random_seed);

    Object * universe = (Object*)malloc(sizeof(Object) * num_objects);

    ofstream MyFile("init_config.txt");  // open file
    
    MyFile << argcv[3] << argcv[4] << argcv[0];
    MyFile << endl;

    // populate
    for (int i = 0; i < num_objects; i++){
        universe[i] = Object(dis(gen64), dis(gen64), dis(gen64), d(gen64));
        // write to file
        MyFile << universe[i].px << universe[i].py << universe[i].pz << universe[i].vx << universe[i].vy << universe[i].vz << universe[i].m;
        MyFile << endl;
    }

    MyFile.close();
    return universe;
}


void objectCollision(Object a, Object b){
    // merge objects into a
    a.m = a.m + b.m;
    a.vx = a.vx + b.vx;
    a.vy = a.vy + b.vy;
    a.vz = a.vz + b.vz;

    delete &b;
    return;
}


void updatePosition(Object obj, int num){
    /*
    acceleration calculation
    */

    Force acc(0, 0, 0);
    for(int k = 0; k < num_objects; k++){
        Force f = forces[num][k];
        acc.x += f.x;
        acc.y += f.y;
        acc.z += f.z;
    }
    acc.x /= obj.m;
    acc.y /= obj.m;
    acc.z /= obj.m;

    /* 
    velocity calculation
    */
    double vx = obj.vx + acc.x * time_step;
    double vy = obj.vy + acc.y * time_step;
    double vz = obj.vz + acc.z * time_step;

    obj.vx = vx;
    obj.vy = vy;
    obj.vz = vz;

    /*
    position calculation
    */
    obj.px += vx * time_step;
    obj.py += vz * time_step;
    obj.py += vz * time_step;

    /*
    rebound effect
    */

    // p <= 0
    if(obj.px <= 0){
        obj.px = 0;
        obj.vx = - obj.vx;
    }

    if(obj.py <= 0){
        obj.py = 0;
        obj.vy = - obj.vy;
    }

    if(obj.pz <= 0){
        obj.pz = 0;
        obj.vz = - obj.vz;
    }

    // p >= size_enclosure
    if(obj.px >= size_enclosure){
        obj.px = size_enclosure;
        obj.vx = - obj.vx;
    }

    if(obj.py >= size_enclosure){
        obj.py = size_enclosure;
        obj.vy = - obj.vy;
    }

    if(obj.pz >= size_enclosure){
        obj.pz = size_enclosure;
        obj.vz = - obj.vz;
    }

    return;
}


Force forceComputation(Object a, Object b){
    // distance
    double dx = b.px - a.px;
    double dy = b.py - a.py;
    double dz = b.pz - a.pz;
    double distance = sqrt(dx*dx + dy*dy + dz*dz);

    if(distance < 1){
        objectCollision(a,b);
        return;
    }
    
    Force f(
        (g * a.m * b.m * dy) / abs(pow(dy, 3)),
        (g * a.m * b.m * dy) / abs(pow(dy, 3)),
        (g * a.m * b.m * dz) / abs(pow(dz, 3))
    );

    return f;
}


int main(int argc, const char ** argcv){

    // check arguments
    if (argc != 5 || num_objects < 0 || num_iterations < 0 || random_seed < 0 || size_enclosure < 0 || time_step < 0){
        cout << "sim-aos invoked with " << argc << "parameters." << endl << "Arguments: "<< endl << " num_objects: " << argcv[0] << endl << " num_iterations: " << argcv[1] << endl << " random_seed: " << argcv[2] << endl << " size_enclosure: " << argcv[3] << endl << " time_step: " << argcv[4] << endl ;
    }
    
    auto * universe = parametersGeneration(argcv);

    Force ** forces = new Force forces *[num_objects];
    for (int i = 0; i < num_objects; ++i) {
        forces[i]= new Force forces[num_objects];
    }
    //vector<vector<Force>> f(num_objects, num_objects);

    for(int iteration; iteration < num_iterations; iteration++){
        if(num_objects == 0){
            break;
        }
        for(int i = 0; i < num_objects; i++){
            for(int j = ++i; j < num_objects; j++){
                Object a = universe[i];
                Object b = universe[j];
                
                Force fa = forceComputation(a, b);
                Force fb(- fa.x, -fa.y, -fa.z);

                forces[i][j] = fa;
                forces[j][i] = fb;
            } 
            updatePosition(universe[i], i);
        }
    }
}