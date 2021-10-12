//
// Created by ivan on 9/10/21.
//
using namespace std;

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <cassert>
#include <fstream>

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

// Global variables
int num_objects;
int num_iterations;
int random_seed;
int size_enclosure;
int time_step;

auto parametersGeneration(int argc, const char * argcv[]){
    // TODO: Pablo
    num_objects = (int) argcv[0];
    num_iterations = (int) argcv[1];
    random_seed = (int) argcv[2];
    size_enclosure = (int) argcv[3];
    time_step = (int) argcv[4];

    mt19937_64 gen64;
    uniform_real_distribution<> dis(0.0, size_enclosure);
    normal_distribution<> d{10^21,10^15};

    gen64.seed(random_seed);

    Object * universe = (Object*)malloc(sizeof(Object) * num_objects);

    ofstream MyFile("init_config.txt");
    
    MyFile << argcv[3] << argcv[4] << argcv[0];
    MyFile << endl;

    for (int i = 0; i < num_objects; i++){
        universe[i] = Object(dis(gen64), dis(gen64), dis(gen64), d(gen64));
        MyFile << universe[i].px << universe[i].py << universe[i].pz << universe[i].vx << universe[i].vy << universe[i].vz << universe[i].m;
        MyFile << endl;
    }

    MyFile.close();
    return universe;
}

 void objectCollision(Object * a, Object * b){
    // TODO: Ivan
    a->m = a->m +b->m;
    a->vx = a->vx + b->vx;
    a->vy = a->vy + b->vy;
    a->vz = a->vz + b->vz;
    delete b;
    return;
}

void forcesComputation(Object * a, Object * b){
    // TODO: Luisda
    int force = 0;


    if(force == 0){
        objectCollision(a, b);
    }

    return;
}

int main(int argc, const char * argcv[]){

    if (argc != 5 || argcv[0] < 0 || argcv[1] < 0 || argcv[2] < 0 || argcv[3] < 0 || argcv[4] < 0){
        cout << "sim-aos invoked with " << argc << "parameters." << endl << "Arguments: "<< endl << " num_objects: " << argcv[0] << endl << " num_iterations: " << argcv[1] << endl << " random_seed: " << argcv[2] << endl << " size_enclosure: " << argcv[3] << endl << " time_step: " << argcv[4] << endl ;
    }
    
    auto * universe = parametersGeneration(argc, argcv);

    for(int iteration; iteration < num_iterations; iteration++){
        if(num_objects == 0){
            break;
        }
        for(int i = 0; i < num_objects; i++){
            for(int j = ++i; j < num_objects; j++){
                Object * a = universe[i];
                Object * b = universe[j];
                forcesComputation(a, b);
            }   
        }
    }
}