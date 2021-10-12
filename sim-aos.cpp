//
// Created by ivan on 9/10/21.
//

#include <iostream>

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

    Object * universe[num_objects];

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
    
    auto * universe = parametersGeneration(argc, argcv);

    for(int iteration; iteration < num_iterations; iteration++){
        if(num_objects == 0){
            break;
        }
        for(int i = 0; i < num_objects; i++){
            for(int j = ++i; j < num_objects; j++){
                Object * a = universe[i];
                Object * b = universe[i];
                forcesComputation(a, b);
            }   
        }
    }
}