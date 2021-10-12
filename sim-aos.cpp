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
        //std::vector<double> f(3, 0);
};

// constants

const double g = 6.674 * pow(10, -11);

// Global variables
int num_objects;
int num_iterations;
int random_seed;
int size_enclosure;
int time_step;
double forces[1][1];  // TODO: quitar este apaño necesario para acceder en updatePosition()

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


void updatePosition(Object * obj, int num){
    double a[3] = {0};
    // acceleration calculation
    for(int k = 0; k <3; k++){
        for(int l = 0; l < num_objects; l++){
            a[k] += forces[3*num + k][3*l + k];
            //a[k] += obj->forces[3*l + k];
        }
        a[k] /= obj->m;
    }

    // velocity calculation
    double vx = obj->vx + a[0] * time_step;
    double vy = obj->vy + a[1] * time_step;
    double vz = obj->vz + a[2] * time_step;

    obj->vx = vx;
    obj->vy = vy;
    obj->vz = vz;

    // position calculation
    obj->px += vx * time_step;
    obj->py += vz * time_step;
    obj->py += vz * time_step;

    // rebound effect

    // p <= 0
    if(obj->px <= 0){
        obj->px = 0;
        obj->vx = - obj->vx;
    }

    if(obj->py <= 0){
        obj->py = 0;
        obj->vy = - obj->vy;
    }

    if(obj->pz <= 0){
        obj->pz = 0;
        obj->vz = - obj->vz;
    }

    // p >= size_enclosure
    if(obj->px >= size_enclosure){
        obj->px = size_enclosure;
        obj->vx = - obj->vx;
    }

    if(obj->py >= size_enclosure){
        obj->py = size_enclosure;
        obj->vy = - obj->vy;
    }

    if(obj->pz >= size_enclosure){
        obj->pz = size_enclosure;
        obj->vz = - obj->vz;
    }

    return;
}


double * forceComputation(Object * a, Object * b){
    if((a->px == b->px) && (a->py == b->py) && (a->pz == b->pz)){
        // TODO: cambiar esto, es posible que con el timestep no coincidan exactamente
        objectCollision(a, b);
        return 0;  // TODO: cambiar esta mierda a vectores y eliminar el vector en forces
    }
    double dx = b->px - a->px;
    double dy = b->py - a->py;
    double dz = b->pz - a->pz;
    
    double f[3];
    f[0] = (g * a->m * b->m * dx) / abs(pow(dx, 3));
    f[1] = (g * a->m * b->m * dy) / abs(pow(dy, 3));
    f[2] = (g * a->m * b->m * dz) / abs(pow(dz, 3));

    return f;
}


int main(int argc, const char * argcv[]){

    if (argc != 5 || argcv[0] < 0 || argcv[1] < 0 || argcv[2] < 0 || argcv[3] < 0 || argcv[4] < 0){
        cout << "sim-aos invoked with " << argc << "parameters." << endl << "Arguments: "<< endl << " num_objects: " << argcv[0] << endl << " num_iterations: " << argcv[1] << endl << " random_seed: " << argcv[2] << endl << " size_enclosure: " << argcv[3] << endl << " time_step: " << argcv[4] << endl ;
    }
    
    auto * universe = parametersGeneration(argc, argcv);

    double forces[3 * num_objects][3 * num_objects];  // TODO: update to use 2d vectors

    for(int iteration; iteration < num_iterations; iteration++){
        if(num_objects == 0){
            break;
        }
        for(int i = 0; i < num_objects; i++){
            for(int j = ++i; j < num_objects; j++){
                Object * a = universe[i];
                Object * b = universe[j];
                double * result = forceComputation(a, b);
                forces[3*i][3*j] = result[0];
                forces[3*j][3*i] = - result[0];
                forces[3*i + 1][3*j + 1] = result[1];
                forces[3*j + 1][3*i + 1] = - result[1];
                forces[3*i + 2][3*j + 2] = result[2];
                forces[3*j + 2][3*i + 2] = - result[2];
            } 
            updatePosition(universe[i], i);
        }
    }
}