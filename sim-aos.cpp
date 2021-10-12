//
// Created by ivan on 9/10/21.
//

#include <iostream>
#include <vector>
#include <math.h>

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

    Object * universe[num_objects];

    return universe;
}


void objectCollision(Object * a, Object * b){
    // TODO: Ivan

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
    
    auto * universe = parametersGeneration(argc, argcv);

    double forces[3 * num_objects][3 * num_objects];  // TODO: update to use 2d vectors

    for(int iteration; iteration < num_iterations; iteration++){
        if(num_objects == 0){
            break;
        }
        for(int i = 0; i < num_objects; i++){
            for(int j = ++i; j < num_objects; j++){
                Object * a = universe[i];
                Object * b = universe[i];
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