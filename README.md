# CA-Project1
By team 89-02: Luis Daniel Casais Mezquida, Ivan Darío Cersósimo & Pablo Ruiz Fernández

## Project description
This project has as a fundamental goal to help students to acquire concerns on performance of sequential programs and to familiarize with their optimization. In addition, they will get familiarize with performance evaluation techniques. <br/>
More specifically, the project focuses in the development of sequential software in the C++ programming language (including the latest versions). <br/>
<br/>
To carry out this project you must implement in C++ a gravitation simulation application for a set of objects. The program must only use sequential programming techniques without making use of any kind of concurrency or parallelism.<br/>
Two versions of the program will be implemented using the arrays of structures (aos) and structures of arrays (soa) techniques.
For a successful execution of the project, the following sections of this statement provides the requirements that must be satisfied by the developed programs. Besides, a binary fle with a version of the executable program to make easier results comparison.<br/>

### The gravitational interaction problem
The problem you must solve consists of the simulation over time of the positions of a set of objects in an enclosure. This computation is performed during the program execution at regular time intervals with length ∆t. <br/>
The program must perform the simulation in a 3D space. We consider the space a closed enclosure. Consequently, when an object impacts a plane in the border of the enclosure, a rebound happens, changing the direction of the object motion. <br/>
Finally, in every time increment ∆t possible collisions between different objects are checked. <br/>
In this way, it is expected that students simulate the behavior of N different objects moving in a 3D space. <br/>

### Program execution and input parameters
Teams will develop two programs named ``sim-soa`` (structure of arrays version) and ``sim-aos`` (array of structures version). <br/>
The program will take the following parameters for its correct execution:
* ``num_objects``: integer number, greater or equal than 0, indicating the number of objects that will be simulated.
* ``num_iterations``: integer number, greater or equal than 0, indicating the number of iterations (time steps) that will be simulated.
* ``seed``: positive integer number to be used as a seed for random distribution generator functions. Note: Two simulations with the same parameters but different seed will lead to different scenarios
* ``size_enclosure``: real positive number indicating the size of the simulation enclosure. The enclosure is considered to be a perfect cube with a vertex in the coordinates origin and with side equal to size_enclosure.
* ``time_step``: real positive number indicating the time increment for each iteration of the simulation.

When the number of parameters is wrong, the program will terminate with code error ``-1`` and an error message to the standard error output. <br/>
```
user@machine:~$ ./sim-soa 100 2000 1
sim-soa invoked with 3 parameters.
Arguments:
num_objects: 100
num_iterations: 2000
random_seed: 1
size_enclosure: ?
time_step: ?
```
When any parameter is wrong, the program will terminate with error code ``-2`` and an error message to the standard error output.
```
user@machine:~$ ./sim-soa 100 2000 1 -1000000 0.1
Error: Invalid box size
sim-soa invoked with 5 parameters.
Arguments:
num_objects: 100
num_iterations: 2000
random_seed: 1
size_enclosure: -1000000
time_step: 0.1
```
Once all input parameters have been processed and initial positions and velocities have been processed for all simulation elements, the program will write this information to a file named ``init_conf.txt``, with the following format:
1. A header line including:
    1. Size of the enclosure.
    2. Time step.
    3.  Number of initial objects.
2. One line per object where the position (x, y, z coordinates), velocity (x, y, z components) and mass must be written. This seven values must be separated by blanks. All real values will be printed with tree decimal places in fractional part.

### Simulation parameters generation
Object position in the 3D space (with coordinates in floating point double precision), must be generated following a random distribution with the seed that has been provided as program argument.
* A single pseudo-random number generator engine must be used for all distributions. It must be a 64 bits Mersenne-Twister generator (``std::mt19937_64``). This engine will be initiated with a random seed received as a parameter by the program. Information is available at https://en.cppreference.com/w/cpp/numeric/random/mersenne_twister_engine.
* All position values will be generated using a uniform distribution with lower bound equal to 0 and upper bound equal to size_enclosure, using an object of type ``uniform_real_distribution``. Information is available at https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution.
* All mass values will be generated using a normal distribution with mean equal to 1e21 and standard deviation equal to 1e15.

### Attraction forces computation
__Gravitational force computation:__ Gravitational constant (_G_), object masses of both objects, _mi_ & _mj_, and position vectors, _pi_ & _pj_, will be used.<br/>
![img1](https://user-images.githubusercontent.com/78721925/136654279-7a3a0f32-0df7-4935-a5a7-87213725f0e8.png) <br/>
The universal gravitational constant will be the value 6.674e−11. <br/>
For an object _i_, this force is applied with a positive value and for object _j_ with negative value. <br/>
An object does not apply any force on itself. <br/>
__Acceleration vector computation:__ Acceleration (_veca_) is determined using its mass (_m_) and all the force contributions. <br/>
![img2](https://user-images.githubusercontent.com/78721925/136654334-aa1fbef8-26fd-4182-9eb3-519a5dfab67e.png) <br/>
__Velocity vector computation:__ Velocity (_v_) is determined using its acceleration (_a_) and the time step (∆t). <br/>
![img3](https://user-images.githubusercontent.com/78721925/136654366-5c92c836-0466-4d4a-9040-aac3a5c1307a.png) <br/>
__Position vector computation:__ Position is determined using the velocity and the time step. <br/>
![img4](https://user-images.githubusercontent.com/78721925/136654390-645e5d97-700a-4483-b3d6-8ec7ffbc5f3f.png) <br/>
__Rebound effect:__ If the position passes one of the enclosure bounds:
1. The object is placed in the corresponding bound plane:
```
px ≤ 0 ⇒ px ← 0
py ≤ 0 ⇒ py ← 0
pz ≤ 0 ⇒ pz ← 0
px ≥ tam ⇒ px ← tam
py ≥ tam ⇒ py ← tam
pz ≥ tam ⇒ pz ← tam
```
2. Besides the corresponding velocity component changes its sign.

### Object collisions
When two objects collide they collapse into an object with larger mass. Two objects are considered to collide when their distance is less than the unit. <br/>
In this case both objects (_a_ and _b_) lead to a single object (_c_) with the following properties:
* New mass results from mass addition (_mc_ = _ma_ + _mb_).
* New velocity results from velocity addition (_vc_ = _va_ + _vb_).

### Data storage
Once all simulation iterations have finished the program shall store final data into a file named ``final_config.txt``, with the same format and contents than the file ``init_config.txt``.

##  Sequential version development
This task consists in the development of the sequential version of the described application in C++17 (or C++20). <br/>
All your source files must compile without problems and shall not emit any compiler warning. In particular, the following compiler warning specific flags must be enabled:
```
-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors
```
Keep also in mind that you will have to perform all evaluations with compiler optimizations enabled (compiler options ``-O3`` and ``-DNDEBUG``). You can easily get this with ``-DCMAKE_BUILD_TYPE=Release``. <br/>
You are allowed to use additional compiler flags as long as you document them in the project report and justify its use. Such flags must be properly included in the CMake configuration file.

##  Sequential performance evaluation
This task consists of the performance evaluation of the sequential application. <br/>
To carry out the performance evaluation you must measure the application execution time. You are expected to represent graphically your results. Keep in mind the following considerations:
* Run each experiment a given number of times and take the average value. You are recommended to perform each experiment at least 10 or more executions so that you can give a confidence interval.
* Study results for different object populations. Consider cases with 1000, 2000 and 4000
* Study results for different number of iterations: Consider cases with 50, 100 and 200.

Represent graphically total execution times. Represent graphically average time per iteration. <br/>
Include in your report the conclusions you may infer from results. Do not limit to describing data. You must search for a convincing explanation of results.
