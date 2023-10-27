# Fms sketch

## Requirements:

- OS: Ubuntu 16.04 LTS
- GCC version >= 7.5.0
- cmake version >= 3.25.2

## To compile FMS sketch

Download this repository and run the following commands in cmd.

```cpp
mkdir build
cd build
cmake ..
make
./FMS 
```

Afterwards, we obtained three types of sketch errors(ll sketch, fm sketch, fms sketch) and saved the results in the build folder.

- fms_m.txt: the RMSD under different M without noisy
- fms_noisy_m.txt: the RMSD under different M with noisy
- fms_noisy_epsilon.txt: the RMSD under different epsilon with noisy

It should be noted that we have added three forms of noise separately,  which are identified by numbers in the file:

(gaussian_noise(1) CDP(2) Binomial mechanism(3)) 

At the same time, the cardinality size of two sets is both $10^5$ and the intersection cardinality size is $5000$ by default. If you want to change the cardinality size, you can run the following commands in cmd after the make operation as an alternative.

```cpp
// n1: car of set S1;  n2: car of set S2;  n3: car of set S1 union S2;  n4: car of set S1 inter S2
./FMS n1 n2 n4
// such as the default value
./FMS 100'000 100'000 5000
```

 Of course, you can also simply modify the parameters of the function batch_xx_est to change the cardinality size or other parameters,

## Code introduction

- build：Create a build directory to store the compilation products, which can avoid mixing the compilation products with the code files, generate a makefile file in the build directory, and then call the compiler to actually compile and link the project. (If you build in a different way, you can delete the folder.) 

- noise：Noise generation file. in addition to the generation of Laplace noise, includes various other types funtions to generate noise such as the Gaussian mechanism, the Distributed Discrete Gaussian mechanism, and the MPCDP-CP.

- MurmurHash3：MurmurHash3 library functions.

- FmsSketch: Our Fms sketch method, including multi-threaded processing part of the code.
- batch_est:  Batch processing functions for cardinality estimation.

- core.cpp: The code runs into the portal, including the main function.
