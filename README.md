# DP-DICE-Bino

This repository contains the codes of manuscript "Flexible Privacy-Preserving Cardinality Query for Massive Distributed Datasets". 

The manuscript aims to solve the private distributed cardinality estimation (PDCE) problem and the private distributed intersection cardinality estimation (PDICE) problem, which are formalized in our manuscript. We propose a novel protocol ***DP-DICE-Bino***, which is computationally efficient and differentially private compared to previous state-of-the-art protocol MPC-FM. Experiments show that our DP-DICE-Bino achieves orders of magnitude speedup and reduces the estimation error by several times in comparison with state-of-the-arts under the same security requirements. 

The repository consists four parts: the codes of DP-DICE-Bino, DP-DICE, MPC-FM and the FMS sketch, which are all detailed in our manuscript. The implementation are developed mainly in C++. This repository implements our DP-DICE-Bino (and DP-DICE / MPC-FM) protocol based on [SPDZ-2](https://github.com/bristolcrypto/SPDZ-2). We evaluate the communication cost and time consumption based on this implementation and the results are displayed in Section 5.2 in our manuscript. 

## DP-DICE-Bino / DP-DICE / MPC-FM

### Requirements:
- OS: Ubuntu 16.04 LTS
- GCC  7
- MPIR library  3.0.0
- libsodium library  1.0.11
- NTL library  11.0.0
- valgrind, version 3.14.0
- farmhash 1.1
- openssl 1.1.1

### To compile DP-DICE-Bino (or DP-DICE / MPC-FM)
Download this repository and run the following commands in terminal: 
```bash
DP-DICE-Bino:
cd DP_DICE_Bino/
make clean
make all

DP-DICE:
cd DP_DICE/
make clean
make all

MPC-FM:
cd MPC_FM/
make clean
make all
```
For each of the folders (DP_DICE_Bino, DP_DICE or MPC_FM), after compilation, we get executable files `DP.x`, `CP.x`,  `pairwise-offline.x` and `Server.x`, where `pairwise-offline.x` will be used for random number generation in the offline preparation phase; `Server.x`,  ` DP.x` and `CP.x` will be used for protocol execution in the online phase. 

The complete process of the DP-DICE-Bino (or DP-DICE / MPC-FM) protocol includes the offline preparation phase and the online phase. 

### Offline preparation
In the offline preparation phase, generating the secret random numbers (**Triple, Rand, Rand2, and RandExp**) used in the online step is necessary. The following is the description of the parameters for random number generation. After compiling DP-DICE-Bino (or DP-DICE / MPC-FM), we can get the executable file `pairwise-offline.x` and use it to generate random numbers. `pairwise-offline.x` can be executed with the following parameters:

```shell
./pairwise-offline.x  -o {Write generated results to files} -h {host of CP0 (default: localhost)} -x  {number of threads }  -f {the modulus used for the SPDZ framework} -N {number of CPs} -p {Current CP sequence number(From 0 to N-1)} -n  {number of Triples to be generated} -nr {number of  Rands to be generated} -nr2  {number of Rand2s to be generated} -nrx{number of  RandExp s to be generated} -E {the bit length of the input plaintext for SPDZ}
```

For a more detailed parameter description of `pairwise-offline.x`, you can execute `./pairwise-offline.x` directly without any parameter to obtain.

The generated random numbers in shared form and MAC-Key are stored in DP_DICE_Bino/Player-Data (or DP_DICE/Player-Data, or MPC_FM/Player-Data)

#### Example for Offline preparation phase

When there are 2 CPs to generate 2000000 **Triples**.
```bash
# On CP 0
./pairwise-offline.x -o -f 60 -N 2 -p 0 -n 2000000
# On CP 1
./pairwise-offline.x -o -f 60 -N 2 -p 1 -n 2000000
```
When there are 2 CPs to generate 2000000 **Rand**.
```bash
# On CP 0
./pairwise-offline.x -o -f 60 -N 2 -p 0 -nr 2000000 -r
# On CP 1
./pairwise-offline.x -o -f 60 -N 2 -p 1 -nr 2000000 -r
```
When there are 2 CPs to generate 2000000 **Rand2**.
```bash
# On CP 0
./pairwise-offline.x -o -f 60 -N 2 -p 0 -nr2 20000 -r2
# On CP 1
./pairwise-offline.x -o -f 60 -N 2 -p 1 -nr2 20000 -r2
```
When there are 2 CPs to generate 2000000 **RandExp30**.
```bash
# On CP 0
./pairwise-offline.x -o -f 60 -N 2 -p 0 -nrx 20000 -rx -E 30
# On CP 1
./pairwise-offline.x -o -f 60 -N 2 -p 1 -nrx 20000 -rx -E 30
```
When there are 2 CPs to generate 2000000 **RandExp18**.
```bash
# On CP 0
./pairwise-offline.x -o -f 60 -N 2 -p 0 -nrx 20000 -rx -E 18
# On CP 1
./pairwise-offline.x -o -f 60 -N 2 -p 1 -nrx 20000 -rx -E 18
```
### Online 
#### Example for Online phase
To facilitate testing the execution of the DP-DICE-Bino (or DP-DICE / MPC-FM) protocol in the online phase, we give a simple test example. Note that the necessary offline data needs to be generated before executing the online phase protocol.
In this example all the DHs together hold about $ 10^6 $ unique data items, with FMS sketch parameter  $m=1024$.

##### To establish the communication network
The command is running on nameserver terminal to establish the communication network of each computing party
```bash
./Server.x 2 5000 1
```
#####  Running the CPs
Open two new commands and run the CPs
```bash
./CP.x  -lgp 60 -np 2 -p 0 -sec 40 -ndp 1 -oM 1000 -uN 1000000 -lan 1 -nt 1 
./CP.x  -lgp 60 -np 2 -p 1 -sec 40 -ndp 1 -oM 1000 -uN 1000000 -lan 1 -nt 1
```
#####  Running the DHs
Run the DH in the Server terminal
```bash
./DP.x -ndp 1 -ncp 2 -p 1 -lgp 60 -oM 1000 -uN 1000000 -nt 1
```
##### parameters for nameserver
```bash
./Server.x {number of computation parties} 5000 {number of threads}
```

##### parameters for CP
```bash
./CP.x  -lgp {length of the prime} -np {number of CPs} -p {party number} -sec {security parameter} -ndp {number of DHs} -oM {number of FMS sketch} -uN { number of unique items} -lan {whether in lan environment} -nt {number of threads} -h {ip of DH}
```

##### parameters for DH
```bash
./DP.x -ndp {number of DHs} -ncp {number of CPs} -p {party number} -lgp {length of the prime} -oM {number of FMS sketch} -uN { number of unique items} -ip {ip files}
```



## FMS sketch

### Requirements:

- OS: Ubuntu 16.04 LTS
- GCC version >= 7.5.0
- cmake version >= 3.25.2

### To compile FMS sketch

Run the following commands in terminal: 

```cpp
cd FMS
mkdir build
cd build
cmake ..
make
./FMS 
```

Afterwards, we will obtain the estimation errors of our FMS sketch and the classic LL sketch and FM sketch. The results are saved in the build folder.

```
fms_m.txt:   the RMSD under different M without noisy
fms_noisy_m.txt:   the RMSD under different M with noisy
fms_noisy_epsilon.txt:   the RMSD under different epsilon with noisy
```

Three mechanisms are used to preserve differential privacy. They are identified by numbers in the file:

```
(1) Gaussian mechanism 

(2) Distributed Discrete Gaussian mechanism 

(3) MPCDP-CP 
```

The cardinality size of two sets is both  $10^5$ and the intersection cardinality size is $5000$ by default. If you want to change the cardinality size, you can run the following commands in terminal after the make operation as an alternative.

```cpp
// n1: car of set S1;  n2: car of set S2;  n3: car of set S1 union S2;  n4: car of set S1 inter S2
./FMS n1 n2 n4
// such as the default value
./FMS 100'000 100'000 5000
```

### To compile other sketch (LL sketch, FM sketch)

Our project also provides the implementation of FM sketch and LL sketch, if you want to get the error of the other two sketches as a comparison, you can run the following commands in terminal:

```cpp
// Output the error of three types of sketch
./FMS ll fm
// you can change the cardinality size while getting the error of three types of sketch 
./FMS ll fm c 100'000 100'000 5000
./FMS c 100'000 100'000 5000 ll fm
```

### File structure

- build：Create a build directory to store the compilation products, which can avoid mixing the compilation products with the code files, generate a makefile file in the build directory, and then call the compiler to actually compile and link the project. (If you build in a different way, you can delete the folder.) 

- noise：Noise generation file. in addition to the generation of Laplace noise, includes various other types funtions to generate noise such as the Gaussian mechanism, the Distributed Discrete Gaussian mechanism, and the MPCDP-CP.
- MurmurHash3：MurmurHash3 library functions.
- FMsketch: Implementation of FM sketch(Flajolet-Martin sketch).
- LLsketch: Implementation of LogLog sketch.
- FMSsketch: Our FMS sketch method, including multi-threaded processing part of the code.
- batch_est:  Batch processing functions for cardinality estimation.
- core.cpp: The entrance of the project.







