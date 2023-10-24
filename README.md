# DP-DICE-Bino

This repository contains the relevant codes of manuscript "Flexible Privacy-Preserving Cardinality Query for Massive Distributed Datasets". 

Counting the number of distinct elements (cardinality) distributed over multiple data holders (DHs) is fundamental with many applications ranging from crowd counting to network monitoring. Counting the cardinality of the intersection set between two groups of DHs also has practical significance. Although many space and computationally efficient methods (e.g., Flajolet-Martin (FM) sketch and HyperLogLog sketch) have been proposed to solve the problems, they are insecure when considering the privacy of each DH's personal dataset. Despite a recently proposed protocol that implements the FM sketch on a secret-sharing based multiparty computation (MPC) framework for solving the problem of private distributed cardinality estimation (PDCE), we observe that the protocol is not differentially private and is computationally expensive. To address the issues, we propose a novel protocol ***DP-DICE-Bino***, which is computationally efficient and differentially private for solving the PDCE problem. The DP-DICE-Bino is flexible in that it can compute the cardinality of any group of DHs without bothering them every time. Besides, our DP-DICE-Bino can also handle the private distributed intersection cardinality estimation (PDICE) problem. Experiments show that our DP-DICE-Bino achieves orders of magnitude speedup and reduces the estimation error by several times in comparison with state-of-the-arts under the same security requirements.

The repository consists four parts: the codes of DP-DICE-Bino, DP-DICE, MPC-FM and the FMS sketch. 

## DP-DICE-Bino





## DP-DICE





## MPC-FM





## FMS sketch









