# Greedy search
# for the best Minimally Complex spin Models (MCM)

This program allows to **find community structures in binary data**, taking into account possible **high order patterns** of data in the detection of the communities (i.e., possible high-order correlations between the variables). It is complementary to the program available [here](https://github.com/clelidm/MinCompSpin), which performs an exhaustive search for the best community. This program performs a greedy search for the best community, which is more appropriate for large systems (more than `15` spin variables).


The idea of the algorithm is based on performing statistical inference with a family of spin models (maximum entropy models for binary data) that have minimal information theoretic complexity. These models are called Minimally Complex Models (MCM). Details can be found in Ref. [1].

The idea of the algorithm is based on performing statistical inference with a family of spin models (maximum entropy models for binary data) with minimal information theoretic complexity. Details can be found in Ref. [1].

----

This repository contains a code developed for the paper Ref. [1] on *Statistical Inference of Minimally Complex Models* available in [arXiv:2008.00520](https://arxiv.org/abs/2008.00520). The code performs a greedy search for the best Minimally Complex Spin Model (MCM) on a basis provided by the user. This greedy approach is useful for systems with a large number of variables, `n > 15`. 

The code performs an hierarchical merging procedure to find an optimal MCM on the basis provided by the user.

To do so, the code start by re-writing the data in the new basis provided. It then start `r` communities

The code performs an hierarchical merging procedure to find an optimal MCM. It starts with `r` communities, where `r` is the number of basis operators provided.

is the each of them having 


We start from the IM based on the basis operators b∗ identified above, which is an MCM with n ICC of rank r=1. We merge two ICCs Ma and Ma′ in all possible ways. Among these combinations, we identify the pair that yields a maximal increase of the evidence in Eq. (8) and merge the corresponding ICCs. This procedure generates an approximation of the MCM that achieves a maximal value of the evidence along the hierarchical merging process, as the number of ICCs varies from n to 1.

The code go through all possible MCMs of a given rank `r`, where an MCM is defined by a given partition of the `r` basis operators provided (see paper). The comparison between models is based on their evidence (posterior probability that the model generates the data, integrated over its parameter values). The selected model is the one with the largest evidence.

One big advantage of this family of models (the MCMs) is that the computation of the evidence doesn’t require fitting the parameters of the model, which allows to significantly accelerate the comparison between models. The limiting factor of this code is the exhaustive search for the best model, as the space of models is exponentially increasing with `r`.

For a given number `r` of basis elements, the number of possible partitions of this set of `r` elements is equal to the Bell number of `r`, which grows roughly exponentially with `r`. Besides, the running time of the program also grows linearly with the number of different states observed in the dataset (so approximatively linearly with the number of datapoints in the dataset). For these reasons, this code is good for use on small systems, typically with `r <~15` variables. Note that for `r~15`, the code will take several days to perform the exhaustive search.

To efficiently generate all possible set partitions of these `r` operators, we used the algorithm described in Ref. [2] and in Ref. [3] (Algorithm E) that generates the partitions in Gray-code order. Thus, the code goes through all possible partitions by changing only one bit for each new partition. 

[1]  C. de Mulatier, P. P. Mazza, M. Marsili, *Statistical Inference of Minimally Complex Models*, [arXiv:2008.00520](https://arxiv.org/abs/2008.00520)

[2]  Ehrlich, Gideon. *Loopless algorithms for generating permutations, combinations, and other combinatorial configurations.* Journal of the ACM (JACM) 20.3 (1973): 500-513.

[3]  D.E. Knuth, *The Art of Computer Programming*, Volume 4, Combinatorial Algorithms: Part 1.693 (Addison-Wesley Professional), (2011).

## Requirements

The code uses the C++11 version of C++.

**To compile:** `g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp Basis_Choice.cpp RandomPart.cpp info_quant.cpp Generate_data_exact.cpp`

**To execute:** `./a.out`

## Examples

All the useful functions that can be called from `int main()` are declared at the beginning of the `main.cpp` file, and described in the sections below. For hands-on and simple tests of the program, check the examples in the function `int main()` of the `main.cpp` file. In the input folder, we provided the binary dataset `Dataset_Shapes_n9_N1e5.dat`, which is the dataset used in the section "Boolean description of a binary dataset" of Ref. [1]. 

## License

This code is an open source project under the GNU GPLv3.
