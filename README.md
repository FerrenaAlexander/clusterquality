# clusterquality

This repo serves as a report for a class project performed in Spring 2019 aimed at hyperparameter optimization during scRNAseq clustering via the Louvain method. Namely, it is a k-fold cross-validation strategy involving perturbation of subsampled data followed by assesment of perturbation effects to final cluster label identification.

Benefits include that this is an attempt at mitigating the need for a priori "guessing" at optimal hyperparameters by providing a evidence for a scored best choice. Drawbacks include that this method is computational expensive and time-consuming to perform. Code optimization and parallelization in the future may mitigate these concerns.
