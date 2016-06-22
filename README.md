#Cardiac Modeling

This repository is jointly curated by members of the [RADICAL team](http://radical.rutgers.edu/) at Rutgers University and the [CBL team](https://people.rit.edu/lxwast/Linwei_Wangs_HomePage/CBL.html) at Rochester Institute of Technology.

The objective is to reduce the runtime of the cardiac model calculations. This will be done by refining to the underlying algorithms, identifying and increasing parallelism in existing code, and extending the single machine code to run on HPCs 

###TODOs###
* Identifying possible bottlenecks in code which can benefit from parallelism
* Remove dependencies on Windows-dependent libraries so that the code can run on Linux environments
* Understanding the structure of the existing code and how Ensemble Toolkit can increase concurrency
