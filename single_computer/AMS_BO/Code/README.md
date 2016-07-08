### Test Environment - Ming
* Windows 7
* MATLAB r2015a (32-int)
* .NET Framework 4.6

### Required Environment
* Windows 7+
* MATLAB
* .NET Framework 4
* Windows SDK 7.1

### Instructions to set up Running Environment
1. First load the entire folder of code into MATLAB, found in directory `<some_dir>/cardiac_modeling/single_computer/AMES_BO/Code`,
referred to as **base_dir** going forward
2. Check to see if MATLAB can compile C++ code. In order to compile the C++ files in order to run the script, check to see if MATLAB has a C++ compiler.
  1. Check using MATLAB terminal with command `mex -setup C++` (or some variant of C++ [`c++, cpp`])
  2. If there is no such compiler, download `Windows SDK 7.1`. This requires `.NET Framework **4**`. If you have a newer version,
    uninstall them and install the specified framework.
  3. Uninstall any `Microsoft Visual C++ 2010 Redistributale` packages currently installed
  4. Install `Windows SDK 7.1`
3. In MATLAB, navigate to directory `**base_dir**/GP_BO_optimization/GP_BO/optimizer/bobyqa`, 
and run the command `mex(strcat('-I"',pwd,'"'), 'bobyqa_alg.cpp')`
