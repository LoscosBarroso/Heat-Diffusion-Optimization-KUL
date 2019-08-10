#Project Mathematical Engineering: Heat Diffusion

Year: 2018-2019
Author 1: Daniel, Loscos Barroso
Author 2: Onurcan, Aktürk
Author 3: Eva Karagiannakou

Project Description:
* C++ implementation of a PDE solver and optimizer for a given heat diffusion problem. The PDE was modelized using a finite differences scheme and solved using the transpose method. The cado also includes several matlab scripts that we used in the design process and to print the results.
Needed libraries:
* Eigen (included in zip), NlOpt (installed in the computer labs),  
Instructions:
* All main programs generate a tempsX.m which can be run in MATLAB to visualize the results where X is the mesh size specified in main.cpp.
* To run the main program: Go to src/ command 'make main' then ./main
* To run the main program in debug version which prints out most of the results: Go to src/ command 'make maind' then ./maind
* To run the finite differences debug mode which prints both gradients computed by the adjoint method and finite differences: Go to src/ command 'make mainfd' then ./mainfd
* To see that the numerical scheme works, go to /PDEexamples, compile solvePDE.cpp by commanding 'make solvePDE' then run with the command './solvePDE name' which reads the metal concentration matrix in file name.txt then solves the system then generates a name.m file which can be run to see the results visually.
* You can use or modify our code for demonstrations concerning the course. Or for any other purpose contemplated in the license file.
