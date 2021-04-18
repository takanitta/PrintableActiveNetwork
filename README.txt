These codes are for a computer simulation of printable active networks described in a paper entitled ÅgA printable active network actuator built from an engineered biomolecular motorÅh published in Nature Materials (DOI: 10.1038/s41563-021-00969-6). The simulation method is described in the paper. The codes are written in Fortran 90 with openACC. These codes can simulate contraction of printable active networks in various conditions described in the paper by turning on and commenting out lines in NetworkNodeSimulation_Main.f90 and by changing parameter values in NetworkNodeSimulation_Parameters.f90 as described below.


Computation environment:
Simulations were run on PCs with Ubuntu OS installed with NVIDIA TITAN V or GeForce GTX 1080. The codes can be also run without GPU by commenting out all the lines with Åg! GPUÅh and turning on all the lines with Åg! CPUÅh in NetworkNodeSimulation_Main.f90.


An addition code needed to run the simulation:
We used Mersenne Twister to generate random numbers with CPU. Please download a Mersenne Twister from the website (http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/VERSIONS/FORTRAN/fortran.html).


Brief descriptions of parameters are given in NetworkNodeSimulation_Parameters.f90.


To simulate contraction of a circular active network:
1. Set ÅgNumNodeÅh in NetworkNodeSimulation_Parameters.f90 to set the number of nodes consisting of a circular active network.
2. Set ÅgRadiusÅh in NetworkNodeSimulation_Parameters.f90 to set the initial radius of a circular active network.
3. At ÅgSet Initial Distribution of NodesÅh in NetworkNodeSimulation_Main.f90, choose ÅgCALL InitialDistributionCircular(X,Y)Åh to generate s an initial distribution of nodes, and comment out others.
4. At ÅgSet contraction modeÅh in NetworkNodeSimulation_Main.f90, choose the line with ÅgNo contraction reversalÅh, and comment out the others.
5. Make sure that the nodes are not assigned as anchored (see below for more details).

To simulate contraction of a rectangular random active network fixed both ends:
1. Set ÅgNumNodeÅh in NetworkNodeSimulation_Parameters.f90 to set the number of nodes consisting of a rectangular active network.
2. Set ÅgX_SizeÅh and ÅgY_SizeÅh in NetworkNodeSimulation_Parameters.f90 to set the initial size of the rectangular active network.
3. At ÅgSet Initial Distribution of NodesÅh in NetworkNodeSimulation_Main.f90, choose ÅgCALL InitialDistributionRectangular(X,Y)Åh to generate an initial distribution of nodes, and comment out others.
4. At ÅgSet contraction modeÅh in NetworkNodeSimulation_Main.f90, choose the line with ÅgNo contraction reversalÅh, and comment out the others.
5. Turn on all the lines with Åg! Fixed AnchorÅh at the end of the lines, and comment out all the lines with Åg! Flexible AnchorÅh at the ends of the lines in NetworkNodeSimulation_Main.f90.
6. Set movable and anchored nodes. Nodes satisfying the condition in the if-statement in SUBROUTINE AssignNodeTypesSort in NetworkNodeSimulation_Func.f90 are assigned as movable, otherwise nodes are assigned as anchored.

To simulate contraction of a rectangular grid active network fixed both ends:
1. Set ÅgNumNodeÅh in NetworkNodeSimulation_Parameters.f90 to set the number of nodes consisting of a rectangular active network.
2. Set ÅgX_SizeÅh and ÅgY_SizeÅh in NetworkNodeSimulation_Parameters.f90 to set the initial size of a rectangular active network.
3. At ÅgSet Initial Distribution of NodesÅh in NetworkNodeSimulation_Main.f90, choose ÅgCALL InitialDistributionRectangularGrid(X,Y)Åh to generate an initial distribution of nodes, and comment out others.
4. At ÅgSet contraction modeÅh in NetworkNodeSimulation_Main.f90, choose the line with ÅgNo contraction reversalÅh, and comment out the others.
5. Turn on all the lines with Åg! Fixed AnchorÅh at the end of the lines, and comment out all the lines with Åg! Flexible AnchorÅh at the ends of the lines in NetworkNodeSimulation_Main.f90.
6. Set movable and anchored nodes. Nodes satisfying the condition in the if-statement in SUBROUTINE AssignNodeTypesSort in NetworkNodeSimulation_Func.f90 are assigned as movable, otherwise nodes are assigned as anchored.

To simulate repeated cycles of contraction/expansion of a rectangular random active network with a flexible lever:
1. Set ÅgNumNodeÅh in NetworkNodeSimulation_Parameters.f90 to set the number of nodes consisting of a rectangular active network.
2. Set ÅgX_SizeÅh and ÅgY_SizeÅh in NetworkNodeSimulation_Parameters.f90 to set the initial size of a rectangular active network.
3. Set ÅgReversalPeriodÅh in NetworkNodeSimulation_Parameters.f90 to set the period each contraction or expansion proceeds.
4. Set ÅgkLeverÅh in NetworkNodeSimulation_Parameters.f90 to set the spring constant of the flexible lever.
5. At ÅgSet Initial Distribution of NodesÅh in NetworkNodeSimulation_Main.f90, choose ÅgCALL InitialDistributionRectangular(X,Y)Åh to an initial distribution of nodes, and comment out others.
6. At ÅgSet contraction modeÅh in NetworkNodeSimulation_Main.f90, choose ÅgPeriodic contraction reversalÅh, and comment out others.
7. Turn on all the lines with Åg! Flexible AnchorÅh at the end of the lines, and comment out all the lines with Åg! Fixed AnchorÅh at the end of the lines in NetworkNodeSimulation_Main.f90.
8. Set movable and anchored nodes. Nodes satisfying the condition in the if-statement in SUBROUTINE AssignNodeTypesSort in NetworkNodeSimulation_Func.f90 are assigned as movable, otherwise nodes are assigned as anchored.

To simulate contraction of a rectangular random active network fixed both ends after repeated cycles of contraction/expansion:
1. Set ÅgNumNodeÅh in NetworkNodeSimulation_Parameters.f90 to set the number of nodes consisting of a rectangular active network.
2. Set ÅgX_SizeÅh and ÅgY_SizeÅh in NetworkNodeSimulation_Parameters.f90 to set the initial size of a rectangular active network.
3. Set ÅgReversalPeriodÅh in NetworkNodeSimulation_Parameters.f90 to set the period each contraction or expansion proceeds.
4. Set ÅgTSContractionONÅh in NetworkNodeSimulation_Parameters.f90 to set the onset of contraction.
5. At ÅgSet Initial Distribution of NodesÅh in NetworkNodeSimulation_Main.f90, choose ÅgCALL InitialDistributionRectangular(X,Y)Åh to generate an initial distribution of nodes, and comment out others.
6. At ÅgSet contraction modeÅh in NetworkNodeSimulation_Main.f90, choose ÅgPeriodic contraction reversal and then contractionÅh, and comment out others.
7. Turn on all the lines with Åg! Fixed AnchorÅh at the end of the lines, and comment out all the lines with Åg! Flexible AnchorÅh at the ends of the lines in NetworkNodeSimulation_Main.f90.
8. Set movable and anchored nodes. Nodes satisfying the condition in the if-statement in SUBROUTINE AssignNodeTypesSort in NetworkNodeSimulation_Func.f90 are assigned as movable, otherwise nodes are assigned as anchored.


Contact:
If you have any question, please contact me (Takahiro Nitta) by nittat@gifu-u.ac.jp.
