# porousRedoxFoam
heterogeneous reverisble reactions of non-stoichiometric solid metal oxide

This folder contains files and programs created under 
GNU GPL v3 License
by Mario Goddy Zuber 2024

(c) 2024 ETH Zurich, Mario Goddy Zuber

if you use any part of this work please cite the scientific contribution:
	
	M. Zuber, Dry Redox Reforming with Concentrated Solar Energy. 
	Thesis. ETH Zurich. Doctor of Sciences, 2024.

	and any subsequent paper by M.Zuber on the use of the solver


in addition, as the solver is largely based on the work by P. Zuk,
please cite the scientific contribution:

	OpenFOAM solver for thermal and chemical conversion in porous media
	Pawel Jan Zuk, Bartosz Tuznik, Tadeusz Rymarz, Kamil Kwiatkowski,
	Marek Dudy?ski, Flavio C. C. Galeazzo, Guenther C. Krieger Filho
	Submitted to Computer Physics Communications
	
	

The solver is implemented in OpenFOAM 9 (openfoam.org), which can be obtained from repository:
https://github.com/mgzuber-eth/porousRedoxFoam
	
There is no dedicated user manual on this program. Please refer to the thesis 
of M. Zuber, the paper of P. Zuk, and the porousGasificationFoam files 
(https://github.com/btuznik/porousGasificationFoam). In brief, to use the solver, save 
the solver and src files in respective user directories, and compile them via 'wmake'.
	
Simulation case files are included in the run directory. Run the simulations via './Allrun'.
The case files for a variation of reduction and oxidation at Tsp=966C are included. 
	
Note that the solver is primarily a thermodynamic solver for solving solid-gas 
heterogeneous reactions. The gas flows through a porous fixed-bed. Apparent kinetics
are implemented. The majority of the thermodynamics are hard-coded into the src files.
