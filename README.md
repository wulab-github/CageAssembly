# CageAssembly
# Rigid Body based diffusion reaction simulation for protein cage assembly

FORTRAN source code:

CageAssemblyAlignFlex_main.f>>            main program of the diffusion reaction algorithm

diag.f>>                                  subroutine of matrix diagonization

rs.f>>                                    subroutines for the matrix operations

rotatefit.f>>                             subroutine for the rigid-body superposition



INPUT files:

CageKineticParameter_34.dat>>             sample input of the interaction rate parameters

CAGE.pdb>>                                the template of a fully assembled aprotein cage

OUTPUT files:

CageAssAlignFlex_VL100_output.dat>>       sample output of simulation record for size of assembled oligomers

CageAssAlignFlex_traj_VL100.pdb>>         sample output of simulation trjaectories in PDB format
