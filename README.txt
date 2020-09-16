
                                                                    
                         Discrete Dipole Approximation (DDA) in MATLAB:


    This package developed in MATLAB software to describe the optical properties of the plasmonic nanoparticles using DDA. Accelerative methods (Fast Fourier Transform (FFT) and Biconjugate Gradient (BCG)) have been applied to this toolbox to reduce the computational cost. Also, this package has the option to use the parallel computing toolbox in MATLAB to reduce the time of calculation more by running the codes in GPU. However, the code can simply be executed in CPU as well.


%=============================================================================================================================================================================%
                           
                         Executing the codes:


A) System requirement to run the codes in GPU: Any system that is equipped with GPU hardware, MATLAB software, and supports GPU toolbox of the MATLAB. 

B) To execute the program, download all files and put all of the codes in one folder. 

C) To calculate absorption, scattering and extinction efficiencies, open the Absorption_Spectra.m file and run it. However, before running this code, you need to upload the bulk refractive index of the metal at different wavelengths (https://refractiveindex.info/ is an excellent reference to obtain the refractive index of bulk metal). You can simply save the bulk refractive indices in an excel sheet (first column wavelengths, second column real part of the refractive indices and third column imaginary part of the refractive indices) and copy the address link of the file in the line 19 of Absorption_Spectra.m file. As an example, we provided bulk refractive indices of Au in the repository.

D) To obtain electric field enhancement, open the E_field_Enhancement.m file and run it. It should be noted that for removing spikes in the surface plot we used midfilt2 syntax of MATLAB. However, only CPU mode of the toolbox can support this syntax.

E). To create different shape nanoparticles with nanocubic voxels, open the Target_Visualization.m file and run it.
  
%===============================================================================================================================================================================% 
 

                          Describing the files in the package:

1) Absortion_Spectra.m: simulates absorption, scattering, and extinction cross-sections of different shape plasmonic nanoparticles in monomeric and dimeric nanostructures. 

2) Biconjugate_Gradient.m: runs the iterative method (BCG) to modify the polarization of each dipole of the plasmonic nanostructures.

3) Coordinates.m: finds coordinate of each dipole inside the smallest rectangular block that encompasses the plasmonic nanostructures.

4) E_field_Enhancement.m: simulates the electric field enhancement around plasmonic nanostructures and plots its surface plot.

5) E_total.m: calculates the total electric field at the position of each dipole.

6) Excluding_NPs.m: excludes the coordinates of the nanocubes inside the plasmonic nanoparticle(s). This file obtains the coordinates of the nanocubes outside the  nanoparticle(s) boundary, which the electric field needs to be calculated.

7) FFT_Interaction.m: computes the 1D FFT of the six tensor blocks of the interaction matrix in the DDA.

8) Target_Visualization.m: Can be used to visualize different shape nanoparticles made of nanocubic voxels.

9) Incident_Field.m: calculates the incident light field at the position of each dipole inside the nanoparticle(s).

10) INDEX_INSIDE_NP.m: obtains coordinates of the nanocubic voxels inside nanoparticle(s).

11) Interaction_Matrix.m: computes the six tensor blocks of the interaction matrix in the DDA. 

12) Inverse_FFT.m: calculates the 1D FFT of the components of the polarization vector in each iteration, performs elementwise multiplication of the FFT of the six tensor blocks of the interaction matrix to the FFT of the three components of the polarization vector, and computes the inverse FFT (IFFT) of each.
 
13) Nps_parameters.m: obtains dimension of the nanoparticle(s).

14) Polarizability.m: calculates the polarizability of each nanocube.

15) RijRij.m: is required to calculate the interaction matrix and depends on the inter-dipoles distance.

16). Two_D_plane.m: chooses an extended 2D plane for calculating E-field.

17). Gaussian_Beam.m: Simulate a Gaussian beam. 


%==============================================================================================================================================================================%

                                   
 

