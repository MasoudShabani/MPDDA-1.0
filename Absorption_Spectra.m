%=========================================================================%
%============= Calculating Qext, Qabs and Qscat using DDA ================%
%=========================================================================%
% The package is able to calculate extinction, absorption and scattering
% ... efficiencies of 4 different shape Nps for monomeric and dimeric
% ...structures
% Qext= Extinction efficiency
% Qabs= Absorption efficiency
% Qscat= Scattering efficiency
% "GpuArray" will transfer data from CPU to GPU
% "gather" will transfer data from GPU to CPU
%=========================================================================%

clear all
clc
%========= loading Wavelength & Bulk Refractive index of the metal =======%
%=========================================================================%
% if input data are saved in a excell sheet it can load as below
%Data=xlsread('/Users/masoudshabaninezhadnavrood/Desktop/DDA All, 09.15.2020/MPDDA- 1.0 copy/Au_Bulk_RI.xlsx');
% if input data are m.file, it can be loaded: ...
Data=load('Copy the address link of the initial data here');

Wavelength=Data(:,1);
Re_n=Data(:,2);                 % Real part of the refractive index
Im_n=Data(:,3);                 % Imaginary part of the refractive index

eps=(Re_n+1i*Im_n).^2;          % Dielectric function of the bulk metal

CT=10^(-5);                     % Convergence threshold value in interative solver

%============================== Input Parameters =========================%
%=========================================================================%

fprintf('\nThe package can be run both in CPU and GPU.');
GPU=input('\nEnter 1 for running in GPU, and 0 for running in CPU:');
clc
nb=input('\nEnter the refractive index of the surrounding medium:');
clc
fprintf('\nEffective radius is defined as the radius of a sphere per equal volume of particle.');
r_eff=input('\nEnter the effective radius of the particle in nm:');
clc
d=input('\nEnter the cubical voxel size in nm:');
clc
fprintf('\nSelecting shape of the target');
fprintf('\nThe package can calculate optical properties for 4 different shapes:');
fprintf('\n1. spherical.');
fprintf('\n2. ellipsoid, which is oriented in z direction.');
fprintf('\n3. rod with semi-sphere caps, which is oriented in z direction.');
fprintf('\n4. rectangular block, which is oriented in z direction.');
Np_shape=input('\n\nType the name of shape inside double quotation, \nEx. "spherical"  or "ellipsoid " or "rod" or "rec_block" :'); 
clc
fprintf('\nSelecting structure of the problem');
fprintf('\nThe package can calculate optical properties for 2 different structure:');
fprintf('\n1. monomeric.');
fprintf('\n2. dimeric.');
Structure=input('\n\nType the name of structure inside double quotation, \nEx. "monomeric" or "dimeric" :');
clc
if Structure=="dimeric"
    fprintf('\nThree kind of arragments for neighboring NPs:');
    fprintf('\n1. For head to tail orientation in z-direction, type "z_orient"  \n');
    fprintf('2. For side by side orientation in x-direction, type "x_orient"  \n');
    fprintf('3. For side by side orientation in y-direction,type "y_orient"  \n');
    arrangement=input('\nType the orientation of the neighboring NPs:');
    clc
else
     arrangement="monomeric";
end
clc
r_np=r_eff*(10^(-9));               % effective radius of the nanoparticle
W=2*pi*3*10^17./Wavelength;         % incident radiation frequency
k=2*pi./Wavelength*nb;              % wave vector of light in first layer


fprintf('\nChoosing polarization and propogation direction of the incident light:');
E01=input('\nEnter polarization direction, \nex. [0 0 1] for z direction:');
K01=input('\nEnter propogation direction, \nex. [1 0 0] for x direction:');
clc

%=============== Obtaining Modified dielectric function ==================%
%=========================================================================%
fprintf('\nModifying the electric permitivity of NPs \nby considering size effec in collision frequency');
Np_name=input('Choose the name of the Np,\nEx. "Au" or "Ag" or "Cu" or "other" :');
clc
if Np_name=="Au"
    Wp=8.9*1.5186*(10^15);              % plasma frequency of Au
    L0=0.07/6.58211951440*10^(-16);     % Collision freguency of Au in bulk medium
    Ap=0.5;                             % damping correction factor
    Vf=1.4*(10^6);                      % Fermi Velocity
elseif Np_name=="Ag"
    L0=3.22*(10^13);  % Collision freguency of Ag in bulk medium or bulk damping constant, Jonhson Paper 1972
    Vf=1.39*(10^6);    % Fermi velocity of Ag atoms in bulk medium
    Ap=0.25;
    Wp=1.393*(10^16); %Jonhson Paper 197
elseif Np_name=="Cu"
    L0=1.45*(10^14);  % Collision freguency of Cu in bulk medium or bulk damping constant, Jonhson Paper 1972
    Ap=0.5;
    Wp=1.344*(10^16); %Jonhson Paper 1972
    Vf=1.59*(10^6);    % Fermi velocity of Cu atoms in bulk medium
else
    Wp=input('\nEnter the plasma frequency of the selected Np:');
    L0=input('\nEnter the bulk collision frequency of the target:');
    Vf=input('\nEnter the Fermi Velocity of the target:');
    Ap=input('\nEnter damping correction factor, \nit is usually between 0.2 till 2:');
    clc
end
L=L0+Ap*Vf/r_np;                         % Modified damping frequency
eps_nps=eps+(Wp^2)./(W.^2+1i*L0.*W)-(Wp.^2)./(W.^2+1i*W.*L); % Modified dielectric function


%====== Finding dimension of the Rectangular block emcopasses NPs ========%
%=========================================================================%
[Lx,Ly,Lz]= Nps_parameters(r_eff,Np_shape);
%=========================================================================%
d_eff=2*r_eff;
volume=4*pi/3*(r_eff^3);
epsb=nb^2;                % Dielectric function of background medium

if GPU==1 
    E0=gpuArray(E01);
    K0=gpuArray(K01);   
elseif GPU==0
    E0=E01;               % Incident electric field
    K0=K01;               % unit vector in direction of wave vector
end
%=============================%
Refractive=eps_nps.^(0.5);
Re_Refractive=real(Refractive);
Im_Refractive=imag(Refractive);

if GPU==1
    ep_nps_eb=gpuArray(eps_nps./epsb);  % Ratio of metal-to-medium dielectric function
elseif GPU==0
    ep_nps_eb=eps_nps./epsb;  % Ratio of metal-to-medium dielectric function
end
%=========================================================================%
fprintf('\n choosing the Meshing: \n');
fprintf('\1. Standard Meshing: The particle size is set to be equal to the distance between outer edges of two opposite boundary dipoles\n\n');
fprintf('\2. Nonstandard Meshing:The particle size is set to be equal to the distance between the center of two opposite boundary dipoles\n');
Meshing=input('\n Type 1 if it is standard meshing, otherwise type 2:');

    % Finding coordinate of the dipoles inside the rectangular block %
[Max_x,Max_y,Max_z,N,Nx,Ny,Nz,r_block,X,Y,Z,d_inter]=Coordinates1(Meshing,GPU,d,Lx,Ly,Lz,...
    d_eff,Structure,arrangement); 

fprintf('\nSelecting the incident light:')
fprintf('\nThe package can support two incident light:');
fprintf('\n1. plane wave.');
fprintf('\n2. gaussian.');
IB=input('\n\nEnter the name of incident beam inside double quotation, \nEx. "plane wave" or "gaussian" :'); 
clc
if IB=="plane wave" 
    z0=0;              % Focus point of the Gaussian beam
    Waist_r=100;       % ratio of waist raduis of beam to wavelength
elseif IB=="gaussian"
   fprintf('\nA Gaussian beam has been choosen.');
   z0=input('\nEnter focus point of beam, \nEx. could be -Max_z/2 (Max_z is the length of structure in z-direction, center is at 0):');
   clc
   fprintf('\nChoosing the waist raduis of beam:');
   fprintf('\nIf the waist raduis is much bigger than the NPs size, the results will be like plane wave results');
   Waist_r=input('\nEnter the RATIO of waist raduis of beam to wavelength, \nEx. 0.1:');
   clc
end
 
tic    % Start time of the calculation
%Finding index of dipoles inside NPs & ignoring outside ones%
fprintf('\nCalculation has been started ...');
fprintf('\n\nName of the Nps:');
fprintf(Np_name);
fprintf('\nNanoparticle Shape:');
fprintf(Np_shape);
fprintf('\nStructure:');
fprintf(Structure);
fprintf('\nTotal number of voxels:');
disp(N);
fprintf('Size of the cubical voxel:')
disp(d);
fprintf('Effective raduis of the NP:')
disp(r_eff);
fprintf('The refractive index of the medium:')
disp(nb);
fprintf('\nThe incident beam:');
fprintf(IB);

[INDEX_INSIDE]=INDEX_INSIDE_NP(GPU,X,Y,Z,N,Np_shape,Lx,Ly,Lz,Structure,d_inter,arrangement);

INDEX_IN=reshape(INDEX_INSIDE,[Nx,Ny,Nz]); % INDEX_IN contain zeros and ones, zeros for dipole outside of NP & one for inside  
%=========================================================================%


%======Calculating RijRij-I3 and 3RijRij-I3 in interaction matrix A=======%
%=========================================================================%
[rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,rjkrjk31_I,rjkrjk32_I,...
    rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,RJK]=RijRij(r_block);
%=========================================================================%
toc;
t0=toc;
t0_initial=t0;

counting2=0;
for J=1:length(Wavelength)
    
    eps_NP_eb=ep_nps_eb(J);
    Lambda=Wavelength(J);
    
    if GPU==1
        kvec=gpuArray(k(J)*K0);
    elseif GPU==0
        kvec=k(J)*K0;
    end
    
    
    %========= Calculating Incident electric filed components ============%
    %=====================================================================%
    [E_x,E_y,E_z,E_vector]=Incident_Field(Lambda,IB,nb,r_block,kvec,K0,INDEX_INSIDE,Nx,Ny,Nz,E0,z0,Waist_r);
    %=====================================================================%
    
    
    %======Obtaining inverse of P of each dipole at different Lambda======%
    %=====================================================================%
    [Inverse_Alpha]=Polarizability(GPU,kvec,eps_NP_eb,INDEX_IN,d,E0);
    %=====================================================================%
    
    
    %========= Calculating 6 tensor blocks of interaction matrix =========%
    %=====================================================================%
    Exp_ikvec_rjk=exp(1i*norm(kvec)*RJK)./RJK;
    ikvec_rjk=(1i*norm(kvec)*RJK-1)./(RJK.^2);   % ikvec_rjk=(1i*norm(kvec)*rjk-1)/rjk^2
    
    [Axx,Axy,Axz,Ayy,Ayz,Azz]=Interaction_Matrix(kvec,Exp_ikvec_rjk,...
        ikvec_rjk, rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,...
        rjkrjk31_I,rjkrjk32_I,rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,Nx,Ny,Nz);
    %=====================================================================%
    
    
    %===== Calculating FFT of six tensor blocks of interaction matrix ====%
    %=====================================================================%
    [FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ]=FFT_Interaction(GPU,Axx...
            ,Axy,Axz,Ayy,Ayz,Azz,Nx,Ny,Nz);
    %=====================================================================%
    
    
             % Iterative Method, Biconjugate gradient & inverse FFT%
    %=== Applying Biconjugate gradient & inverse FFT to obtainPx,Py,Pz ===% 
    %=====================================================================%
    [px,py,pz]=Biconjugate_Gradient(E_x,E_y,E_z,Nx,Ny,Nz,N,Inverse_Alpha,...
             INDEX_IN,E_vector,FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ,CT);
    %=====================================================================%
    
    
                    % Deleting unnessesary data %
    %=====================================================================%                
    clear FFT_AXX FFT_AXY FFT_AXZ FFT_AYY FFT_AYZ FFT_AZZ 
    %=====================================================================%
    
    
             % Ignoring polarizibality of dipole outside NPs %
    %=====================================================================%        
    px=px.*INDEX_IN;
    py=py.*INDEX_IN;
    pz=pz.*INDEX_IN;
    %=====================================================================%   
    
    PX_vector=reshape(px,[N,1]);
    PY_vector=reshape(py,[N,1]);
    PZ_vector=reshape(pz,[N,1]);
    Inv_Alpha=reshape(Inverse_Alpha,[N,1]);
    Inv_Alpha_vec=[Inv_Alpha;Inv_Alpha;Inv_Alpha];
    counting2=counting2+1;
    
    
    %======== Calculating Cabs, Cscat & Cext for each wavelength =========%
    %=====================================================================%
    P_vector=[PX_vector;PY_vector;PZ_vector];
    Cabs(J)=4*pi*norm(kvec)/sum(abs(E0.^2))*((imag(dot((conj(P_vector)),conj(P_vector.*Inv_Alpha_vec)))-2/3*norm(kvec)^3*(norm(P_vector).^2)));
    Cext(J)=4*pi*norm(kvec)/sum(abs(E0.^2))*imag(dot((E_vector),P_vector));
    Cscat(J)=Cext(J)-Cabs(J);
    %=====================================================================c
    Time_each(J)=toc-t0;
    t0=toc;
    
                    % Deleting unnessesary data %
    %=====================================================================% 
    clear PX_vector PY_vector PZ_vector P_vector px py pz a_CM_Nps anr_Nps ...
        aLDR_Nps  a_CM_Matrix anr_Matrix aLDR_Matrix Ex Ey Ez E_x E_y E_z ...
        E_vector Exp_ikvec_rjk ikvec_rjk Inverse_Alpha
    %=====================================================================%
 
end
Total_Time=toc;


%========== Calculating Qabs, Qscat & Qext for all wavelengths ===========%
%=========================================================================%
Q_ABS=gather(Cabs)/(pi*(r_eff^2));   % Transfering Data from GPU to CPU
Q_EXT=gather(Cext)/(pi*(r_eff^2));   % Transfering Data from GPU to CPU
Q_SCAT=gather(Cscat)/(pi*(r_eff^2)); % Transfering Data from GPU to CPU
%=========================================================================%


%======================== Plot & Save the Result =========================%
%=========================================================================%
plot(Wavelength,Q_EXT,Wavelength,Q_ABS,Wavelength,Q_SCAT)

fprintf('\nExtinction, Absorption and Scattering efficiencies have been calculated.\n\n');
fprintf('Total time for calculating whole spectra in SECONDS:');
disp(Total_Time);
fprintf('\nCritririon for stopping iterative solver,refer line 26:');
disp(CT);

% fprintf('\n\nI. Extinction efficiency:\n');
% disp(Q_EXT);
% fprintf('\n\nII. Absorption efficiency:\n');
% disp(Q_ABS);
% fprintf('\n\nIII. Scattering efficiency:\n');
% disp(Q_SCAT);
fprintf('\nTHE RESULTS CAN BE FOUND IN  Result_Spectra.mat FILE');

save Result_Spectra20 d N r_eff t0_initial Total_Time Time_each Wavelength  Q_EXT Q_ABS Q_SCAT 
%=========================================================================%



