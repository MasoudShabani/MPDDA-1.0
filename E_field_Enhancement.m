%=========================================================================%
%==== Calculating enhanced E-field at monomeric and dimeric structures ===%
%=========================================================================%
% The package is generic and is able to calculate electric field around and
% ... inside an excited plasmonic NP. The package can calculate E-field for
% ... spherical, ellipsoid, rod shape, and rectangular block in monomeric and
% ... dimeric structure.

% "GpuArray" will transfer data from CPU to GPU
% "gather" will transfer data from GPU to CPU
%=========================================================================%
clear all
clc

CT=10^(-5);                     % Convergence threshold value in interative solver

%============ LSPR wavelength and bulk RI of metal at LSPR ===============%
%=========================================================================%
fprintf('\nThe package can be run both in CPU and GPU.');
GPU=input('\nEnter 1 for running in GPU, and 0 for running in CPU:');
clc
Lambda=input('\nEnter the incident light wavelength:');
Re_n=input('\nEnter the real part of bulk refractive index of target:');
Im_n=input('\nEnter the imaginary part of bulk refractive index of target:');
clc

eps=(Re_n+1i*Im_n).^2;
Re_eps=real(eps);
Im_eps=imag(eps);
%=========================================================================%


%========================== Input Parameters =============================%
%=========================================================================%
nb=input('\nEnter the refractive index of the surrounding medium:');
epsb=nb^2;         % Dielectric function of background medium
clc
fprintf('\nEffective radius is defined as the radius of a sphere per equal volume of particle.');
r_eff=input('\nEnter the effective radius of the particle in nm:');
d_eff=2*r_eff;
volume=4*pi/3*(r_eff^3);
clc
d=input('\nEnter the cubical voxel size in nm:');
clc
fprintf('\nSelecting shape of the target');
fprintf('\nThe package can E-field around and inside of4 different shapes:');
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
W=2*pi*3*10^17./Lambda;         % incident radiation frequency
k=2*pi./Lambda*nb;              % wave vector of light in first layer

%========= Modifying dielectric function by adding size effect  ==========%
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
Refractive=eps_nps.^(0.5);
Re_Refractive=real(Refractive);
Im_Refractive=imag(Refractive);

if GPU==1
    ep_nps_eb=gpuArray(eps_nps./epsb);  % Ratio of metal-to-medium dielectric function
elseif GPU==0
    ep_nps_eb=eps_nps./epsb;  % Ratio of metal-to-medium dielectric function
end
%=========================================================================%


%====== Finding dimension of the Rectangular block emcopasses NPs ========%
%=========================================================================%
[Lx,Ly,Lz]= Nps_parameters(r_eff,Np_shape);
%=========================================================================%

%============ Selectring the 2D plane for calculating E-field =============%
fprintf('\nChoosing the 2D plane to calculate E-field.\n1. xy plane.\n2.xz plane. \n3. yz plane.');
plane_2D=input('\nEnter 2D plane, Ex. "xy" or "xz" or "yz" :');
clc
%=========================================================================%
%==== Finding coordinate of the dipoles inside the rectangular block =====%
%=========================================================================%
[Max_x,Max_y,Max_z,N,Nx,Ny,Nz,r_block,X,Y,Z,d_inter]=Coordinates(GPU,d,Lx,Ly,Lz,...
    d_eff,Structure,arrangement); 
Nx_target=Nx;
Ny_target=Ny;
Nz_target=Nz;
Xrectan=X;
Yrectan=Y;
Zrectan=Z;
%=========================================================================%

%===================== Selecting the incident light ======================%
%=========================================================================%
fprintf('\nChoosing polarization and propogation direction of the incident light:');
E01=input('\n\nEnter polarization direction, \nex. [0 0 1] for z direction:');
K01=input('\n\nEnter propogation direction, \nex. [1 0 0] for x direction:');
%==============%
if GPU==1 
    E0=gpuArray(E01);
    K0=gpuArray(K01);   
elseif GPU==0
    E0=E01;               % Incident electric field
    K0=K01;               % unit vector in direction of wave vector
end
%=============================%
clc
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
 
%=========== Printing general information reqarding the problem ==========%
tic    % Start time of the calculation
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
%=========================================================================%


%===== Finding index of dipoles inside NPs & ingoring other elements =====%
%=========================================================================%
[INDEX_INSIDE]=INDEX_INSIDE_NP(GPU,X,Y,Z,N,Np_shape,Lx,Ly,Lz,Structure,d_inter,arrangement);

INDEX_IN=reshape(INDEX_INSIDE,[Nx,Ny,Nz]);% INDEX_IN contains zeros and ones, 
                  % ... zeros for dipole outside of NP, and ones for inside                                     
%=========================================================================%


%===== Calculating RijRij-I3 and 3RijRij-I3 in interaction matrix A ======%
%=========================================================================%
[rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,rjkrjk31_I,rjkrjk32_I,...
    rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,RJK]=RijRij(r_block);
%=========================================================================%


%====== Obtaining inverse of polarizibilty of each nanocube at LSPR ======%
%=========================================================================%
eps_NP_eb=ep_nps_eb;
kvec=k*K0;
[Inverse_Alpha]=Polarizability(GPU,kvec,eps_NP_eb,INDEX_IN,d,E0);
%=========================================================================%


%============== Calculating Componests of Incident E-filed ===============%
%=========================================================================%
[E_x,E_y,E_z,E_vector]=Incident_Field(Lambda,IB,nb,r_block,kvec,K0,INDEX_INSIDE,Nx,Ny,Nz,E0,z0,Waist_r);
%=========================================================================%


%======= Calculating 6 tensor blocks: Axx, Axy, Axz, Ayy, Ayz & Azz ======%
%=========================================================================%
Exp_ikvec_rjk=exp(1i*norm(kvec)*RJK)./RJK;
ikvec_rjk=(1i*norm(kvec)*RJK-1)./(RJK.^2);   % ikvec_rjk=(1i*norm(kvec)*rjk-1)/rjk^2

[Axx,Axy,Axz,Ayy,Ayz,Azz]=Interaction_Matrix(kvec,Exp_ikvec_rjk,...
    ikvec_rjk, rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,...
    rjkrjk31_I,rjkrjk32_I,rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,Nx,Ny,Nz);
%=========================================================================%


%===== Calculating FFT of six tensor blocks of the interaction matrix ====%
%=========================================================================%
[FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ]=FFT_Interaction(GPU,Axx...
    ,Axy,Axz,Ayy,Ayz,Azz,Nx,Ny,Nz);
%=========================================================================%


              % Biconjugate Gradient(BCG), iterative method %
%============= Applying BCG to obtain Px,Py,Pz by iteration ==============%
%=========================================================================%
[px,py,pz]=Biconjugate_Gradient(E_x,E_y,E_z,Nx,Ny,Nz,N,Inverse_Alpha,...
    INDEX_IN,E_vector,FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ,CT);
%=========================================================================%


                   % Deleting unnessesary data %
%=========================================================================%                
clear FFT_AXX FFT_AXY FFT_AXZ FFT_AYY FFT_AYZ FFT_AZZ 
%=========================================================================%
    
    
             % Ignoring polarizibality of dipole outside NPs %
%=========================================================================%        
px=px.*INDEX_IN;
py=py.*INDEX_IN;
pz=pz.*INDEX_IN;

PX_vector=reshape(px,[N,1]);
PY_vector=reshape(py,[N,1]);
PZ_vector=reshape(pz,[N,1]);
P_vector=[PX_vector;PY_vector;PZ_vector];
%=========================================================================%   

                    % Deleting unnessesary data %
%=========================================================================%  
clear  a_CM_Nps anr_Nps aLDR_Nps  a_CM_Matrix anr_Matrix aLDR_Matrix ...
    Ex Ey Ez E_x E_y E_z ...
    E_vector Exp_ikvec_rjk ikvec_rjk r_block 
%=========================================================================% 




%=========================================================================%
%======== E-field outside of the NPs can be calculated in any plane ======% 
% Ex. I am calculating E-field in xz plane at y=0
% The electric field outside of the NPs has been calculated
%=========================================================================% 

%= Obtaining coordinate of voxels inside & outside of NPs in extended system =%
%=========================================================================%
[r_block,r_plane,Nx,Ny,Nz,N_plane,x_plane,y_plane,z_plane,SIZE,Xex,Yex,Zex]=...
    Two_D_plane(GPU,d,Max_x,Max_y,Max_z,Lx,Ly,Lz,plane_2D);
%=========================================================================%

%==== Obtaining E_incident at position of each dipole on the 2D plane ====%
%=========================================================================%
if IB=="plane wave"    % Incident beam is plane wave
    kr=kvec(1)*r_plane(:,1)+kvec(2)*r_plane(:,2)+kvec(3)*r_plane(:,3);
    expikr=exp(1i*kr);
    
    Ex_incident=(E0(1)*expikr);
    Ey_incident=(E0(2)*expikr);
    Ez_incident=(E0(3)*expikr);
elseif IB=="gaussian"   % Incident beam is Gaussian
    
    if K0(1)==1
        r=(r_plane(:,2).^2+r_plane(:,3).^2).^(0.5);
        z=r_plane(:,1);
    elseif K0(2)==1
        r=(r_plane(:,1).^2+r_plane(:,3).^2).^(0.5);
        z=r_plane(:,2);
    else
        r=(r_plane(:,1).^2+r_plane(:,2).^2).^(0.5);
        z=r_plane(:,3);
    end
  
    W0=Waist_r*Lambda;                    %Waist radius of the beam
    K=2*pi*nb/Lambda;                     %wave number
    zR=pi*(W0.^2).*nb./Lambda;
    Qz=atan((z-z0)./zR);                    % Is the Gouy phase at z
    Wz=W0.*((1+((z-z0)./zR).^2).^0.5);      % the spot size parameter
    Rz_inverse=z./((z-z0).^2+zR.^2);        % inverse of radius of curvature
    
    E=(W0./Wz).*exp(-(r./Wz).^2-1i*(K.*(z-z0)+K.*(r.^2).*Rz_inverse/2-Qz));
    Ex_incident=E0(1)*E;
    Ey_incident=E0(2)*E;
    Ez_incident=E0(3)*E;
end
%=========================================================================%


%============ Calculating total E-field Outside of the traget ============%
%=========================================================================%
[Ex_out,Ey_out,Ez_out]=E_total(GPU,N_plane,Nx,Ny,Nz,INDEX_IN,Inverse_Alpha,Nx_target,...
    Ny_target,Nz_target,x_plane,y_plane,z_plane,r_block,Np_shape,plane_2D,...
    Lx,Ly,Lz,kvec,d_inter,Structure,Ex_incident,Ey_incident,Ez_incident,px,py,pz,arrangement);

E_out=sqrt(Ex_out.*(Ex_out)+Ey_out.*(Ey_out)+Ez_out.*(Ez_out));

%medfilt2(I,[m n]) performs median filtering, where each output pixel contains the median value in the m-by-n neighborhood around the corresponding pixel in the input image.
%Function medfilt2 discards points that differ considerably from their surroundings ...
%... and replaces them by the median of that point and a specified number of neighboring points. 
%=========================================================================%
Et_out=abs(reshape(E_out,[SIZE(1,1),SIZE(1,2)]));
if GPU==0 
    % GPU can not support medfilt2
    Et_out=medfilt2(Et_out,[5 5]);
elseif GPU==1
    fprintf('\GPU can not support medfilt2 syntax for removing spike in surface plot.');
end
I_out=Et_out.^2;
I_out=gather(I_out);
%=========================================================================%


%============== Calculating total E-field Inside the target ==============%
%=========================================================================%

if plane_2D=="xy"
    Index_in=find(Zrectan==0);
    Xex=gather(Xex);
    Yex=gather(Yex);
elseif plane_2D=="xz"
    Index_in=find(Yrectan==0);
    Xex=gather(Xex);
    Zex=gather(Zex);
elseif plane_2D=="yz"
    Index_in=find(Xrectan==0);
    Yex=gather(Yex);
    Zex=gather(Zex);
end

Xtarget=Xrectan(Index_in);
Ytarget=Yrectan(Index_in);
Ztarget=Zrectan(Index_in);
Px_in=PX_vector(Index_in);
Py_in=PY_vector(Index_in);
Pz_in=PZ_vector(Index_in);
Inverse_Polar=Inverse_Alpha(Index_in);

Ex_in=Inverse_Polar.*Px_in;
Ey_in=Inverse_Polar.*Py_in;
Ez_in=Inverse_Polar.*Pz_in;
E_in=(Ex_in.*(Ex_in)+Ey_in.*(Ey_in)+Ez_in.*(Ez_in)).^(0.5);
%=========================================================================%


%================= Surface plot of the electric field ====================%
%=========================================================================%
if plane_2D=="xy"
    X_in=reshape(Xtarget,[Nx_target,Ny_target]);
    Y_in=reshape(Ytarget,[Nx_target,Ny_target]);
    Et_in=abs(reshape(E_in,[Nx_target,Ny_target]));
    I_in=Et_in.^2;
    surf(Xex,Yex,Et_out,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    shading interp;
    % waitforbuttonpress
    axis equal
    hold on
%     surf(X_in,Y_in,Et_in,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
%     shading interp;
%     %waitforbuttonpress
%     axis equal
elseif plane_2D=="xz"
    X_in=reshape(Xtarget,[Nx_target,Nz_target]);
    Z_in=reshape(Ztarget,[Nx_target,Nz_target]);
    Et_in=abs(reshape(E_in,[Nx_target,Nz_target]));
    I_in=Et_in.^2;
    surf(Xex,Zex,Et_out,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    shading interp;
    % waitforbuttonpress
    axis equal
    hold on
%     surf(X_in,Z_in,Et_in,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
%     shading interp;
    %waitforbuttonpress
%     axis equal
    
elseif plane_2D=="yz"
    Y_in=reshape(Ytarget,[Ny_target,Nz_target]);
    Z_in=reshape(Ztarget,[Ny_target,Nz_target]);
    Et_in=abs(reshape(E_in,[Ny_target,Nz_target]));
    I_in=Et_in.^2;
    surf(Yex,Zex,Et_out,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    shading interp;
    % waitforbuttonpress
    axis equal
%     hold on
%     surf(Y_in,Z_in,Et_in,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
%     shading interp;
%     waitforbuttonpress
%     axis equal
end

Time=toc;
%=========================================================================%
fprintf('\nElectric field calculation has been finished.\n\n');
fprintf('\nThe electric field has been calculated in :');
fprintf(plane_2D);
fprintf('\nTotal time for calculating electric field in SECONDS:');
disp(Time);
fprintf('\nCritririon for stopping iterative solver,refer line 15:');
disp(CT);
fprintf('\nTHE RESULTS CAN BE FOUND IN  Result_E_enhanced.mat FILE');
%=========================================================================%
save Result_E_enhanced Time Lambda X Y Z Et_out Et_in I_out I_in 
