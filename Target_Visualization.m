
%============== Making a 3D nanoparticle with cubic voxels ===============%
clear all
clc
tic

%============================ Input variables  ===========================%
fprintf('\nEffective radius is defined as the radius of a sphere per equal volume of particle.');
r_eff=input('\nEnter the effective radius of the particle in nm:');
clc
d=input('\nEnter the cubical voxel size in nm:');
clc
fprintf('\nSelecting shape of the target');
fprintf('\nThe package can calculate optical properties for 4 different shapes:');
fprintf('\n1. spherical.');
fprintf('\n2. ellipsoid.');
fprintf('\n3. rod with semi-sphere caps.');
fprintf('\n4. rod-nanoshell.');
Np_shape=input('\n\nType the name of shape inside double quotation, \nEx. "spherical"  or "ellipsoid " or "rod" or "nanoshell_rod" :'); 
clc
dx=d;          % x-length of each cube
dy=d;          % y-length of each cube
dz=d;          % z-length of each cube
%=========================================================================%


%================== Coordinate of the center of Np =======================%
%=========================================================================%
X0=0;
Y0=0;
Z0=0;
%=========================================================================%

if Np_shape=="spherical"       % Nps is Sphere
    d_eff=2*r_eff;     % Effective diameter
    Lx=d_eff;
    Ly=d_eff;
    Lz=d_eff;
    
elseif Np_shape=="ellipsoid"   % Nps is Ellipsoid
    ARyx=input('\n\nEnter the ratio of y-semi axis to x-semi axis:');
    ARzx=input('\nEnter the ratio of z-semi axis to x-semi axis:');
    %%%%% if a=b and a<c Nps is prolate (incident light in z-direction) %%%
    %%%%% if a=b and a>c Nps is oblate  (incident light in z-direction) %%%
    d_eff=2*r_eff;
    a=((1/(ARyx*ARzx))^(1/3))*r_eff;    % Semi- minor axis in x- direction
    b=ARyx*a;                           % Semi- minor axis in y- direction
    c=ARzx*a;                           % Semi- major axis in z- direction
    
    Lx=2*a;
    Ly=2*b;
    Lz=2*c;
    
elseif Np_shape=="rod"  % Np is Rod with caps
    d_eff=2*r_eff;
    volume=4*pi/3*(r_eff^3);
    AR=2;         % Aspect Ratio: ratio of the long axis to short one
    r=(volume/(pi*(2*(AR-1)+4/3)))^(1/3); % Raduis of the Rod
    diameter=2*r;                         % Diameter of the Rod
    higth=AR*diameter;                    % Higth of the Rod
    
    Lx=diameter;
    Ly=diameter;
    Lz=higth;
    
elseif Np_shape=="nanoshell_rod"%Np is nanoshell (Rod with caps)
    d_eff=2*r_eff;
    volume=4*pi/3*(r_eff^3);
    ARyx=input('\n\nEnter the ratio of y-semi axis to x-semi axis:');
    ARzx=input('\nEnter the ratio of z-semi axis to x-semi axis:');
    a=(volume/(ARyx*ARzx))^(1/3);     % First side of rectangular block
    b=ARyx*a;                        % Second side of rectangular block
    c=ARzx*a;                        % Third side of rectangular block
    
    Lx=a;
    Ly=b;
    Lz=c;
end

%=========================================================================%


%===== Obtaining Coordinates of nanocubs or nanocells within the NPS =====%
%=========================================================================%
Max_x=Lx;
Max_y=Ly;
Max_z=Lz;



ix=-round(Max_x/(2*dx)):round(Max_x/(2*dx));
iy=-round(Max_y/(2*dy)):round(Max_y/(2*dy));
iz=-round(Max_z/(2*dz)):round(Max_z/(2*dz));
[y,x,z]=meshgrid(iy,ix,iz);
Nx=length(ix);                   % Number of the dipoles in the x-direction 
Ny=length(iy);                   % Number of the dipoles in the y-direction 
Nz=length(iz);                   % Number of the dipoles in the z-direction

N=Nx*Ny*Nz;                      % Total number of the dipoles 

% Converting each component of the dipoles coordinates to a vector 
X=(reshape(x,[N,1]))*dx;         % X-coordinates of the dipoles
Y=(reshape(y,[N,1]))*dy;         % Y-coordinates of the dipoles
Z=(reshape(z,[N,1]))*dz;         % Z-coordinates of the dipoles

r_block=[X Y Z];                 % Position of the each dipoles inside Nps

%=========================================================================%


%================ Finding Dipoles inside the nanoparticle ================%
%=========================================================================%
if Np_shape=="spherical"      %Nps is a Sphere
    Index_in=find(sqrt((X-X0).^2+(Y-Y0).^2+(Z-Z0).^2)<=(Lx/2-d/2));
    
elseif Np_shape=="ellipsoid"  %Nps is an Ellipsoid, vertically oriented
    Index_in=find(sqrt((X-X0).^2/((Lx/2)^2)+(Y-Y0).^2/((Ly/2)^2)+(Z-Z0).^2/((Lz/2)^2))<=1);
    
elseif Np_shape=="rod"  %Nps is Rod with caps, vertically oriented
    Index_in=find((sqrt((X-X0).^2+(Y-Y0).^2+(abs(Z-Z0)-Lz/2+Lx/2).^2)<=(Lx/2)...
        & abs(Z-Z0)>(Lz/2-Lx/2))|((sqrt((X-X0).^2+(Y-Y0).^2)<=(Lx/2))...
        & abs(Z-Z0)<=(Lz/2-Lx/2)));
    
elseif Np_shape=="nanoshell_rod"                %Nps is Rod-nanoshell with caps, vertically oriented
    % Drawing the cubes of a nanoshell, remove one quadratant of shape
    shell=input('\nEnter shell thickness :');
    Lx_core=Lx-shell;          % Diameter of the core region
    Ly_core=Ly-shell;          % Diameter of the core region
    Lz_core=Lz-shell;          % Hight of the core region
    
    % Finding dipoles inside the core region
    Index_core=find((sqrt((X-X0).^2+(Y-Y0).^2+(abs(Z-Z0)-Lz_core/2+Lx_core/2).^2)<=(Lx_core/2)...
        & abs(Z-Z0)>(Lz_core/2-Lx_core/2))|((sqrt((X-X0).^2+(Y-Y0).^2)<=(Lx_core/2))...
        & abs(Z-Z0)<=(Lz_core/2-Lx_core/2)));
    % Finding dipoles inside the whole particle
    Index_total=find((sqrt((X-X0).^2+(Y-Y0).^2+(abs(Z-Z0)-Lz/2+Lx/2).^2)<=(Lx/2)...
        & abs(Z-Z0)>(Lz/2-Lx/2))|((sqrt((X-X0).^2+(Y-Y0).^2)<=(Lx/2))...
        & abs(Z-Z0)<=(Lz/2-Lx/2)));
    
    % Finding dipoles inside the first quadrant of the core region
    Index_core_quad=find((X>=0 & Y>=0 & Z>=0 & (sqrt((X-X0).^2+(Y-Y0).^2+...
        (abs(Z-Z0)-Lz_core/2+Lx_core/2).^2)<=(Lx_core/2) ...
        & abs(Z-Z0)>(Lz_core/2-Lx_core/2)))| ...
        (X>=0 & Y>=0 & Z>=0 &(sqrt((X-X0).^2+(Y-Y0).^2)<=(Lx_core/2))...
        & abs(Z-Z0)<=(Lz_core/2-Lx_core/2)));
    
    % Finding dipoles inside the first quadrant of the whole Np
    Index_total_quad=find((X>=0 & Y>=0 & Z>=0 &(sqrt((X-X0).^2+(Y-Y0).^2+...
        (abs(Z-Z0)-Lz/2+Lx/2).^2)<=(Lx/2) & abs(Z-Z0)>(Lz/2-Lx/2)))| ...
        (X>=0 & Y>=0 & Z>=0 &(sqrt((X-X0).^2+(Y-Y0).^2)<=(Lx/2))& abs(Z-Z0)<=(Lz/2-Lx/2)));
    
    % Removing the contribution of the core region
    Multiply_core=zeros(N,1);
    Multiply_total=zeros(N,1);
    Multiply_core(Index_core,1)=Index_core;
    Multiply_core(Index_core_quad,1)=0;
    
    Multiply_total(Index_total,1)=Index_total;
    Multiply_total(Index_total_quad,1)=0;
    
    INDEX_Shell=Multiply_total-Multiply_core;
    
    Index_in=nonzeros(INDEX_Shell);
    Index_core=nonzeros(Multiply_core);
end
%=========================================================================%

X_in=X(Index_in,1);
Y_in=Y(Index_in,1);
Z_in=Z(Index_in,1);

for i=1:length(X_in)
    % Obtaining coordinate of the corners of the each nanocube
    nodeCoordinates=[0+X_in(i), 0+Y_in(i), 0+Z_in(i); 0+X_in(i), 0+Y_in(i), d+Z_in(i);...
        d+X_in(i), 0+Y_in(i), d+Z_in(i); d+X_in(i), 0+Y_in(i), 0+Z_in(i);...
        0+X_in(i), d+Y_in(i), 0+Z_in(i); 0+X_in(i), d+Y_in(i), d+Z_in(i);...
        d+X_in(i), d+Y_in(i), d+Z_in(i); d+X_in(i), d+Y_in(i), 0+Z_in(i)];
    
    elementNodes = [1 4 3 2; 5 8 7 6; 1 2 6 5; 3 4 8 7; 2 3 7 6; 1 5 8 4];
    patch('Faces', elementNodes,'EdgeColor','g','FaceColor','y', 'Vertices', nodeCoordinates)
    light               % create a light
    lighting gouraud    % preferred method for lighting curved surfaces
    ii=i
    axis equal
    hold on
end
