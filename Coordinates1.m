
function [Max_x,Max_y,Max_z,N,Nx,Ny,Nz,r_block,X,Y,Z,d_inter]=Coordinates1(Meshing,GPU,d,Lx,Ly,Lz,...
                                                         d_eff,Structure,arrangement)
dx=d;
dy=d;
dz=d;

if Structure=="monomeric" 
    d_inter=0;
    Max_x=Lx;
    Max_y=Ly;
    Max_z=Lz;
elseif Structure=="dimeric"           
    fprintf('\nA dimeric structure has been choosen.');
    d_ratio=input('\nEnter the ratio of (interparticle distance)/(effective diameter):');
    clc
    d_inter=d_ratio*d_eff;
                % For nanoparticles with head-tail orientation in z-direction
    if arrangement=="x_orient"  % side by side arrangement in x-direction
        Max_x=2*Lx+d_inter;
        Max_y=Ly;
        Max_z=Lz;
    elseif arrangement=="y_orient" % side by side arrangement in x-direction
        Max_x=Lx;
        Max_y=2*Ly+d_inter;
        Max_z=Lz;
    elseif arrangement=="z_orient"
        Max_x=Lx;    % head to tail arrangement in z-direction
        Max_y=Ly;
        Max_z=2*Lz+d_inter;
    end            
    
end

%======= Obtaining Coordinates of nanocubs or nanocells within NPS =======%
%=========================================================================%
if Meshing==1 % Standard Meshing
    
    if GPU==1
        ix=gpuArray(-round(Max_x/(2*dx))+0.5:round(Max_x/(2*dx))-0.5);
        iy=gpuArray(-round(Max_y/(2*dy))+0.5:round(Max_y/(2*dy))-0.5);
        iz=gpuArray(-round(Max_z/(2*dz))+0.5:round(Max_z/(2*dz))-0.5);
    elseif GPU==0
        ix=-round(Max_x/(2*dx))+0.5:round(Max_x/(2*dx))-0.5;
        iy=-round(Max_y/(2*dy))+0.5:round(Max_y/(2*dy))-0.5;
        iz=-round(Max_z/(2*dz))+0.5:round(Max_z/(2*dz))-0.5;
    end
    
elseif Meshing==2 % Nonstandard Meshing
    if GPU==1
        ix=gpuArray(-round(Max_x/(2*dx)):round(Max_x/(2*dx)));
        iy=gpuArray(-round(Max_y/(2*dy)):round(Max_y/(2*dy)));
        iz=gpuArray(-round(Max_z/(2*dz)):round(Max_z/(2*dz)));
    elseif GPU==0
        
        ix=-round(Max_x/(2*dx)):round(Max_x/(2*dx));
        iy=-round(Max_y/(2*dy)):round(Max_y/(2*dy));
        iz=-round(Max_z/(2*dz)):round(Max_z/(2*dz));
    end
end
[y,x,z]=meshgrid(iy,ix,iz);
Nx=length(ix);                   % Number of the dipoles in the x-direction 
Ny=length(iy);                   % Number of the dipoles in the y-direction 
Nz=length(iz);                   % Number of the dipoles in the z-direction 

N=Nx*Ny*Nz;

X=reshape(x,[N,1])*dx; % X-coordinates of the dipoles
Y=reshape(y,[N,1])*dy; % Y-coordinates of the dipoles
Z=reshape(z,[N,1])*dz; % Z-coordinates of the dipoles

r_block=[X Y Z];         % Position of the each dipoles inside rectangular block
%=========================================================================%

end