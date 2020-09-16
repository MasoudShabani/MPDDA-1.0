
%===== Calculating FFT of six tensor blocks of the interaction matrix=====%

function [FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ]=FFT_Interaction(GPU,...
    Axx,Axy,Axz,Ayy,Ayz,Azz,Nx,Ny,Nz)

if GPU==1
    AXX=zeros(2*Nx-1,2*Ny-1,2*Nz-1,'gpuArray');
    AXY=zeros(2*Nx-1,2*Ny-1,2*Nz-1,'gpuArray');
    AXZ=zeros(2*Nx-1,2*Ny-1,2*Nz-1,'gpuArray');
    AYY=zeros(2*Nx-1,2*Ny-1,2*Nz-1,'gpuArray');
    AYZ=zeros(2*Nx-1,2*Ny-1,2*Nz-1,'gpuArray');
    AZZ=zeros(2*Nx-1,2*Ny-1,2*Nz-1,'gpuArray');
elseif GPU==0
    AXX=zeros(2*Nx-1,2*Ny-1,2*Nz-1);
    AXY=zeros(2*Nx-1,2*Ny-1,2*Nz-1);
    AXZ=zeros(2*Nx-1,2*Ny-1,2*Nz-1);
    AYY=zeros(2*Nx-1,2*Ny-1,2*Nz-1);
    AYZ=zeros(2*Nx-1,2*Ny-1,2*Nz-1);
    AZZ=zeros(2*Nx-1,2*Ny-1,2*Nz-1);
end

AXX(1:Nx,1:Ny,1:Nz)=Axx;
AXY(1:Nx,1:Ny,1:Nz)=Axy;
AXZ(1:Nx,1:Ny,1:Nz)=Axz;
AYY(1:Nx,1:Ny,1:Nz)=Ayy;
AYZ(1:Nx,1:Ny,1:Nz)=Ayz;
AZZ(1:Nx,1:Ny,1:Nz)=Azz;


%======= Expanding dimension of each line in X-direction to 2Nx-1 ========%
%=========================================================================%
AXX(Nx+1:2*Nx-1,1:Ny,1:Nz)=Axx(Nx:-1:2,1:Ny,1:Nz);
clear Axx
AXY(Nx+1:2*Nx-1,1:Ny,1:Nz)=-Axy(Nx:-1:2,1:Ny,1:Nz);
clear Axy
AXZ(Nx+1:2*Nx-1,1:Ny,1:Nz)=-Axz(Nx:-1:2,1:Ny,1:Nz);
clear Axz
AYY(Nx+1:2*Nx-1,1:Ny,1:Nz)=Ayy(Nx:-1:2,1:Ny,1:Nz);
clear Ayy
AYZ(Nx+1:2*Nx-1,1:Ny,1:Nz)=Ayz(Nx:-1:2,1:Ny,1:Nz);
clear Ayz
AZZ(Nx+1:2*Nx-1,1:Ny,1:Nz)=Azz(Nx:-1:2,1:Ny,1:Nz);
clear Azz

%========== Calculating FFT in X-direction 
AXX(:,1:Ny,1:Nz)=fft(AXX(:,1:Ny,1:Nz));
AXY(:,1:Ny,1:Nz)=fft(AXY(:,1:Ny,1:Nz));
AXZ(:,1:Ny,1:Nz)=fft(AXZ(:,1:Ny,1:Nz));
AYY(:,1:Ny,1:Nz)=fft(AYY(:,1:Ny,1:Nz));
AYZ(:,1:Ny,1:Nz)=fft(AYZ(:,1:Ny,1:Nz));
AZZ(:,1:Ny,1:Nz)=fft(AZZ(:,1:Ny,1:Nz));
%=========================================================================%



%======= Expanding dimension of each line in Y-direction to 2Ny-1 ========%
%=========================================================================%
AXX(:,Ny+1:2*Ny-1,1:Nz)=AXX(:,Ny:-1:2,1:Nz);
AXY(:,Ny+1:2*Ny-1,1:Nz)=-AXY(:,Ny:-1:2,1:Nz);
AXZ(:,Ny+1:2*Ny-1,1:Nz)=AXZ(:,Ny:-1:2,1:Nz);
AYY(:,Ny+1:2*Ny-1,1:Nz)=AYY(:,Ny:-1:2,1:Nz);
AYZ(:,Ny+1:2*Ny-1,1:Nz)=-AYZ(:,Ny:-1:2,1:Nz);
AZZ(:,Ny+1:2*Ny-1,1:Nz)=AZZ(:,Ny:-1:2,1:Nz);

%========== Calculating FFT in Y-direction 
% FFT in Matlab calculate DFT in x-direction, so for performing FFT ...
% ... in y-direction, first we permute x and y lines, peform the ...
% ... fft and again perform permutation to return back to first ...
% ... matrices dimension

AXX=fft(AXX(:,:,1:Nz),[],2); % Calculating FFT in y-direction
AXY=fft(AXY(:,:,1:Nz),[],2); % Calculating FFT in y-direction
AXZ=fft(AXZ(:,:,1:Nz),[],2); % Calculating FFT in y-direction
AYY=fft(AYY(:,:,1:Nz),[],2); % Calculating FFT in y-direction
AYZ=fft(AYZ(:,:,1:Nz),[],2); % Calculating FFT in y-direction
AZZ=fft(AZZ(:,:,1:Nz),[],2); % Calculating FFT in y-direction
%=========================================================================%


%======== Expanding dimension of each line in Z-direction to 2Nz-1 =======%
%=========================================================================%
AXX(:,:,Nz+1:2*Nz-1)=AXX(:,:,Nz:-1:2);
AXY(:,:,Nz+1:2*Nz-1)=AXY(:,:,Nz:-1:2);
AXZ(:,:,Nz+1:2*Nz-1)=-AXZ(:,:,Nz:-1:2);
AYY(:,:,Nz+1:2*Nz-1)=AYY(:,:,Nz:-1:2);
AYZ(:,:,Nz+1:2*Nz-1)=-AYZ(:,:,Nz:-1:2);
AZZ(:,:,Nz+1:2*Nz-1)=AZZ(:,:,Nz:-1:2);

%========== Calculating FFT in Z-direction 
FFT_AXX=fft(AXX,[],3); % Calculating FFT in z-direction
clear AXX
FFT_AXY=fft(AXY,[],3); % Calculating FFT in z-direction
clear AXY
FFT_AXZ=fft(AXZ,[],3); % Calculating FFT in z-direction
clear AXZ
FFT_AYY=fft(AYY,[],3); % Calculating FFT in z-direction
clear AYY
FFT_AYZ=fft(AYZ,[],3); % Calculating FFT in z-direction
clear AYZ
FFT_AZZ=fft(AZZ,[],3); % Calculating FFT in z-direction
clear AZZ
%=========================================================================%

end