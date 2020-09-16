%=========================================================================%
%===== Performing inverse FFT to obtain matrix-vector multiplication =====%
%=========================================================================%

function [Aqkx,Aqky,Aqkz]=Inverse_FFT(qkx,qky,qkz,Nx,Ny,Nz,INDEX_IN,Inverse_Alpha,...
              FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ)
          
%==== Aqk=A*qk ====%

%==== Calculate FFT in y-direction for px (qx), py (qy), and pz (qz)======%
%=========================================================================%
FFT_qkx_Y=fft(qkx,2*Ny-1,2); % Calculating FFT in y direction
FFT_qky_Y=fft(qky,2*Ny-1,2); % Calculating FFT in y direction
FFT_qkz_Y=fft(qkz,2*Ny-1,2); % Calculating FFT in y direction


%===== Calculate FFT in x-direction of the vectors in previous step ======%
%=========================================================================%
FFT_qkx_x=fft(FFT_qkx_Y,2*Nx-1,1); % Calculating FFT in x direction
clear FFT_qkx_Y
FFT_qky_x=fft(FFT_qky_Y,2*Nx-1,1); % Calculating FFT in x direction
clear FFT_qky_Y
FFT_qkz_x=fft(FFT_qkz_Y,2*Nx-1,1); % Calculating FFT in x direction
clear FFT_qkz_Y


%===== Calculate FFT in z-direction of the vectors in previous step ======%
%=========================================================================%
FFT_qkx_z=fft(FFT_qkx_x,2*Nz-1,3);% Calculating FFT in Z direction
clear FFT_qkx_x
FFT_qky_z=fft(FFT_qky_x,2*Nz-1,3);% Calculating FFT in Z direction
clear FFT_qky_x
FFT_qkz_z=fft(FFT_qkz_x,2*Nz-1,3);% Calculating FFT in Z direction
clear FFT_qkz_x


%==== Element-wise multiplicaition of FFT of matrix and FFT of vector ====%
%=========================================================================%
FFT_APX=FFT_AXX.*FFT_qkx_z+FFT_AXY.*FFT_qky_z+FFT_AXZ.*FFT_qkz_z;
FFT_APY=FFT_AXY.*FFT_qkx_z+FFT_AYY.*FFT_qky_z+FFT_AYZ.*FFT_qkz_z;
FFT_APZ=FFT_AXZ.*FFT_qkx_z+FFT_AYZ.*FFT_qky_z+FFT_AZZ.*FFT_qkz_z;
clear FFT_qkx_z FFT_qky_z FFT_qkz_z


%==================== Performing ifft in z-direction =====================%
%=========================================================================%
IFFT_APX_Z=ifft(FFT_APX,[],3);% Calculating IFFT in Z direction
clear FFT_APX
IFFT_APY_Z=ifft(FFT_APY,[],3);% Calculating IFFT in Z direction
clear FFT_APY
IFFT_APZ_Z=ifft(FFT_APZ,[],3);% Calculating IFFT in Z direction
clear FFT_APZ


%=================== Performing ifft in x-direction ======================%
%=========================================================================%
IFFT_APX_X=ifft(IFFT_APX_Z(1:2*Nx-1,1:2*Ny-1,1:Nz));
clear IFFT_APX_Z
IFFT_APY_X=ifft(IFFT_APY_Z(1:2*Nx-1,1:2*Ny-1,1:Nz));
clear IFFT_APY_Z
IFFT_APZ_X=ifft(IFFT_APZ_Z(1:2*Nx-1,1:2*Ny-1,1:Nz));
clear IFFT_APZ_Z


%==================== Performing ifft in y-direction =====================%
%=========================================================================%
IFFT_APX=ifft(IFFT_APX_X(1:Nx,1:2*Ny-1,1:Nz),[],2);
clear IFFT_APX_X
IFFT_APY=ifft(IFFT_APY_X(1:Nx,1:2*Ny-1,1:Nz),[],2);
clear IFFT_APY_X
IFFT_APZ=ifft(IFFT_APZ_X(1:Nx,1:2*Ny-1,1:Nz),[],2);
clear IFFT_APZ_X


%==== Negating contribution of the dipoles outside of the NPs boundary====%
%=========================================================================%

Aqkx=IFFT_APX(1:Nx,1:Ny,1:Nz).*INDEX_IN+Inverse_Alpha.*qkx.*INDEX_IN;
clear IFFT_APX
Aqky=IFFT_APY(1:Nx,1:Ny,1:Nz).*INDEX_IN+Inverse_Alpha.*qky.*INDEX_IN;
clear IFFT_APY
Aqkz=IFFT_APZ(1:Nx,1:Ny,1:Nz).*INDEX_IN+Inverse_Alpha.*qkz.*INDEX_IN;
clear IFFT_APZ

end
