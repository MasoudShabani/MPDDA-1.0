
%=========================================================================%
%================ Finding total E_field at each nanocube =================%
%=========================================================================%
function [Ex_out,Ey_out,Ez_out]=E_total(GPU,N_plane,Nx,Ny,Nz,INDEX_IN,Inverse_Alpha,Nx_target,...
    Ny_target,Nz_target,x_plane,y_plane,z_plane,r_block,Np_shape,plane_2D,...
    Lx,Ly,Lz,kvec,d_inter,Structure,Ex_incident,Ey_incident,Ez_incident,px,py,pz,arrangement)

%========== Defining new extended Rectangular box ===============%
if GPU==1
    qkx=zeros(Nx,Ny,Nz,'gpuArray');
    qky=zeros(Nx,Ny,Nz,'gpuArray');
    qkz=zeros(Nx,Ny,Nz,'gpuArray');
    INDEX_IN_extended=zeros(Nx,Ny,Nz,'gpuArray');
    Inverse_Alpha_extended=zeros(Nx,Ny,Nz,'gpuArray');
elseif GPU==0
    qkx=zeros(Nx,Ny,Nz);
    qky=zeros(Nx,Ny,Nz);
    qkz=zeros(Nx,Ny,Nz);
    INDEX_IN_extended=zeros(Nx,Ny,Nz);
    Inverse_Alpha_extended=zeros(Nx,Ny,Nz);
end

if plane_2D=="xy"
    Inverse_Alpha_extended((Nx-Nx_target)/2+1:(Nx+Nx_target)/2,(Ny-Ny_target)/2+1:(Ny+Ny_target)/2,1:Nz_target)=Inverse_Alpha;
    Inverse_Alpha=[];
    Inverse_Alpha=Inverse_Alpha_extended;
    
    INDEX_IN_extended((Nx-Nx_target)/2+1:(Nx+Nx_target)/2,(Ny-Ny_target)/2+1:(Ny+Ny_target)/2,1:Nz_target)=INDEX_IN;
    
    INDEX_IN=[];
    INDEX_IN=1;
    qkx((Nx-Nx_target)/2+1:(Nx+Nx_target)/2,(Ny-Ny_target)/2+1:(Ny+Ny_target)/2,1:Nz_target)=px;
    qky((Nx-Nx_target)/2+1:(Nx+Nx_target)/2,(Ny-Ny_target)/2+1:(Ny+Ny_target)/2,1:Nz_target)=py;
    qkz((Nx-Nx_target)/2+1:(Nx+Nx_target)/2,(Ny-Ny_target)/2+1:(Ny+Ny_target)/2,1:Nz_target)=pz;

elseif plane_2D=="xz"
    Inverse_Alpha_extended((Nx-Nx_target)/2+1:(Nx+Nx_target)/2,1:Ny_target,(Nz-Nz_target)/2+1:(Nz+Nz_target)/2)=Inverse_Alpha;
    Inverse_Alpha=[];
    Inverse_Alpha=Inverse_Alpha_extended;
    
    INDEX_IN_extended((Nx-Nx_target)/2+1:(Nx+Nx_target)/2,1:Ny_target,(Nz-Nz_target)/2+1:(Nz+Nz_target)/2)=INDEX_IN;
    
    INDEX_IN=[];
    INDEX_IN=1;
    qkx((Nx-Nx_target)/2+1:(Nx+Nx_target)/2,1:Ny_target,(Nz-Nz_target)/2+1:(Nz+Nz_target)/2)=px;
    qky((Nx-Nx_target)/2+1:(Nx+Nx_target)/2,1:Ny_target,(Nz-Nz_target)/2+1:(Nz+Nz_target)/2)=py;
    qkz((Nx-Nx_target)/2+1:(Nx+Nx_target)/2,1:Ny_target,(Nz-Nz_target)/2+1:(Nz+Nz_target)/2)=pz;

elseif plane_2D=="yz"
    Inverse_Alpha_extended(1:Nx_target,(Ny-Ny_target)/2+1:(Ny+Ny_target)/2,(Nz-Nz_target)/2+1:(Nz+Nz_target)/2)=Inverse_Alpha;
    Inverse_Alpha=[];
    Inverse_Alpha=Inverse_Alpha_extended;
    
    INDEX_IN_extended(1:Nx_target,(Ny-Ny_target)/2+1:(Ny+Ny_target)/2,(Nz-Nz_target)/2+1:(Nz+Nz_target)/2)=INDEX_IN;
    
    INDEX_IN=[];
    INDEX_IN=1;
    qkx(1:Nx_target,(Ny-Ny_target)/2+1:(Ny+Ny_target)/2,(Nz-Nz_target)/2+1:(Nz+Nz_target)/2)=px;
    qky(1:Nx_target,(Ny-Ny_target)/2+1:(Ny+Ny_target)/2,(Nz-Nz_target)/2+1:(Nz+Nz_target)/2)=py;
    qkz(1:Nx_target,(Ny-Ny_target)/2+1:(Ny+Ny_target)/2,(Nz-Nz_target)/2+1:(Nz+Nz_target)/2)=pz;
    
end

[Outside_Index]=Excluding_NPs(x_plane,y_plane,z_plane,N_plane,Np_shape,Lx,Ly,Lz,d_inter,Structure,arrangement);


[rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,rjkrjk31_I,...
    rjkrjk32_I,rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,RJK]=RijRij(r_block);


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


%=== Performing inverse FFT to obtain matrix-vector multiplication ===%
%=====================================================================%
[Aqkx,Aqky,Aqkz]=Inverse_FFT(qkx,qky,qkz,Nx,Ny,Nz,INDEX_IN,Inverse_Alpha,...
    FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ);

if plane_2D=="xy"
    Apkx=Aqkx(1:Nx,1:Ny,(Nz+1)/2);
    Apky=Aqky(1:Nx,1:Ny,(Nz+1)/2);
    Apkz=Aqkz(1:Nx,1:Ny,(Nz+1)/2);
elseif plane_2D=="xz"
    Apkx=Aqkx(1:Nx,(Ny+1)/2,1:Nz);
    Apky=Aqky(1:Nx,(Ny+1)/2,1:Nz);
    Apkz=Aqkz(1:Nx,(Ny+1)/2,1:Nz);
elseif plane_2D=="yz"
    Apkx=Aqkx((Nx+1)/2,1:Ny,1:Nz);
    Apky=Aqky((Nx+1)/2,1:Ny,1:Nz);
    Apkz=Aqkz((Nx+1)/2,1:Ny,1:Nz);
end

Apx=reshape(Apkx,[N_plane,1]);
Apy=reshape(Apky,[N_plane,1]);
Apz=reshape(Apkz,[N_plane,1]);

Ex_scat=-(Apx);
Ey_scat=-(Apy);
Ez_scat=-(Apz);

clear Axx Axy Axz Ayy Ayz Azz rhat1_I rhat2_I rhat3_I rhat4_I rhat5_I rhat6_I ...
    rhat31_I rhat32_I rhat33_I rhat34_I rhat35_I rhat36_I

%=========================================================================%
%====== excluding the contribuation of the dipoles inside the NPs ========%
%=========================================================================%
Ex_out=(Ex_incident+Ex_scat).*Outside_Index;
Ey_out=(Ey_incident+Ey_scat).*Outside_Index;
Ez_out=(Ez_incident+Ez_scat).*Outside_Index;

%=========================================================================%
end


