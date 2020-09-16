%=========================================================================%           
                 % Iterative Method, Biconjugate gradient %
            % Applying Biconjugate gradient to obtainPx,Py,PZ %
%=========================================================================%
function [px,py,pz]=Biconjugate_Gradient(E_x,E_y,E_z,Nx,Ny,Nz,N,Inverse_Alpha,...
    INDEX_IN,E_vector,FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ,CT)


%================ Initial amounts of the PX, PY and PZ ===================%
px=0;
py=0;
pz=0;

Apkx=0;
Apky=0;
Apkz=0;

%============================ rk=E1-A*Pk;
rkx=(E_x-Apkx);
rky=(E_y-Apky);
rkz=(E_z-Apkz);

%============================= qk=rk;
qkx=rkx;
qky=rky;
qkz=rkz;
%============================= qk_bar=conj(qk);

Error=1;
counting=0;

%=================== Applying complex conjugate gradient =================%
while Error>CT  % CT is convergence threshold 
    
    %=== Performing inverse FFT to obtain matrix-vector multiplication ===%
    %=====================================================================%
    [Aqkx,Aqky,Aqkz]=Inverse_FFT(qkx,qky,qkz,Nx,Ny,Nz,INDEX_IN,Inverse_Alpha,...
              FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ);
    %=====================================================================%
    
    rkx_vector=reshape(rkx,[N,1]);
    rky_vector=reshape(rky,[N,1]);
    rkz_vector=reshape(rkz,[N,1]);
    
    Aqkx_vector=reshape(Aqkx,[N,1]);
    Aqky_vector=reshape(Aqky,[N,1]);
    Aqkz_vector=reshape(Aqkz,[N,1]);
    
    qkx_vector=reshape(qkx,[N,1]);
    qky_vector=reshape(qky,[N,1]);
    qkz_vector=reshape(qkz,[N,1]);
    
    
    %==================== alphak=(rk')*rk/(qk'*Aqk)=======================%
    %=====================================================================%
    alphak=(transpose(rkx_vector)*rkx_vector+transpose(rky_vector)*rky_vector...
        +transpose(rkz_vector)*rkz_vector)/(transpose(qkx_vector)*Aqkx_vector...
        +transpose(qky_vector)*Aqky_vector+transpose(qkz_vector)*Aqkz_vector);
    %=====================================================================%
    
               % Deleting unnessesary variables, matrices ...%
    clear Aqkx_vector Aqky_vector Aqkz_vector
    %=====================================================================%
    
    
    %========================= Pk=Pk+alphak*qk  ==========================%
    %=====================================================================%
    rk0x=rkx;
    clear rkx
    rk0y=rky;
    clear rky
    rk0z=rkz;
    clear rkz
    
    rk0x_vector=rkx_vector;
    clear rkx_vector
    rk0y_vector=rky_vector;
    clear rky_vector
    rk0z_vector=rkz_vector;
    clear rkz_vector
    
    rkx=rk0x-alphak*Aqkx;
    clear rk0x Aqkx
    rky=rk0y-alphak*Aqky;
    clear rk0y Aqky
    rkz=rk0z-alphak*Aqkz;
    clear rk0z Aqkz
    
    rkx_vector=reshape(rkx,[N,1]);
    rky_vector=reshape(rky,[N,1]);
    rkz_vector=reshape(rkz,[N,1]);
    
    px=(px+alphak*qkx).*INDEX_IN;
    py=(py+alphak*qky).*INDEX_IN;
    pz=(pz+alphak*qkz).*INDEX_IN;
    
    rk_vector=[rkx_vector;rky_vector;rkz_vector];
    Error=norm(rk_vector)/norm(E_vector);
    clear rk_vector
    %=====================================================================%
    
    
    %================== betak=rk_bar'*rk/(rk0_bar'*rk0)===================%
    %=====================================================================%
    betak=(transpose(rkx_vector)*rkx_vector+transpose(rky_vector)*rky_vector...
        +transpose(rkz_vector)*rkz_vector)/(transpose(rk0x_vector)*rk0x_vector...
        +transpose(rk0y_vector)*rk0y_vector+transpose(rk0z_vector)*rk0z_vector);
    clear rkx_vector rky_vector rkz_vector rk0x_vector rk0y_vector rk0z_vector
    %=====================================================================%
    
    
    %========================== qk=rk+betak*qk ===========================%
    %================== qk_bar=rk_bar+conj(betak)*qk_bar =================%
    %=====================================================================%
    qkx=(rkx+betak*qkx).*INDEX_IN;
    qky=(rky+betak*qky).*INDEX_IN;
    qkz=(rkz+betak*qkz).*INDEX_IN;
    
    counting=counting+1;
    %=====================================================================%
    
end
end