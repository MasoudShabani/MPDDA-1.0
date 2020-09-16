%=========================================================================%
%===== Finding index of dipoles inside NPs & ingoring other elements =====%
%=========================================================================%

function [INDEX_INSIDE]=INDEX_INSIDE_NP(GPU,X,Y,Z,N,Np_shape,Lx,Ly,Lz,Structure,d_inter,arrangement)

%========================== Monomeric structure ==========================%
%=========================================================================%
if Structure=="monomeric"
    X0=0;
    Y0=0;
    Z0=0;
    
    if Np_shape=="spherical"    % NPs are sphere
        % Index of cubes inside NP
        Index_in=find((sqrt((X-X0).^2+(Y-Y0).^2+(Z-Z0).^2)<=(Lx/2)));
        
    elseif Np_shape=="ellipsoid" % NPs are ellipsoid, head-tail orientation in z-direction
        % Index of cubes inside NP
        Index_in=find((sqrt((X-X0).^2/((Lx/2)^2)+(Y-Y0).^2/((Ly/2)^2)+...
            (Z-Z0).^2/((Lz/2)^2))<=1));
        
        
    elseif Np_shape=="rod" % NPs are Rod with caps, vertically oriented
        % Index of cubes inside NP
        Index_in=find((sqrt((X-X0).^2+(Y-Y0).^2+(abs(Z-Z0)-Lz/2+Lx/2).^2)<=(Lx/2)...
            & abs(Z-Z0)>(Lz/2-Lx/2))|((sqrt((X-X0).^2+(Y-Y0).^2)<=(Lx/2))...
            &abs(Z-Z0)<=(Lz/2-Lx/2)));
        
    elseif Np_shape=="rec_block" % NPs are Rectangular block, vertically oriented
        % Index of cubes inside NP
        Index_in=find((abs(X-X0)<=Lx/2 & abs(Y-Y0)<=Ly/2 & abs(Z-Z0)<=Lz/2));
    end
    %=========================================================================%
    
elseif Structure=="dimeric"
    if arrangement=="x_orient"
        X10=-(d_inter/2+Lx/2);
        Y10=0;
        Z10=0;
        X20=(d_inter/2+Lx/2);
        Y20=0;
        Z20=0;
    elseif arrangement=="y_orient"
        X10=0;
        Y10=-(d_inter/2+Ly/2);
        Z10=0;
        X20=0;
        Y20=(d_inter/2+Ly/2);
        Z20=0;
    elseif arrangement=="z_orient"
        X10=0;
        Y10=0;
        Z10=-(d_inter/2+Lz/2);
        X20=0;
        Y20=0;
        Z20=(d_inter/2+Lz/2);
    end
    
    if Np_shape=="spherical"       % NPs are sphere
        % Index of cubes inside NP
        Index_in=find((sqrt((X-X10).^2+(Y-Y10).^2+(Z-Z10).^2)<=(Lx/2))|...
            (sqrt((X-X20).^2+(Y-Y20).^2+(Z-Z20).^2)<=(Lx/2)));
        
    elseif Np_shape=="ellipsoid"   %NPs are ellipsoid and have head to tail orientation in z-direction
        % Index of cubes inside NP
        Index_in=find((sqrt((X-X10).^2/((Lx/2)^2)+(Y-Y10).^2/((Ly/2)^2)+...
            (Z-Z10).^2/((Lz/2)^2))<=1)|(sqrt((X-X20).^2/((Lx/2)^2)+...
            (Y-Y20).^2/((Ly/2)^2)+(Z-Z20).^2/((Lz/2)^2))<=1));
        
    elseif Np_shape=="rod"   %NPs are Rod with caps, vertically oriented
        % Index of cubes inside NP
        Index_in=find((sqrt((X-X10).^2+(Y-Y10).^2+(abs(Z-Z10)-Lz/2+Lx/2).^2)...
            <=(Lx/2)& abs(Z-Z10)>(Lz/2-Lx/2))|...
            ((sqrt((X-X10).^2+(Y-Y10).^2)<=(Lx/2))&abs(Z-Z10)<=(Lz/2-Lx/2))|...
            (sqrt((X-X20).^2+(Y-Y20).^2+(abs(Z-Z20)-Lz/2+Lx/2).^2)<=(Lx/2)...
            & abs(Z-Z20)>(Lz/2-Lx/2))|((sqrt((X-X20).^2+(Y-Y20).^2)<=(Lx/2))...
            &abs(Z-Z20)<=(Lz/2-Lx/2)));
        
    elseif Np_shape=="rec_block"     %NPs are Rectangular block, vertically oriented
        % Index of cubes inside NP
        Index_in=find((abs(X-X10)<=Lx/2 & abs(Y-Y10)<=Ly/2 & abs(Z-Z10)<=Lz/2)|...
            (abs(X-X20)<=Lx/2 & abs(Y-Y20)<=Ly/2 & abs(Z-Z20)<=Lz/2));
        
    end
end
%=========================================================================%
if GPU==1
    INDEX_INSIDE=zeros(N,1,'gpuArray'); % Defining a zero array in GPU
elseif GPU==0
    INDEX_INSIDE=zeros(N,1);            % Defining a zero array in CPU
end
% Considering dipoles inside the Np and ignoring contribution of the other terms
INDEX_INSIDE(Index_in,1)=1;
end