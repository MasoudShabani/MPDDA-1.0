%=========================================================================%
% Finding index of the dipoles outside of the NPs and excluding inside NPs%
%=========================================================================%
function [Outside_Index]=Excluding_NPs(x_plane,y_plane,z_plane,N_plane,Np_shape,Lx,Ly,Lz,d_inter,Structure,arrangement)

if Structure=="monomeric" % Monomer structure
    X0=0; 
    Y0=0;
    Z0=0;
    
    if Np_shape=="spherical"       % NPs are sphere
        % Index of cubes inside NP
        Index_in=find((sqrt((x_plane-X0).^2+(y_plane-Y0).^2+(z_plane-Z0).^2)<=(Lx/2)));
        
    elseif Np_shape=="ellipsoid" %NPs are ellipsoid, head-tail orientation in z-direction
        % Index of cubes inside NP
        Index_in=find((sqrt((x_plane-X0).^2/((Lx/2)^2)+(y_plane-Y0).^2/((Ly/2)^2)+...
            (z_plane-Z0).^2/((Lz/2)^2))<=1));
        
        
    elseif Np_shape=="rod"    %NPs are Rod with caps, vertically oriented
        % Index of cubes inside NP
        Index_in=find((sqrt((x_plane-X0).^2+(y_plane-Y0).^2+(abs(z_plane-Z0)-Lz/2+Lx/2).^2)<=(Lx/2)...
            & abs(z_plane-Z0)>(Lz/2-Lx/2))|((sqrt((x_plane-X0).^2+(y_plane-Y0).^2)<=(Lx/2))...
            &abs(z_plane-Z0)<=(Lz/2-Lx/2)));
        
    elseif Np_shape=="rec_block" %NPs are Rectangular block, vertically oriented
        % Index of cubes inside NP
        Index_in=find((abs(x_plane-X0)<=Lx/2 & abs(y_plane-Y0)<=Ly/2 & abs(z_plane-Z0)<=Lz/2));
    end
    
    %=====================================================================%
elseif Structure=="dimeric" % Dimer, structure=2
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
    
    if Np_shape=="spherical"         % NPs are sphere
        % Index of cubes inside NP
        Index_in=find((sqrt((x_plane-X10).^2+(y_plane-Y10).^2+(z_plane-Z10).^2)<=(Lx/2))|...
            (sqrt((x_plane-X20).^2+(y_plane-Y20).^2+(z_plane-Z20).^2)<=(Lx/2)));
        
    elseif Np_shape=="ellipsoid" %NPs are ellipsoid, head-tail orientation in z-direction
        % Index of cubes inside NP
        Index_in=find((sqrt((x_plane-X10).^2/((Lx/2)^2)+(y_plane-Y10).^2/((Ly/2)^2)+...
            (z_plane-Z10).^2/((Lz/2)^2))<=1)|(sqrt((x_plane-X20).^2/((Lx/2)^2)+...
            (y_plane-Y20).^2/((Ly/2)^2)+(z_plane-Z20).^2/((Lz/2)^2))<=1));
        
        
    elseif Np_shape=="rod"    %NPs are Rod with caps, vertically oriented
        % Index of cubes inside NP
        Index_in=find((sqrt((x_plane-X10).^2+(y_plane-Y10).^2+(abs(z_plane-Z10)-Lz/2+Lx/2).^2)<=(Lx/2)...
            & abs(z_plane-Z10)>(Lz/2-Lx/2))|((sqrt((x_plane-X10).^2+(y_plane-Y10).^2)<=(Lx/2))...
            &abs(z_plane-Z10)<=(Lz/2-Lx/2))|...
            (sqrt((x_plane-X20).^2+(y_plane-Y20).^2+(abs(z_plane-Z20)-Lz/2+Lx/2).^2)<=(Lx/2)...
            & abs(z_plane-Z20)>(Lz/2-Lx/2))|((sqrt((x_plane-X20).^2+(y_plane-Y20).^2)<=(Lx/2))...
            &abs(z_plane-Z20)<=(Lz/2-Lx/2)));
        
    elseif Np_shape=="rec_block" %NPs are Rectangular block, vertically oriented
        % Index of cubes inside NP
        Index_in=find((abs(x_plane-X10)<=Lx/2 & abs(y_plane-Y10)<=Ly/2 & abs(z_plane-Z10)<=Lz/2)|...
            (abs(x_plane-X20)<=Lx/2 & abs(y_plane-Y20)<=Ly/2 & abs(z_plane-Z20)<=Lz/2));
    end
end

Multiply_Nps=zeros(N_plane,1);
Multiply_Nps(Index_in,1)=1;

Mult=ones(N_plane,1);
Outside_Index=Mult-Multiply_Nps;
clear Multiply_Nps Index_in Multiply_coeff

end