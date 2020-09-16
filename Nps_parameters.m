%=========================================================================%
%===== Finding dimension of the Rectangular block that encompass NP ======%
%=========================================================================%
function [Lx,Ly,Lz]= Nps_parameters(r_eff,Np_shape)

if Np_shape=="spherical"      
    d_eff=2*r_eff;     % Effective diameter
    Lx=d_eff;
    Ly=d_eff;
    Lz=d_eff;
    
elseif Np_shape=="ellipsoid" 
    ARyx=input('\n\nEnter the ratio of y-semi axis to x-semi axis:');
    ARzx=input('\nEnter the ratio of z-semi axis to x-semi axis:');
    clc
    a=((1/(ARyx*ARzx))^(1/3))*r_eff;    % Semi- minor axis in x- direction
    b=ARyx*a;                           % Semi- minor axis in y- direction
    c=ARzx*a;                           % Semi- major axis in z- direction
    
    Lx=2*a;                             % Length of Np in x-direction
    Ly=2*b;                             % Length of Np in y-direction
    Lz=2*c;                             % Length of Np in z-direction
    
elseif Np_shape=="rod"
    volume=4*pi/3*(r_eff^3);
    AR=input('\nEnter the aspect ratio, ratio of the long axis to short one:');
    clc
    r=(volume/(pi*(2*(AR-1)+4/3)))^(1/3); % Raduis of the Rod
    diameter=2*r;                         % Diameter of the Rod
    higth=AR*diameter;                    % Higth of the Rod
    
    Lx=diameter;                       % Length of Np in x-direction
    Ly=diameter;                       % Length of Np in y-direction
    Lz=higth;                          % Length of Np in z-direction
    
elseif Np_shape=="rec_block"              
    volume=4*pi/3*(r_eff^3);
    ARyx=input('\nEnter the ratio of the y-axis to x one:');
    ARzx=input('\nEnter the ratio of the z-axis to x one:');
    clc
    a=(volume/(ARyx*ARzx))^(1/3);    % First side of rectangular block
    b=ARyx*a;                        % Second side of rectangular block
    c=ARzx*a;                        % Third side of rectangular block
    
    Lx=a;                            % Length of Np in x-direction
    Ly=b;                            % Length of Np in y-direction
    Lz=c;                            % Length of Np in z-direction
    
end
end