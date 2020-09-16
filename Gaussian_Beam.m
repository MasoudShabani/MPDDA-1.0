% Guassian Incident light

clear all 
clc 

E0=1;
nb=1.34;    % Refractive index of medium
Lamda=525;  % incident light wavelength 
rr=-20*Lamda:Lamda/10:20*Lamda;    % the radial distance from the center axis of the beam
zz=-100*Lamda:Lamda/10:100*Lamda;    % the axial distance from the beam's focus  
[r,z]=meshgrid(rr,zz);

W0=2*Lamda;                       %Waist radius of the beam
k=2*pi*nb/Lamda;                     %wave number
zR=pi*(W0.^2).*nb./Lamda;    
Qz=atan(z./zR);                    % Is the Gouy phase at z
Wz=W0.*((1+(z./zR).^2).^0.5);      % the spot size parameter   
Rz_inverse=z./(z.^2+zR.^2);         % inverse of radius of curvature 

Egaus=E0.*(W0./Wz).*exp(-(r./Wz).^2-1i*(k.*z+k.*(r.^2).*Rz_inverse/2-Qz));
Igaus=abs(Egaus).^2;

surf(r/Lamda,z/Lamda,Igaus,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
shading interp;
waitforbuttonpress
%axis equal
hold on
