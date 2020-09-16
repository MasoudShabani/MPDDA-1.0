%=========================================================================%
%============ Calculating Incident electric filed components =============%
%=========================================================================%

function [E_x,E_y,E_z,E_vector]=Incident_Field(Lambda,IB,nb,r_block,kvec,K0,INDEX_INSIDE,Nx,Ny,Nz,E0,z0,Waist_r)

kr=kvec(1)*r_block(:,1)+kvec(2)*r_block(:,2)+kvec(3)*r_block(:,3);

%======== Here we are gonna chose the type of the incident light =========%
%=========================================================================%
if IB=="plane wave"            % Incident beam is plane wave
    expikr=exp(1i*kr);
    Ex=E0(1)*expikr.*INDEX_INSIDE;
    Ey=E0(2)*expikr.*INDEX_INSIDE;
    Ez=E0(3)*expikr.*INDEX_INSIDE;
    
    E_x=reshape(Ex,[Nx,Ny,Nz]);
    E_y=reshape(Ey,[Nx,Ny,Nz]);
    E_z=reshape(Ez,[Nx,Ny,Nz]);
    
    E_vector=[Ex;Ey;Ez];
    
elseif IB=="gaussian"        % Incident beam is Gaussian 
    
    if K0(1)==1
        r=(r_block(:,2).^2+r_block(:,3).^2).^(0.5);
        z=r_block(:,1);
    elseif K0(2)==1
        r=(r_block(:,1).^2+r_block(:,3).^2).^(0.5);
        z=r_block(:,2);
    else
        r=(r_block(:,1).^2+r_block(:,2).^2).^(0.5);
        z=r_block(:,3);
    end
    
    
    W0=Waist_r*Lambda;                     %Waist radius of the beam
    K=2*pi*nb/Lambda;                      %wave number
    zR=pi*(W0.^2).*nb./Lambda;
    Qz=atan((z-z0)./zR);                    % Is the Gouy phase at z
    Wz=W0.*((1+((z-z0)./zR).^2).^0.5);      % the spot size parameter
    Rz_inverse=z./((z-z0).^2+zR.^2);        % inverse of radius of curvature
    
    E=(W0./Wz).*exp(-(r./Wz).^2-1i*(K.*(z-z0)+K.*(r.^2).*Rz_inverse/2-Qz));
    Ex=E0(1)*E.*INDEX_INSIDE;
    Ey=E0(2)*E.*INDEX_INSIDE;
    Ez=E0(3)*E.*INDEX_INSIDE;
    
    E_x=reshape(Ex,[Nx,Ny,Nz]);
    E_y=reshape(Ey,[Nx,Ny,Nz]);
    E_z=reshape(Ez,[Nx,Ny,Nz]);
    
    E_vector=[Ex;Ey;Ez];
end

end