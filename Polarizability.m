%=========================================================================%
%%Obtaining inverse of polarizibilty of each nanocube at different Lambda%%
%=========================================================================%

function [Inverse_Alpha]=Polarizability(GPU,kvec,eps_NP_eb,INDEX_IN,d,E0)
k0 = 2*pi;
b1 = -1.891531;
b2 = 0.1648469;
b3 = -1.7700004;
dcube = d^3;
a_hat = kvec/norm(kvec);
e_hat = E0/norm(E0);
S = 0;

for j = 1:3
    S = S + (a_hat(j)*e_hat(j))^2;
end

if GPU==1
    S=gpuArray(S);
end


a_CM =3*dcube/(4*pi)*(eps_NP_eb - 1)./(eps_NP_eb + 2); % Clausius-Mossotti
anr=a_CM./(1 + (a_CM/dcube).*(b1+eps_NP_eb*b2+eps_NP_eb*b3*S)*((norm(kvec)*d)^2));

if GPU==1     % running in GPU
    aLDR=gpuArray(anr./(1-2/3*1i*(anr/dcube)*((norm(kvec)*d)^3)));
elseif GPU==0 % running in CPU
    aLDR=(anr./(1-2/3*1i*(anr/dcube)*((norm(kvec)*d)^3)));
end

Inverse_Alpha=1/aLDR*INDEX_IN;
end
