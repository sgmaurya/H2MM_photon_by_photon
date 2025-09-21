function [Rho] = Rho_power_fast(A_dt0 ,Rho0 ,power)
N=size(A_dt0,1);
A0=eye(N,N);
Rho=Rho0;
for i=1:power-1
    A0=A0*A_dt0;
    [Rho] = Rho_product_fast(A0 ,A_dt0 , Rho, Rho0);
end
end
