function [Rho12] = Rho_product_fast(A_dt1 ,A_dt2 ,Rho1 ,Rho2)
N=size(Rho1,1);
Rho12=zeros(N,N,N,N);
for m = 1:N
    for n =1:N
        Rho12(:,:,m,n) = ...
            reshape(reshape(Rho1(:,:,m,n),[N N])*A_dt2+A_dt1*reshape(Rho2(:,:,m,n),[N N]),[1 1 N N]);
    end
end
end


% The old version is slower but simpler to understand. 
% function [Rho12] = Rho_product(A_dt1 ,A_dt2 ,Rho1 ,Rho2)
% N=size(Rho1,1);
% for i = 1:N
%     for m = 1:N
%         for n =1:N
%             for k = 1:N
%                 Rho12(i,m,n,k) = reshape(Rho1(i,m,n,:),[1 N])*reshape(A_dt2(:,k),[N 1])+...
%                     reshape(A_dt1(i,:),[1 N])*reshape(Rho2(:,m,n,k),[N 1]);
%             end
%         end
%     end
% end
% end