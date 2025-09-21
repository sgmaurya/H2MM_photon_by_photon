function M = CalculatePowerOfTransMatrices( R, TransMat )
% CalculatePowerOfTransMatrices is used to calculate the transition matrix
% over all time intervals represented in R.
% INPUTS:
% R

%|------------------------------------------------------------------------|
%| Calculating the matrices from R (R is the result of Factorization.m)   |
%|------------------------------------------------------------------------|
N=size(TransMat,1);
M=zeros(N,N,size(R,2));
M(:,:,1) = TransMat;
for ii = 2: length(R)
    Mtemp = eye(length(TransMat)); % Identity matrix
    for iii = 1 : size(R{1,ii} ,1)
        M(:,:,ii) = Mtemp *...
            M(:,:,R{1,ii}(iii, 1)) ^ double(R{1,ii}(iii, 2));
        Mtemp = M(:,:,ii);
    end
end

end