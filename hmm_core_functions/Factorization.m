function [ R, Augmented_T ] = Factorization( T, TAU )
% Factorization computes an efficient representation of the sparse time set,T.
% i.e. each time point, T(i) is given as the sum of two numbers,
% if T(1)=TAU, T(i) =R{i}(1,1)*R{i}(1,1)+R{i}(2,1)*R{i}(2,1);
% if T(1)>TAU, T(i) =R{i+1}(1,1)*R{i+1}(1,1)+R{i}(2,1)*R{i+1}(2,1);
%|------------------------------------------------------------------------|
%| INPUT                                                                  |
%| T                  : array of delta_t                                  |
%| TAU                : basic time interval                               |
%|------------------------------------------------------------------------|


%|------------------------------------------------------------------------|
%| OUTPUT                                                                 |
%| R is a cell array of matrices that holds the factorization elements.   |
%| Each matrix in R has 2 columns:
%| The first column is
%|------------------------------------------------------------------------|
R = Main_rutin( T, TAU );

%%%%%%%%%%%%%%%%%%%%%%  improving the R partition so for all i R{i}(k,2)=1
B = [];
temp=0;
Augmented_T=T;
for i=1:length(R)
    temp_ind=[];
    temp_ind=find(R{i}(:,2)>1);%
    if and(~isempty(temp_ind),size(R{i},1)>1)
        for j=1:length(temp_ind)
            temp=temp+1;
            temp_arrival_time=T(R{i}(temp_ind(j),1));
            B(temp,1:2)=[temp_arrival_time R{i}(temp_ind(j),2)];
        end
    end
end

if ~isempty(B)
    [B]= unique(B,'rows');
    T2=B(:,1).*B(:,2);
    T2=[T2' T];
    [Augmented_T indT2]=sort(T2);
    Augmented_T = unique(Augmented_T);
    R = Main_rutin( Augmented_T, TAU );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end
function [R] = Main_rutin(T, TAU)
%|------------------------------------------------------------------------|
%| INPUT                                                                  |
%| T                  : array of delta_t                                  |
%| TAU                : basic time interval                               |
%|------------------------------------------------------------------------|


%|------------------------------------------------------------------------|
%| OUTPUT                                                                 |
%| R is a cell array of matrices that holds the factorization elements.   |
%| Each matrix in R has 2 columns:
%| The first column is
%|------------------------------------------------------------------------|

Ts = unique(T);   % Same data as in T, but with no repetitions and sorted

% The first element of Ts should be TAU (J[1]=1)
if Ts(1) ~= TAU
    Ts = [TAU Ts];
end
%|------------------------------------------------------------------------|
%| Create an array of j = (delta_t)/(TAU)                                 |
%|------------------------------------------------------------------------|

J = int32(Ts / TAU);

%|------------------------------------------------------------------------|
%| First element                                          |
%|------------------------------------------------------------------------|

R{1, 1} = [1 1];

%|------------------------------------------------------------------------|
%| Calculating the i’th element                                           |
%|------------------------------------------------------------------------|

for i = 2 : length(Ts)
    
    d = idivide(J(i), J(i - 1));
    m = mod(J(i), J(i - 1));   
    R{1, i} = [(i-1) d];
    
    while m ~= 0        
        % Check if we already have the result in cell array M
        
        if any(m == J(1:i-1))
            ind=find(J(1:i-1)==m,1,'first');
            R{1, i} = [R{1, i}; ind 1];
            break           
        else
            k  = find(J < m, 1, 'last');
            d = idivide(m, J(k));
            m = mod(m ,J(k));
            
            R{1, i} = [R{1, i}; k d];
        end        
    end    
end

end