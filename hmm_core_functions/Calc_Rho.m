function [Rho] = Calc_Rho(transmat_t,R)
%Notation:
% Ns is the number of states.
% X = [x(1),x(2),...,x(T)] = hidden state sequence.
% Y = [y(1),y(2),...,y(T)] = observation sequence.
% t, represents the full time axis i.e. t=[1,2,3,4,.....].
% to(k) is the time of the k'th observation P(a|b), probability of a given b.
% transmat is the transition probability defined for a single bin time
% interval, i.e. transmat(i,j)=P(x(t+1)=j | x(t)=i).
% dt - time interval between observations.
% sum_l(fun(l)), is the sum of function (fun) over all possible values of l.
% Short description:
%   Calc_Rho calculates the sum of
% P(x(t)=i,x(t)=j, x(to(n+1))=m |x(to(n))=k) (Rho). over all t in the range to(n)<=t<to(n+1).
%INPUTS:
% transmat_t is an Ns x Ns x len_t array of transition probabilities
% defined for a variable time interval,
%   i.e. transmat_t(:,:,dt)=(transmat^dt). R, represents the set of time
%   intervals, i.e. The i'th time interval is given by, dt(i) =
%   sum_k(R{i}(k,1)^R{i}(k,2)).

% Ensure transmat_t is double precision
transmat_t = double(transmat_t);

len_t=size(transmat_t,3);
transmat=transmat_t(:,:,1);
Ns=size(transmat,1);
Rho{len_t}=zeros(Ns,Ns,Ns,Ns);
for t=1:len_t-1% Initialization.
    Rho{t}=zeros(Ns,Ns,Ns,Ns);
end
for i=1:Ns
    for j=1:Ns% if dt=1 than there is only one possible transition. 
        Rho{1}(i,j,i,j)=transmat(i,j);
    end
end
for t=2:len_t
    Rtemp=R{t};
    tempRsize=size(Rtemp,1);
    if tempRsize>1
        for Rind=1:tempRsize-1
            if Rind==1
                transmat_1 = double(reshape(transmat_t(:,:,Rtemp(Rind,1)),[Ns Ns]));
                if Rtemp(Rind,2) == 1
                    Rho_1 = Rho{Rtemp(Rind,1)};
                else
                    Rho_1 = Rho_power_fast(transmat_1, Rho{Rtemp(Rind,1)} ,Rtemp(Rind,2));
                    transmat_1 = transmat_1 ^ double(Rtemp(Rind,2));
                end
                transmat_2 = double(reshape(transmat_t(:,:,Rtemp(Rind+1,1)),[Ns Ns]));
                if Rtemp(Rind+1,2)==1
                    Rho_2=Rho{Rtemp(Rind+1,1)};
                else
                    Rho_2=Rho_power_fast(transmat_2, Rho{Rtemp(Rind+1,1)} ,Rtemp(Rind+1,2));
                    transmat_2 = transmat_2 ^ double(Rtemp(Rind+1,2));
                end
                Rho_temp=Rho_product_fast(transmat_1 ,transmat_2 , Rho_1 ,Rho_2);
            else
                Rho_1=Rho_temp;
                transmat_1=transmat_1*transmat_2;
                transmat_2 = double(reshape(transmat_t(:,:,Rtemp(Rind+1,1)),[Ns Ns]));
                if Rtemp(Rind+1,2)==1
                    Rho_2=Rho{Rtemp(Rind+1,1)};
                else
                    Rho_2=Rho_power_fast(transmat_2, Rho{Rtemp(Rind+1,1)} ,Rtemp(Rind+1,2));
                    transmat_2 = transmat_2 ^ double(Rtemp(Rind+1,2));
                end
                Rho_temp=Rho_product_fast(transmat_1 ,transmat_2 , Rho_1 ,Rho_2);
            end
        end
        Rho{t} = Rho_temp;
    else
        Rho_1 = Rho{Rtemp(1,1)};
        transmat_1 = double(reshape(transmat_t(:,:,Rtemp(1,1)),[Ns Ns]));
        Rho{t} = Rho_power_fast(transmat_1, Rho_1,Rtemp(1,2));
    end
end
end