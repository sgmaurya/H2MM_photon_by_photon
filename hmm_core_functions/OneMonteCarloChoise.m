function [ Choise ] = OneMonteCarloChoise( ChoisProbVec )
%OneMonteCarloChoise uses ChoisProbVec a vector that describes probabilities and  
%uses a rand function to make a random choice outputed as an integer
len=length(ChoisProbVec);
s=sum(ChoisProbVec);
if abs(s-1)>10e-8
    disp('sum(ChoisProbVec)~=0')
    Choise=[];
    return
end
temp_size=size(ChoisProbVec);
if temp_size(1)>temp_size(2)
    ChoisProbVec=ChoisProbVec';
end
pic=cumsum(ChoisProbVec);
rand_num=rand(1);
ischoice=0;
for i=1:length(pic)
    if rand_num<=pic(i);
        Choise=i;
        ischoice=1;
        break
    end
end
if ischoice==0
    1
end

end

