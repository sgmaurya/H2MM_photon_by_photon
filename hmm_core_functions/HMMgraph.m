
%% Load files
f=figure;

load HMM_parm.mat;

% Look for the flag file first in the current directory, then in the parent.
if isfile('HMMtypeFlag.mat')
    load('HMMtypeFlag.mat')
elseif isfile('../HMMtypeFlag.mat')
    load('../HMMtypeFlag.mat')  
end
%% Set graph parameters and labels

numStates = length(Prior);
rateMat = K; %Rate Matrix

% Edit population distribution typ:
population1 = Prior; %Steady state distribution of states (aVNorm)
%populationType='(Steady-state distribution)';
population2 = Prio_sort; % Estimated prior from HMM
populationType='(Distribution: SS|Prior)';

FRETvalues = Obs; %FRET values of states

rateMat(rateMat<0)=0;
% rateMat=round(rateMat);

rateMat=round(rateMat/1000,2,'significant');% In ms^-1
G = digraph(rateMat);

nodeNames ={};

% for three color FRET

if exist('indx')&&indx==3
   
    for i=1:numStates/2
        nodeNames{i}= ['Bright ', num2str(round(FRETvalues(i),2)),newline,'(',num2str(100*round(population1(i),2)),'|'...
            num2str(100*round(population2(i),2)),'%)'];
        
        nodeNames{i+numStates/2}= ['Dark ', num2str(round(FRETvalues(i+numStates/2),2)),newline,'(',num2str(100*round(population1(i+numStates/2),2)),'|'...
            num2str(100*round(population2(i+numStates/2),2)),'%)'];
    end
    
else
    for i=1:numStates
        nodeNames{i}= [num2str(round(FRETvalues(i),2)),newline,'(',num2str(100*round(population1(i),2)),'|'...
            num2str(100*round(population2(i),2)),'%)'];
        
    end
end


G.Nodes.Name = nodeNames';
       
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight)+ 0.000001;
NodeSize = 150*population1';
NodeCol = colororder;

p = plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths,'Layout','force');
p.Marker = 'o';
p.NodeColor = NodeCol(mod(0:numStates-1,size(NodeCol,1))+1,:);
p.MarkerSize = NodeSize;
p.ArrowSize=20;
p.EdgeColor ='k';

p.NodeFontSize=30;
p.EdgeFontSize=30;

name=pwd;
count = length(name);
while name(count)~='\'
    count=count-1;
    if name(count)=='_'
        name(count)=' ';
    end
end
title({name(count+1:end),[populationType, ', rates in ms^{-1}']});
axis off
save('Graphs','f')

