% HMMgraph_function.m
% This is a function-based version of HMMgraph. It accepts all necessary
% data as input arguments, making it robust and independent of the workspace.

function HMMgraph_function(K, Prior, Prio_sort, Obs, indx, guessName, output_folder_path)

    f = figure('Visible', 'off'); % Create figure but keep it hidden until the end

    % --- The core logic is the same, but uses input arguments ---
    numStates = length(Prior);
    rateMat = K; % Rate Matrix from input

    population1 = Prior; % Steady state distribution from input
    population2 = Prio_sort; % Estimated prior from input
    populationType = '(Distribution: SS|Prior)';

    FRETvalues = Obs; % FRET values from input

    rateMat(rateMat<0) = 0;
    rateMat = round(rateMat/1000, 2, 'significant'); % In ms^-1
    G = digraph(rateMat);

    nodeNames = {};

    % Logic for three-color FRET
    if exist('indx','var') && ~isempty(indx) && indx == 3 && numStates > 1
        for i = 1:numStates/2
            nodeNames{i} = ['Bright ', num2str(round(FRETvalues(i),2)), newline, '(', num2str(100*round(population1(i),2)), '|', num2str(100*round(population2(i),2)), '%)'];
            nodeNames{i+numStates/2} = ['Dark ', num2str(round(FRETvalues(i+numStates/2),2)), newline, '(', num2str(100*round(population1(i+numStates/2),2)), '|', num2str(100*round(population2(i+numStates/2),2)), '%)'];
        end
    else % Default logic for 2-state or other models
        for i = 1:numStates
            nodeNames{i} = [num2str(round(FRETvalues(i),2)), newline, '(', num2str(100*round(population1(i),2)), '|', num2str(100*round(population2(i),2)), '%)'];
        end
    end

    G.Nodes.Name = nodeNames';
           
    LWidths = 5 * G.Edges.Weight / max(G.Edges.Weight) + 1e-6;
    % Ensure NodeSize is a column vector and handle single-state case
    NodeSize = 150 * population1(:) + 10; % Add small base size

    p = plot(G, 'EdgeLabel', G.Edges.Weight, 'LineWidth', LWidths, 'Layout', 'force');
    p.Marker = 'o';
    p.NodeColor = 'r';
    p.MarkerSize = NodeSize;
    p.ArrowSize = 20;
    p.EdgeColor = 'k';
    p.NodeFontSize = 24;
    p.EdgeFontSize = 24;

    % Create a clean title from the folder path
    [~, folder_name, ~] = fileparts(output_folder_path);
    folder_name = strrep(folder_name, '_', ' '); % Make it readable
    
    title({folder_name, [populationType, ', rates in ms^{-1}']});
    axis off;
    
    % Make the figure visible now that it's fully drawn
    set(f, 'Visible', 'on');

    % Save the figure as a PNG for easy viewing
    saveas(f, fullfile(output_folder_path, 'HMM_Kinetic_Graph.png'));
    
    % The original save can be kept if needed, but saveas is often better
    try
        save(fullfile(output_folder_path, 'Graphs.mat'), 'f');
    catch ME_save
        warning('Could not save figure handle to .mat file. PNG was saved. Error: %s', ME_save.message);
    end

end