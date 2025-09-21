% --- START OF FILE HMMgraph_local.m ---
function HMMgraph_local()

%% --- Configuration for Graph Appearance ---
graph_config.node_size_method = 'linear_pop_old_style';
graph_config.node_linear_scale_factor = 200; 
graph_config.node_min_size_clamp = 12;
graph_config.node_max_size_clamp = 70;

graph_config.edge_min_width = 2.5;
graph_config.edge_max_width = 8;
graph_config.edge_width_scale_factor = 7;
graph_config.arrow_size = 18;
graph_config.arrow_position = 0.55; 

graph_config.show_edge_labels = true;
graph_config.edge_label_font_size = 11; 
graph_config.edge_label_font_weight = 'bold'; 
graph_config.edge_label_color = [0 0 0];    
graph_config.min_rate_for_label_ms_inv = 1e-5;

graph_config.edge_color_mode = 'fixed';
graph_config.fixed_edge_color = [0 0 0]; 

graph_config.node_text_font_size = 12; 
graph_config.node_text_font_weight = 'bold';
graph_config.node_text_color = [0 0 0]; 

graph_config.node_color_scheme = 'custom_for_2state_example';

graph_config.layout_style = 'force'; 
graph_config.force_layout_iterations = 300;
graph_config.force_layout_weight_effect = 'inverse'; 

%% Load files
f=figure('Name', 'HMM State Graph (Rates, Pop%, Black Arrows)', 'Visible', 'on', 'Color', 'w', 'Position', [200 200 500 450]);

param_file = 'HMM_parm.mat';
if ~isfile(param_file), error('HMMgraph: %s not found in current directory (%s).', param_file, pwd); end
loaded_params = load(param_file);

essential_vars = {'Prior_HMMp', 'Prio_sort_HMMp', 'Obs_fret_values_HMMp', 'K_HMMp'};
missing_vars_check = cellfun(@(v) ~isfield(loaded_params, v), essential_vars);
if any(missing_vars_check)
    fprintf('HMMgraph: Variables in HMM_parm.mat:\n');
    disp(fieldnames(loaded_params));
    error('HMMgraph: Essential HMM parameters missing: %s', strjoin(essential_vars(missing_vars_check), ', '));
end

population1 = loaded_params.Prior_HMMp; 
population2 = loaded_params.Prio_sort_HMMp; 
FRETvalues_nodes = loaded_params.Obs_fret_values_HMMp;
rateMat_s_inv = loaded_params.K_HMMp; 

fprintf('HMMgraph DEBUG: Loaded K_HMMp (rates in s^-1):\n'); disp(rateMat_s_inv);

HMMtypeFlag_found = false;
loaded_flag_data = struct();
if isfile('HMMtypeFlag.mat')
    loaded_flag_data = load('HMMtypeFlag.mat'); HMMtypeFlag_found = true;
elseif isfile(fullfile('..','HMMtypeFlag.mat'))
    loaded_flag_data = load(fullfile('..','HMMtypeFlag.mat')); HMMtypeFlag_found = true;
end
indx_for_graph = 1; 
if HMMtypeFlag_found
    if isfield(loaded_flag_data, 'indx_guessmodel_for_graph'), indx_for_graph = loaded_flag_data.indx_guessmodel_for_graph;
    elseif isfield(loaded_flag_data, 'indx'), indx_for_graph = loaded_flag_data.indx; end
end

%% Validate inputs and Basic Setup
numStates = length(population1);
if numStates == 0, warning('HMMgraph: numStates is 0.'); title('No States'); axis off; return; end
if ~isrow(population1), population1 = population1(:)'; end
if ~isrow(population2), population2 = population2(:)'; end 
if size(FRETvalues_nodes,1) > 1, FRETvalues_nodes = FRETvalues_nodes'; end
if any([length(FRETvalues_nodes), length(population2), size(rateMat_s_inv,1), size(rateMat_s_inv,2)] ~= numStates)
    error('HMMgraph: Dimension mismatch for %d states.', numStates);
end

rateMat_for_labels_ms_inv = rateMat_s_inv * 0.001; 
fprintf('HMMgraph DEBUG: Rates converted to ms^-1 (rateMat_for_labels_ms_inv) BEFORE thresholding:\n'); disp(rateMat_for_labels_ms_inv);

rateMat_for_labels_ms_inv(abs(rateMat_for_labels_ms_inv) < graph_config.min_rate_for_label_ms_inv) = 0;
fprintf('HMMgraph DEBUG: Rates in ms^-1 (rateMat_for_labels_ms_inv) AFTER thresholding (min_rate_for_label_ms_inv = %.0e ms^-1):\n', graph_config.min_rate_for_label_ms_inv); disp(rateMat_for_labels_ms_inv);

formatted_edge_labels = cell(numStates, numStates);
fprintf('HMMgraph DEBUG: Formatting edge labels...\n');
for i_r = 1:numStates
    for j_c = 1:numStates
        val_ms_current = rateMat_for_labels_ms_inv(i_r, j_c);
        formatted_edge_labels{i_r, j_c} = ''; 
        if i_r == j_c || val_ms_current == 0
            % No label
        else
            abs_val = abs(val_ms_current);
            if abs(abs_val - round(abs_val)) < 0.01 && round(abs_val) >= 1
                formatted_edge_labels{i_r, j_c} = sprintf('%d', round(val_ms_current));
            elseif abs_val < 0.001 
                formatted_edge_labels{i_r, j_c} = sprintf('%.1e', val_ms_current);
            elseif abs_val < 0.01
                 formatted_edge_labels{i_r, j_c} = sprintf('%.3f', val_ms_current);
            elseif abs_val < 0.1
                 formatted_edge_labels{i_r, j_c} = sprintf('%.2f', val_ms_current);
            else
                 formatted_edge_labels{i_r, j_c} = sprintf('%.1f', val_ms_current);
            end
        end
        if i_r ~= j_c
            fprintf('HMMgraph DEBUG: Label for transition %d -> %d (rate_ms: %.2e): "%s"\n', i_r, j_c, val_ms_current, formatted_edge_labels{i_r, j_c});
        end
    end
end

rateMat_for_digraph_structure = rateMat_s_inv;
if numStates > 0, rateMat_for_digraph_structure(1:numStates+1:end) = 0; end
fprintf('HMMgraph DEBUG: rateMat_for_digraph_structure (input to digraph()):\n'); disp(rateMat_for_digraph_structure);
G = digraph(rateMat_for_digraph_structure);

fprintf('HMMgraph DEBUG: After G = digraph(...):\n');
fprintf('HMMgraph DEBUG: isempty(G.Edges) = %d\n', isempty(G.Edges));
if ~isempty(G.Edges)
    fprintf('HMMgraph DEBUG: G.Edges table:\n'); disp(G.Edges);
    fprintf('HMMgraph DEBUG: Does G.Edges have a "Weight" column initially? %d\n', ismember('Weight', G.Edges.Properties.VariableNames));
    if ismember('Weight', G.Edges.Properties.VariableNames) && ~isempty(G.Edges.Weight)
        fprintf('HMMgraph DEBUG: Initial G.Edges.Weight values:\n'); disp(G.Edges.Weight);
    elseif ismember('Weight', G.Edges.Properties.VariableNames) && isempty(G.Edges.Weight)
        fprintf('HMMgraph DEBUG: Initial G.Edges.Weight column exists but is EMPTY.\n');
    end
end

if ~isempty(G.Edges)
    fprintf('HMMgraph DEBUG: (Re)Populating G.Edges.Weight with values from rateMat_s_inv...\n');
    new_weights = zeros(height(G.Edges),1);
    for k_edge = 1:height(G.Edges)
        s_node = G.Edges.EndNodes(k_edge,1); t_node = G.Edges.EndNodes(k_edge,2);
        new_weights(k_edge) = rateMat_s_inv(s_node, t_node); 
    end
    G.Edges.Weight = new_weights; 
    fprintf('HMMgraph DEBUG: G.Edges.Weight (re)populated. Values:\n'); disp(G.Edges.Weight);
else
    fprintf('HMMgraph DEBUG: SKIPPED (re)populating G.Edges.Weight because G.Edges is empty.\n');
end

%% Node Names
nodeNames = cell(numStates,1);
if indx_for_graph == 3 && numStates > 1 && mod(numStates,2)==0 
    halfStates = numStates/2;
    for i=1:halfStates 
        fret_val = FRETvalues_nodes(i);
        prec = ifElse(abs(fret_val - round(fret_val,1)) < 0.001 && (round(fret_val,1)~=0||fret_val==0), 1, 2);
        nodeNames{i}= sprintf(['B %.*f',newline,'(%.0f%%|%.0f%%)'], prec, fret_val, 100*population1(i), 100*population2(i));
    end
    for i=1:halfStates 
        idx = i + halfStates;
        fret_val = FRETvalues_nodes(idx);
        prec = ifElse(abs(fret_val - round(fret_val,1)) < 0.001 && (round(fret_val,1)~=0||fret_val==0), 1, 2);
        nodeNames{idx}= sprintf(['D %.*f',newline,'(%.0f%%|%.0f%%)'], prec, fret_val, 100*population1(idx), 100*population2(idx));
    end
else 
    for i=1:numStates
        fret_val_current = FRETvalues_nodes(i);
        if abs(fret_val_current - round(fret_val_current, 1)) < 0.001 && (round(fret_val_current,1) ~= 0 || fret_val_current == 0)
            precision = 1;
        else
            precision = 2;
        end
        nodeNames{i}= sprintf('%.*f\n(%.0f%%|%.0f%%)', ...
                              precision, fret_val_current, ...
                              100*population1(i), ...
                              100*population2(i));
    end
end
if numStates > 0, G.Nodes.Name = nodeNames; end

%% Node Sizes
NodeSizeVec = zeros(numStates, 1);
pop_for_sizing = population1;
if isrow(pop_for_sizing), pop_for_sizing = pop_for_sizing'; end
if strcmpi(graph_config.node_size_method, 'linear_pop_old_style')
    NodeSizeVec = graph_config.node_linear_scale_factor * pop_for_sizing;
else, NodeSizeVec = repmat(30, numStates, 1); end
NodeSizeVec(NodeSizeVec < graph_config.node_min_size_clamp) = graph_config.node_min_size_clamp;
NodeSizeVec(NodeSizeVec > graph_config.node_max_size_clamp) = graph_config.node_max_size_clamp;
NodeSizeVec(isnan(NodeSizeVec) | isinf(NodeSizeVec) | NodeSizeVec <=0) = graph_config.node_min_size_clamp;

%% Node Colors
NodeColorsMat = [];
if numStates > 0
    if numStates == 2 && strcmpi(graph_config.node_color_scheme, 'custom_for_2state_example')
        [~, fret_order] = sort(FRETvalues_nodes, 'descend');
        NodeColorsMat = NaN(numStates,3);
        NodeColorsMat(fret_order(1),:) = [0.8500 0.3250 0.0980]; 
        NodeColorsMat(fret_order(2),:) = [0    0.4470    0.7410]; 
    else
        default_colors = get(groot,'defaultAxesColorOrder');
        num_default_colors = size(default_colors, 1);
        color_indices = mod(0:numStates-1, num_default_colors) + 1;
        NodeColorsMat = default_colors(color_indices, :);
    end
else, NodeColorsMat = [0 0 0]; end

%% Edge Widths
EdgeWidthsVec = [];
if ~isempty(G.Edges) && ismember('Weight', G.Edges.Properties.VariableNames) && ~isempty(G.Edges.Weight)
    weights_for_width_scaling = abs(G.Edges.Weight);
    max_w = max(weights_for_width_scaling);
    if max_w > eps 
        EdgeWidthsVec = graph_config.edge_width_scale_factor * weights_for_width_scaling / max_w + graph_config.edge_min_width;
    else, EdgeWidthsVec = repmat(graph_config.edge_min_width, height(G.Edges), 1); end
    EdgeWidthsVec(EdgeWidthsVec > graph_config.edge_max_width) = graph_config.edge_max_width;
    EdgeWidthsVec(EdgeWidthsVec < graph_config.edge_min_width) = graph_config.edge_min_width;
    EdgeWidthsVec(isnan(EdgeWidthsVec) | isinf(EdgeWidthsVec)) = graph_config.edge_min_width;
else
    fprintf('HMMgraph DEBUG: EdgeWidthsVec calculation skipped or defaulted because G.Edges.Weight is missing/empty.\n');
    EdgeWidthsVec = graph_config.edge_min_width; 
end


%% Plot the graph
custom_layout_applied = false;
plot_options = {};
if numStates == 3 && strcmpi(graph_config.layout_style, 'force')
    [~, fret_order_layout] = sort(FRETvalues_nodes, 'descend');
    x_coords = NaN(1, numStates); y_coords = NaN(1, numStates);
    x_coords(fret_order_layout(1)) = 0;     y_coords(fret_order_layout(1)) = 1;
    x_coords(fret_order_layout(2)) = -0.866;y_coords(fret_order_layout(2)) = -0.5;
    x_coords(fret_order_layout(3)) = 0.866; y_coords(fret_order_layout(3)) = -0.5;
    plot_options = {'XData', x_coords, 'YData', y_coords};
    custom_layout_applied = true;
    fprintf('HMMgraph: Applied custom triangular layout for 3 states.\n');
elseif numStates == 4 && strcmpi(graph_config.layout_style, 'force')
    [~, fret_order_layout] = sort(FRETvalues_nodes, 'descend');
    x_coords = NaN(1, numStates); y_coords = NaN(1, numStates);
    x_coords(fret_order_layout(1)) = -0.7; y_coords(fret_order_layout(1)) = 0.7;
    x_coords(fret_order_layout(2)) = 0.7;  y_coords(fret_order_layout(2)) = 0.7;
    x_coords(fret_order_layout(3)) = 0.7;  y_coords(fret_order_layout(3)) = -0.7;
    x_coords(fret_order_layout(4)) = -0.7; y_coords(fret_order_layout(4)) = -0.7;
    plot_options = {'XData', x_coords, 'YData', y_coords};
    custom_layout_applied = true;
    fprintf('HMMgraph: Applied custom square layout for 4 states.\n');
else
    plot_options = {'Layout', graph_config.layout_style, ...
                    'Iterations', graph_config.force_layout_iterations, ...
                    'WeightEffect', graph_config.force_layout_weight_effect};
end

p = plot(G, plot_options{:});

p.Marker = 'o';
if numStates > 0
    p.NodeColor = NodeColorsMat;
    p.MarkerSize = NodeSizeVec;
end
p.NodeFontSize = graph_config.node_text_font_size;
p.NodeFontWeight = graph_config.node_text_font_weight;
p.NodeLabelColor = graph_config.node_text_color;

fprintf('HMMgraph DEBUG: Just before styling edges (inside plot section):\n');
fprintf('HMMgraph DEBUG: isempty(G.Edges) = %d\n', isempty(G.Edges));
if ~isempty(G.Edges)
    fprintf('HMMgraph DEBUG: ismember("Weight", G.Edges.Properties.VariableNames) = %d\n', ismember('Weight', G.Edges.Properties.VariableNames));
    if ismember('Weight', G.Edges.Properties.VariableNames) && ~isempty(G.Edges.Weight)
        fprintf('HMMgraph DEBUG: numel(G.Edges.Weight) = %d\n', numel(G.Edges.Weight));
        fprintf('HMMgraph DEBUG: Any NaN/Inf in G.Edges.Weight? %d (NaN), %d (Inf)\n', any(isnan(G.Edges.Weight)), any(isinf(G.Edges.Weight)));
        fprintf('HMMgraph DEBUG: G.Edges.Weight values:\n'); disp(G.Edges.Weight);
    elseif ismember('Weight', G.Edges.Properties.VariableNames) && isempty(G.Edges.Weight)
         fprintf('HMMgraph DEBUG: G.Edges.Weight column exists but is EMPTY.\n');
    end
else
    fprintf('HMMgraph DEBUG: G.Edges is EMPTY just before styling (inside plot section).\n');
end

% --- CORRECTED CONDITION FOR STYLING EDGES ---
if ~isempty(G.Edges) && ismember('Weight', G.Edges.Properties.VariableNames) && ~isempty(G.Edges.Weight) && numel(G.Edges.Weight) > 0 
    p.LineWidth = EdgeWidthsVec;
    p.ArrowSize = graph_config.arrow_size;
    p.ArrowPosition = graph_config.arrow_position;
    p.EdgeAlpha = 1.0;

    fprintf('HMMgraph DEBUG: graph_config.show_edge_labels = %d\n', graph_config.show_edge_labels);
    if graph_config.show_edge_labels
        edgeLabelsForPlot = cell(height(G.Edges),1);
        for ie = 1:height(G.Edges)
            s_node = G.Edges.EndNodes(ie,1); t_node = G.Edges.EndNodes(ie,2);
            edgeLabelsForPlot{ie} = formatted_edge_labels{s_node, t_node};
        end
        p.EdgeLabel = edgeLabelsForPlot;
        p.EdgeFontSize = graph_config.edge_label_font_size;
        p.EdgeFontWeight = graph_config.edge_label_font_weight;
        p.EdgeLabelColor = graph_config.edge_label_color; 
        fprintf('HMMgraph DEBUG: Edge labels ASSIGNED to plot object.\n');
    else
        p.EdgeLabel = {};
        fprintf('HMMgraph DEBUG: Edge labels SKIPPED due to graph_config.show_edge_labels = false.\n');
    end
    
    fprintf('HMMgraph DEBUG: graph_config.edge_color_mode = %s\n', graph_config.edge_color_mode);
    fprintf('HMMgraph DEBUG: graph_config.fixed_edge_color = [%.1f %.1f %.1f]\n', graph_config.fixed_edge_color);
    if strcmpi(graph_config.edge_color_mode, 'fixed')
        p.EdgeColor = graph_config.fixed_edge_color; 
        fprintf('HMMgraph DEBUG: p.EdgeColor SET to fixed color: [%.1f %.1f %.1f]\n', p.EdgeColor);
    else
        p.EdgeColor = [0.15 0.15 0.15]; 
        fprintf('HMMgraph DEBUG: p.EdgeColor set to default non-fixed color.\n');
    end
else
    fprintf('HMMgraph DEBUG: MAIN STYLING BLOCK SKIPPED. Conditions:\n');
    fprintf('  isempty(G.Edges): %d\n', isempty(G.Edges));
    if ~isempty(G.Edges)
        fprintf('  ismember("Weight", G.Edges.Properties.VariableNames): %d\n', ismember('Weight', G.Edges.Properties.VariableNames));
        if ismember('Weight', G.Edges.Properties.VariableNames)
             fprintf('  isempty(G.Edges.Weight): %d\n', isempty(G.Edges.Weight));
             fprintf('  numel(G.Edges.Weight) > 0: %d\n', numel(G.Edges.Weight) > 0);
        end
    end
end

if ~custom_layout_applied && numStates == 2 && strcmpi(graph_config.layout_style, 'force')
    xData = p.XData; yData = p.YData;
    [~, fret_order] = sort(FRETvalues_nodes, 'descend');
    top_node_idx = fret_order(1); bottom_node_idx = fret_order(2);

    if abs(yData(top_node_idx) - yData(bottom_node_idx)) < abs(xData(top_node_idx) - xData(bottom_node_idx)) * 0.9 
        p.XData = yData; p.YData = xData;
        if p.YData(top_node_idx) < p.YData(bottom_node_idx), p.YData = -p.YData; end
    else 
        if yData(top_node_idx) < yData(bottom_node_idx), p.YData = -yData; end 
    end
    p.XData = p.XData - mean(p.XData); 
    
    x_diff = abs(p.XData(1) - p.XData(2));
    y_diff = abs(p.YData(1) - p.YData(2));
    if y_diff > x_diff * 5 && x_diff < 0.1 
        nudge = 0.05 + 0.05 * (max(NodeSizeVec)/50); 
        p.XData(top_node_idx) = p.XData(top_node_idx) - nudge;
        p.XData(bottom_node_idx) = p.XData(bottom_node_idx) + nudge;
    end
end

axis tight;
current_xlim = xlim; x_range = range(current_xlim); if x_range < 0.5, x_range=0.5; end
current_ylim = ylim; y_range = range(current_ylim); if y_range < 1, y_range=1; end

if custom_layout_applied
    axis_padding_factor = 0.40; 
else
    axis_padding_factor = 0.25;
end
xlim(current_xlim + [-x_range*axis_padding_factor x_range*axis_padding_factor]);
ylim(current_ylim + [-y_range*axis_padding_factor y_range*axis_padding_factor]);

daspect([1 1 1]);
title(''); axis off; box off; set(f, 'Color', 'w');

output_basename = 'HMM_State_Model_RatesPopBlackArrows';
try
    saveas(f, [output_basename '.fig']);
    print(f, [output_basename '.png'], '-dpng', '-r300');
    fprintf('HMMgraph: Saved graph as %s.fig and %s.png\n', output_basename, output_basename);
catch ME_save_graph
    warning('HMMgraph: Could not save graph. Error: %s', getReport(ME_save_graph, 'basic', 'hyperlinks', 'off'));
end

end

function C = ifElse(condition, A, B)
    if condition, C = A; else, C = B; end
end
% --- END OF FILE HMMgraph_local.m ---