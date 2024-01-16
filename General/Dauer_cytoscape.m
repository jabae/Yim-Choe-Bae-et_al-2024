
clc;clear;
format compact;

%% Data import

cd C:\Users\WORM-05\Documents\MATLAB\excel\

data1 = readmatrix('total_connection_arranged_new.xlsx','Sheet','dataset1','Range','C3:HP224');
data2 = readmatrix('total_connection_arranged_new.xlsx','Sheet','dataset2','Range','C3:HP224');
data3 = readmatrix('total_connection_arranged_new.xlsx','Sheet','dataset3','Range','C3:HP224');
data4 = readmatrix('total_connection_arranged_new.xlsx','Sheet','dataset4','Range','C3:HP224');
data5 = readmatrix('total_connection_arranged_new.xlsx','Sheet','dataset5','Range','C3:HP224');
data6 = readmatrix('total_connection_arranged_new.xlsx','Sheet','dataset6','Range','C3:HP224');
data7 = readmatrix('total_connection_arranged_new.xlsx','Sheet','dataset7','Range','C3:HP224');
data8 = readmatrix('total_connection_arranged_new.xlsx','Sheet','dataset8','Range','C3:HP224');
data_dauer = readmatrix('total_connection_arranged_new.xlsx','Sheet','dauer_221021','Range','C3:HP224');

data_JSH = readmatrix('total_connection_arranged_new.xlsx','Sheet','JSH','Range','C3:HP224');
data_N2U = readmatrix('total_connection_arranged_new.xlsx','Sheet','N2U','Range','C3:HP224');



% SD_cutoff

cutoff = 0.05;
data1_sd_cutoff = cutoff_adj(data1, cutoff);
data2_sd_cutoff = cutoff_adj(data2, cutoff);
data3_sd_cutoff = cutoff_adj(data3, cutoff);
data4_sd_cutoff = cutoff_adj(data4, cutoff);
data5_sd_cutoff = cutoff_adj(data5, cutoff);
data6_sd_cutoff = cutoff_adj(data6, cutoff);
data7_sd_cutoff = cutoff_adj(data7, cutoff);
data8_sd_cutoff = cutoff_adj(data8, cutoff);
data_JSH_sd_cutoff = cutoff_adj(data_JSH, cutoff);
data_N2U_sd_cutoff = cutoff_adj(data_N2U, cutoff);
data_dauer_sd_cutoff = cutoff_adj(data_dauer, cutoff);



%% mini matrix

data1_mini = mini_matrix(data1_sd_cutoff, 'class');
data2_mini = mini_matrix(data2_sd_cutoff, 'class');
data3_mini = mini_matrix(data3_sd_cutoff, 'class');
data4_mini = mini_matrix(data4_sd_cutoff, 'class');
data5_mini = mini_matrix(data5_sd_cutoff, 'class');
data6_mini = mini_matrix(data6_sd_cutoff, 'class');
data7_mini = mini_matrix(data7_sd_cutoff, 'class');
data8_mini = mini_matrix(data8_sd_cutoff, 'class');
data_N2U_mini = mini_matrix(data_N2U_sd_cutoff, 'class');
data_dauer_mini = mini_matrix(data_dauer_sd_cutoff, 'class');


% submini matrix
data1_submini = mini_matrix(data1_sd_cutoff, 'subclass');
data2_submini = mini_matrix(data2_sd_cutoff, 'subclass');
data3_submini = mini_matrix(data3_sd_cutoff, 'subclass');
data4_submini = mini_matrix(data4_sd_cutoff, 'subclass');
data5_submini = mini_matrix(data5_sd_cutoff, 'subclass');
data6_submini = mini_matrix(data6_sd_cutoff, 'subclass');
data7_submini = mini_matrix(data7_sd_cutoff, 'subclass');
data8_submini = mini_matrix(data8_sd_cutoff, 'subclass');
data_N2U_submini = mini_matrix(data_N2U_sd_cutoff, 'subclass');
data_dauer_submini = mini_matrix(data_dauer_sd_cutoff, 'subclass');




%% ------------ Circuit comparing ------------


d_dauer = data_dauer_sd_cutoff;
d_8 = data8_sd_cutoff;



%% ASG circuit


name_cell = {'ASG' 'ASI' 'AIA'};
category = 'class';
[circuit_dauer, circuit_list_dauer] = circuit_compare(d_dauer, name_cell, category);
[circuit_data8, circuit_list_data8] = circuit_compare(d_8, name_cell, 'class');

cytoscape_cell = cytoscape_cutoff_circuit(circuit_data8, circuit_dauer, circuit_list_data8, circuit_list_dauer, category, 0);

% export the result of cytoscape_cell to excel sheet



%% RIC circuit


name_cell = {'RIC' 'AVA' 'AVB' 'IL1'};
category = 'class';
[circuit_dauer, circuit_list_dauer] = circuit_compare(d_dauer, name_cell, category);
[circuit_data8, circuit_list_data8] = circuit_compare(d_8, name_cell, 'class');

cytoscape_cell = cytoscape_cutoff_circuit(circuit_data8, circuit_dauer, circuit_list_data8, circuit_list_dauer, category, 0);

% export the result of cytoscape_cell to excel sheet




%% ------------ Node comparing ------------


%% node 비교하기- in / out diagram 그리기. cytoscape_cutoff하고 x,y도 넣어서

neuron_id = 44;  
category = 'class';
[data8_cell, data8_list] = node_compare(data8_mini, neuron_id, category); 
[dauer_cell, dauer_list] = node_compare(data6_mini, neuron_id, category);

[cytoscape_cell, cytoscape_node] = cytoscape_cutoff_node(data8_cell, dauer_cell, data8_list, dauer_list, category, 0.2728, neuron_id);

% export the result of cytoscape_cell and cytoscape_node to excel sheet
% cytoscape_node : x, y position of nodes










%% ------------ functions ------------


%% basic functions


function edge_list = adj2list(adj)

num_node = size(adj,1);
node_list = [1:num_node];

matrix_sz = [num_node, num_node];
edge_list = [];

for ii = 1 : num_node^2
    
    temp = adj(ii);
    if temp ~= 0
      [pre, post] = ind2sub(matrix_sz,ii);
      
      temp = [pre, post];
      edge_list = [edge_list; temp];

    end
    
end

    
edge_list_1 = edge_list(:,1);
[~, idx] = sort(edge_list_1);
edge_list_pre = edge_list(idx, :);
edge_list = edge_list_pre;

end





function edge_list = adj2list_post(adj)


num_node = size(adj,1);
node_list = [1:num_node];

matrix_sz = [num_node, num_node];
edge_list = [];



for ii = 1 : num_node^2
    
    temp = adj(ii);
    if temp ~= 0
      [pre, post] = ind2sub(matrix_sz,ii);
      
     
      temp = [pre, post];
      edge_list = [edge_list; temp];

    end
    
end

end





function edge_list_weight = adj2list_weight(adj)

num_node = size(adj,1);
node_list = [1:num_node];

matrix_sz = [num_node, num_node];
edge_list_weight = [];

for ii = 1 : num_node^2
    
    temp = adj(ii);
    if temp ~= 0
      [pre, post] = ind2sub(matrix_sz,ii);
      
      temp2 = [pre, post, temp];
      edge_list_weight = [edge_list_weight; temp2];

    end
    
end

    
edge_list_1 = edge_list_weight(:,1);
[~, idx] = sort(edge_list_1);
edge_list_pre = edge_list_weight(idx, :);
edge_list_weight = edge_list_pre;

end



function edge_list_weight = row2list(adj)

num_node = length(adj);
node_list = [1:num_node];
matrix_sz = size(adj);
edge_list_weight = [];

for ii = 1 : num_node
    
    temp = adj(ii);
    if temp ~= 0
      [pre, post] = ind2sub(matrix_sz,ii);

      
      temp2 = [pre, post, temp];
      edge_list_weight = [edge_list_weight; temp2];

    end
    
end


if ~isempty(edge_list_weight)

    edge_list_1 = edge_list_weight(:,1);
    [~, idx] = sort(edge_list_1);
    edge_list_pre = edge_list_weight(idx, :);
    edge_list_weight = edge_list_pre;
    
end


end






function edge_cell = adj2cell(adj)


num_node = size(adj,1);
node_list = [1:num_node];

matrix_sz = [num_node, num_node];
edge_list = [];
edge_cell = {};

for ii = 1 : num_node
    edge_cell{ii} = [];
end

for ii = 1 : num_node^2
    
    temp = adj(ii);
    if temp ~= 0
      [pre, post] = ind2sub(matrix_sz,ii);

      temp = [pre, post];
      edge_list = [edge_list; temp];

      edge_cell{pre} = [edge_cell{pre}, post];   
      
    end
    
   
    
end

end




function edge_cell = adj2cell_post(adj)


num_node = size(adj,1);
node_list = [1:num_node];

matrix_sz = [num_node, num_node];
edge_list = [];
edge_cell = {};

for ii = 1 : num_node
    edge_cell{ii} = [];
end

for ii = 1 : num_node^2
    
    temp = adj(ii);
    if temp ~= 0
      [pre, post] = ind2sub(matrix_sz,ii);
      
      temp = [pre, post];
      edge_list = [edge_list; temp];

      edge_cell{post} = [edge_cell{post}, pre];   
      
    end
    
   
   
end

end



function adj = list2adj(edge_list, num_node)

num_edge = size(edge_list,1);
adj = zeros(num_node);

for ii = 1 : num_edge
    temp = edge_list(ii,:);
    pre = temp(1);
    post = temp(2);
    adj(pre,post) = 1;
end

end



function adj_weight = list2adj_weight(edge_list_weight, num_node)

num_edge = size(edge_list_weight,1);
adj_weight = zeros(num_node);

for ii = 1 : num_edge
    temp = edge_list_weight(ii,:);
    pre = temp(1);
    post = temp(2);
    weight = temp(3);
    adj_weight(pre,post) = adj_weight(pre,post) + weight;
end

end




function adj = cell2adj(edge_cell, num_node)

adj = zeros(num_node);

for ii = 1 : num_node
    temp = edge_cell{ii};
    num_temp = length(temp);
    
    for jj = 1 : num_temp
        post = temp(jj);
        adj(ii,post) = 1;
    end 
    
end

end



function adj = inc2adj(s,t, num_node)

num_edge = length(s);
adj = zeros(num_node);

for ii = 1 : num_edge
    pre = s(ii);
    post = t(ii);
    adj(pre,post) = 1;
end

end



function edge_list_post = postarrange(edge_list)

temp = edge_list(:,2);
[~, idx] = sort(temp);
edge_list_post = edge_list;

for ii = 1 : size(edge_list, 1)
    temp_idx = idx(ii);
    edge_list_post(ii,:) = edge_list(temp_idx,:);
end

end




function weight_list = list2weight(edge_list, sheet_name)


adj = readmatrix('total_synapse_arranged_sd.xlsx','Sheet', sheet_name,'Range','C3:IH242');
weight_list = zeros(size(edge_list, 1), 3);

for ii = 1 : size(edge_list, 1);
    
    pre = edge_list(ii,1);
    post = edge_list(ii,2);
    weight = adj(pre,post);
    weight_list(ii,:) = [pre post weight];

end

end



function weight_arranged = weight_arrange(edge_list)

weight = edge_list(:,3);
[~, idx] = sort(weight, 'descend');
weight_arranged = edge_list(idx, :);


end




% cutoff matrix 

function new_adj = cutoff_adj(adj, cutoff)

    num_node = size(adj, 1);


    temp_cutoff = 1-cutoff;
    list = adj2list_weight(adj);
    weight = list(:,3);
    weight_sort = sort(weight);
    re_length = round(length(weight_sort)* temp_cutoff);
    new_weight = weight_sort(1:re_length);
    sd = std(new_weight,0,'all'); 
    weight_sd = weight/sd;
    
    pre = list(:,1);
    post= list(:,2);
    new_list = [pre post weight_sd];
    
    new_adj = list2adj_weight(new_list, num_node);

end




% make mini matrix by clustering class or subclass

function new_matrix = mini_matrix(adj, category)


cd C:\Users\WORM-05\Documents\MATLAB\excel\
key = readtable('segmentation_key_common_230926.csv');

common_ID = key.Common_ID;

if strcmp(category, 'class')
    code = key.class_code;
elseif strcmp(category, 'subclass')
    code = key.subclass_code;
else
    fprintf("Please write down 'class' or 'subclass'")
end
    
unique_code = unique(code);
new_matrix = zeros(length(unique_code));

for ii = 1 : 222
    comm_pre_id = ii;
    new_pre_id = code(ii);

    for jj = 1 : 222
        comm_post_id = jj;
        new_post_id = code(jj);
        
        temp_weight = adj(ii,jj);
        
        new_matrix(new_pre_id, new_post_id) =  new_matrix(new_pre_id, new_post_id) + temp_weight;    

    end
        

end


end







function name_cell = id2name_mini(edge_list, category)

cd C:\Users\WORM-05\Documents\MATLAB\excel\
key = readtable('segmentation_key_subclass_230724.csv');

if strcmp(category, 'class')
  neuron_id = key.Hobert_class;
elseif strcmp(category, 'subclass')
  neuron_id = key.Hobert_subclass;
else
  fprintf("Please write down 'class' or 'subclass'")
end

name_cell = cell(size(edge_list,1), size(edge_list,2));

for ii = 1: size(edge_list,1)
    
    pre = edge_list(ii,1);
    
    if size(edge_list,2) > 1
        post = edge_list(ii,2);
        if size(edge_list,2) == 3
            weight = edge_list(ii,3);
        end
    end

    pre_id = neuron_id{pre};
    if size(edge_list,2) > 1
        post_id = neuron_id{post};
    end
    
    if size(edge_list,2) == 3
        fprintf('%s %s %d \n', pre_id, post_id, weight);  
        name_cell{ii, 1} = neuron_id{pre};
        name_cell{ii, 2} = neuron_id{post};
        name_cell{ii, 3} = weight;
    elseif size(edge_list,2) == 2
        fprintf('%s %s \n', pre_id, post_id);  
        name_cell{ii, 1} = pre_id;
        name_cell{ii, 2} = neuron_id{post};
    else
        fprintf('%s \n', pre_id);  
        name_cell{ii, 1} = pre_id;        
    end
    
end

end




% name_cell to numbered id

function new_list = name2id_mini(name_cell, category) 

cd C:\Users\WORM-05\Documents\MATLAB\excel\
key = readtable('segmentation_key_subclass_230724.csv');

if strcmp(category, 'class')
  neuron_id = key.Hobert_class;
  code = key.class_code;
elseif strcmp(category, 'subclass')
  neuron_id = key.Hobert_subclass;
  code = key.subclass_code;  
else
  fprintf("Please write down 'class' or 'subclass'")
end


new_list = [];

for ii = 1 : length(code)
    check_id = neuron_id(ii);
  
    for jj = 1 : length(name_cell);
        check_source = name_cell{jj}; 
        
        if strcmp(check_id,check_source)
            new_list = [new_list; ii]; 
        end
    end

end

end





%% functions for compare


function [new_name_cell, new_matrix_list] = circuit_compare(adj, name_cell, category)

% adj to new_matrix
new_matrix = mini_matrix(adj, category);

% name_cell to numbered id
new_list = name2id_mini(name_cell, category) ;


cd C:\Users\WORM-05\Documents\MATLAB\excel\
key = readtable('segmentation_key_subclass_230724.csv');

% if you write down numbered id_list, make adj matrix
if strcmp(category, 'class')
  node_length = length(key.class_code);
  nodes = [1:node_length];
elseif strcmp(category, 'subclass')
  node_length = length(key.subclass_code);
  nodes = [1:node_length];
else
  fprintf("Please write down 'class' or 'subclass'")
end

other_list = setdiff(nodes, new_list);
new_matrix(other_list, :) = 0;
new_matrix(:, other_list) = 0;

% mini matrix to list
new_matrix_list = adj2list_weight(new_matrix);

% id2name(subclass neuron)
new_name_cell = id2name_mini(new_matrix_list, category);

end





% write neuron subclass code -> type_code

function [type_code, type_name] = type_subclass(check_pre)

data_sub = readtable('segmentation_key_common_230926.csv');

for jj = 1 : length(data_sub.Hobert_subclass)
    temp_check = check_pre(1:3);
    temp_compare = data_sub.Hobert_subclass{jj};
    temp_compare = temp_compare(1:3);
    if strcmp(temp_check, temp_compare);
        check = data_sub.Common_ID(jj);
        break
    end   
end


data_type = readtable('230724_neurontype_revised.csv');
neuron_type = data_type.Cook_typecode2;  % muscle categorize


% analysis neuron types 

node_list = [1:222];

tf_s = neuron_type == 1;
tf_i = neuron_type == 2;
tf_m = neuron_type == 3;
tf_mus = neuron_type == 4;
tf_o = neuron_type == 0;

neuron_sensory = node_list(tf_s);
neuron_inter = node_list(tf_i);
neuron_moter = node_list(tf_m);
muscle = node_list(tf_mus);
endorgan = node_list(tf_o);

if sum(neuron_sensory == check, 'all')
    type_code = 1;
    type_name = 'sensory';
elseif sum(neuron_inter == check, 'all')
    type_code = 2;
    type_name = 'inter';
elseif sum(neuron_moter == check, 'all')
    type_code = 3;
    type_name = 'motor';
elseif sum(neuron_inter == check, 'all')
    type_code = 4;
    type_name = 'muscle';
else
    type_code = 0;
    type_name = 'endorgan';
end

end







%% functions for id transformation

% if yopu write list(common id) -> neuron id
% input: n * 1,2 (only id) or n * 3 (with weight) array

function name_cell = id2name(edge_list)

cd C:\Users\WORM-05\Documents\MATLAB\excel\
data = readtable('segmentation_key_common_230926.csv');
neuron_id = data.Neuron;

name_cell = cell(size(edge_list,1), size(edge_list,2));

for ii = 1: size(edge_list,1)
    
    pre = edge_list(ii,1);
    
    if size(edge_list,2) > 1
        post = edge_list(ii,2);
        if size(edge_list,2) == 3
            weight = edge_list(ii,3);
        end
    end

    pre_id = neuron_id{pre};
    if size(edge_list,2) > 1
        post_id = neuron_id{post};
    end
    
    if size(edge_list,2) == 3
        fprintf('%s %s %d \n', pre_id, post_id, weight);  
        name_cell{ii, 1} = neuron_id{pre};
        name_cell{ii, 2} = neuron_id{post};
        name_cell{ii, 3} = weight;
    elseif size(edge_list,2) == 2
        fprintf('%s %s \n', pre_id, post_id);  
        name_cell{ii, 1} = pre_id;
        name_cell{ii, 2} = neuron_id{post};
    else
        fprintf('%s \n', pre_id);  
        name_cell{ii, 1} = pre_id;        
    end
    
end

end





% if you write neuron's alphabet name -> common id
% input : cell array
%{
example:
name_cell = {'ADAL' 'ADAR'; 'VC2' 'RIH'; 'AVAL' 'AVAR';};
edge_list = name2id(name_cell);
%}

function new_list = name2id(name_cell)  

cd C:\Users\WORM-05\Documents\MATLAB\excel\
data = readtable('segmentation_key_common_230926.csv');
neuron_id = data.Neuron;

new_list = NaN(size(name_cell,1), size(name_cell,2));

for ii = 1 : 222
    
    check_id = neuron_id(ii);

    for jj = 1 : size(name_cell,1);
        for kk = 1 : size(name_cell, 2); 
            check_source = name_cell{jj,kk};
            if strcmp(check_id,check_source)
                new_list(jj,kk) = ii;    
            end
        end
    end
   
end

end






function [total_cell, total_list, output_list_cell, input_list_cell, output_list, input_list] = node_compare (adj, neuron_id, style)



% write down numbered_id list > adj matrix

if strcmp(style, 'class')

    row = adj(neuron_id, :);
    col = adj(:, neuron_id);
    output_list = row2list(row);  % output
    input_list = row2list(col);  % input
    
    category = 'class';

    row = adj(neuron_id, :);
    col = adj(:, neuron_id);
    output_list = row2list(row);  % output
    input_list = row2list(col);  % input
    if ~isempty(output_list)
        if ~isempty(input_list)
            output_list(:,1) = neuron_id;
            input_list(:,2) = neuron_id;
            total_list = [output_list; input_list];

            % subclass neuron's id2name
            output_list_cell = id2name_mini(output_list, category);
            input_list_cell = id2name_mini(input_list, category);
            total_cell = id2name_mini(total_list, category);
        
        else
            output_list(:,1) = neuron_id;
            total_list = [output_list; input_list];
            output_list_cell = id2name_mini(output_list, category);
            total_cell = id2name_mini(total_list, category);
        end
        
    else
        if ~isempty(input_list)
            input_list(:,2) = neuron_id;
            total_list = [output_list; input_list];
            input_list_cell = id2name_mini(output_list, category);
            total_cell = id2name_mini(total_list, category);
        else
            total_list = [];
            total_cell = [];
        end
    end
    
elseif strcmp(style, 'subclass')
    
    category = 'subclass';

    row = adj(neuron_id, :);
    col = adj(:, neuron_id);
    output_list = row2list(row);  % output
    input_list = row2list(col);  % input
    if ~isempty(output_list)
        if ~isempty(input_list)
            output_list(:,1) = neuron_id;
            input_list(:,2) = neuron_id;
            total_list = [output_list; input_list];

            % subclass neuron's id2name
            output_list_cell = id2name_mini(output_list, category);
            input_list_cell = id2name_mini(input_list, category);
            total_cell = id2name_mini(total_list, category);
        
        else
            output_list(:,1) = neuron_id;
            total_list = [output_list; input_list];
            output_list_cell = id2name_mini(output_list, category);
            total_cell = id2name_mini(total_list, category);
        end
        
    else
        if ~isempty(input_list)
            input_list(:,2) = neuron_id;
            total_list = [output_list; input_list];
            input_list_cell = id2name_mini(output_list, category);
            total_cell = id2name_mini(total_list, category);
        else
            total_list = [];
            total_cell = [];
        end
    end
    
else
  fprintf("Please write down 'class' or 'subclass'")
end


end












%% functions for cytoscape export



% ex. (data8_cell, dauer_cell, data8_list, dauer_list, 'class', 0.2);

function [cytoscape_cell, node_cell] = cytoscape_cutoff_circuit(nondauer_circuit, dauer_circuit, nondauer_list_full, dauer_list_full, style, cutoff)


nondauer_list = nondauer_list_full(:, 1:2);
dauer_list = dauer_list_full(:, 1:2);
union_list = union(nondauer_list, dauer_list,'rows');



% remove self connections
self_list = [];
for ii = 1 : size(union_list, 1);
    pre = union_list(ii,1);
    post = union_list(ii,2);
    if pre == post
        self_list = [self_list; ii];
    end
    
    
    
    

% remove CEPsh
if strcmp(style, 'whole')
    if post == 220
       self_list = [self_list; ii];
    elseif post == 221
        self_list = [self_list; ii];
    elseif post == 222
        self_list = [self_list; ii];
    elseif post == 223
        self_list = [self_list; ii];
    end
elseif strcmp(style, 'class')
    if post == 85
       self_list = [self_list; ii];
    end
elseif strcmp(style, 'subclass')
    
    cd C:\Users\WORM-05\Documents\MATLAB\excel\
    key = readtable('segmentation_key_common_230926.csv');
    subclass_code = key.subclass_code;

    if post == subclass_code(220)
       self_list = [self_list; ii];
    end
    
else
    fprintf("Please write down 'whole' or 'class' 'subclass'")   ; 
end


end
union_list(self_list, :) = [];

    
    
    
    
    

% common edge or dauer/nondauer specific 
dauer_only_idx = [];
nondauer_only_idx = [];
common_idx = [];
for ii = 1 : size(union_list, 1);
    temp_row = union_list(ii, :);
    [~,~,i_dauer] = intersect(temp_row, dauer_list, 'rows');
    [~,~,i_nondauer] = intersect(temp_row, nondauer_list, 'rows');
    if i_dauer         
        if i_nondauer   
            common_idx = [common_idx; ii];
        else           
            dauer_only_idx = [dauer_only_idx; ii];
        end
    else               
        nondauer_only_idx = [nondauer_only_idx; ii];
    end
end
dauer_idx = union(dauer_only_idx, common_idx);
nondauer_idx = union(nondauer_only_idx, common_idx);


cytoscape_cell = cell(size(union_list, 1), 9);

for ii = 1 : size(union_list, 1);
    for jj = 1 : 9
        cytoscape_cell{ii,jj} = 0;
    end
end




% pre / post name
for ii = 1 : size(union_list, 1);
    temp_row = union_list(ii, :);
    [~,~,i_dauer] = intersect(temp_row, dauer_list, 'rows');
    [~,~,i_nondauer] = intersect(temp_row, nondauer_list, 'rows');
    if i_dauer   
        pre = dauer_circuit{i_dauer,1};
        post = dauer_circuit{i_dauer,2};
        [~, type_name] = type_subclass(pre);
        cytoscape_cell{ii, 1} = type_name;        
        cytoscape_cell{ii, 2} = pre;
        cytoscape_cell{ii, 3} = post;               
    else          
        pre = nondauer_circuit{i_nondauer,1};
        post = nondauer_circuit{i_nondauer,2};
        [~, type_name] = type_subclass(pre);
        cytoscape_cell{ii, 1} = type_name;          
        cytoscape_cell{ii, 2} = pre;
        cytoscape_cell{ii, 3} = post;       
    end
end



for ii = 1 : size(union_list, 1)
    cytoscape_cell{ii, 4} = 'weight_dauer';
end
for ii = 1 : size(union_list, 1)
    cytoscape_cell{ii, 5} = 'weight_nondauer';
end



for ii = 1 : size(union_list, 1);
    temp_row = union_list(ii, :);
    [~,~,i_dauer] = intersect(temp_row, dauer_list, 'rows');
    [~,~,i_nondauer] = intersect(temp_row, nondauer_list, 'rows');
    dauer_weight = dauer_list_full(i_dauer,3);
    nondauer_weight = nondauer_list_full(i_nondauer,3);
    if i_dauer          
        if i_nondauer   
            cytoscape_cell{ii, 6} = dauer_weight;
            cytoscape_cell{ii, 7} = nondauer_weight;
            cytoscape_cell{ii, 8} = dauer_weight/nondauer_weight;
            cytoscape_cell{ii, 9} = nondauer_weight/dauer_weight;   
        else          
            cytoscape_cell{ii, 6} = dauer_weight;
            cytoscape_cell{ii, 8} = 100;     
            cytoscape_cell{ii, 9} = 0; 
        end
    else              
        cytoscape_cell{ii, 7} = nondauer_weight;
        cytoscape_cell{ii, 8} = 0;  
        cytoscape_cell{ii, 9} = 100;  
    end

    
end





% cutoff

weight_dauer = [];
weight_nondauer = [];

for ii = 1 : size(cytoscape_cell,1);
    
    temp_dauer = cytoscape_cell{ii,6};
    temp_nondauer = cytoscape_cell{ii,7};
    if isempty(temp_dauer)
        temp_dauer = 0;
    end
    if isempty(temp_nondauer)
        temp_nondauer = 0;
    end   

    weight_dauer = [weight_dauer; temp_dauer];
    weight_nondauer = [weight_nondauer; temp_nondauer];

end

idx_dauer = find(weight_dauer < cutoff);
idx_nondauer = find(weight_nondauer < cutoff);
idx_remove = intersect(idx_dauer, idx_nondauer);
idx_survival = setdiff([1:size(cytoscape_cell,1)], idx_remove);

cytoscape_cell_new = cell(length(idx_survival), size(cytoscape_cell,2));

for ii = 1 : length(idx_survival)
    temp = idx_survival(ii);
    for jj = 1:size(cytoscape_cell,2)
        cytoscape_cell_new{ii,jj} = cytoscape_cell{temp,jj};
    end
end

cytoscape_cell = cytoscape_cell_new;


union_list(idx_remove, :) = [];


nodes = unique(union_list);



node_cell = cell(length(nodes),2);

for ii = 1 : length(nodes)
    
    check = nodes(ii);

    if strcmp(style, 'whole')
        data_type = readtable('230724_neurontype_revised.csv');
        subtype = data_type.Cook_type;
        node_cell{ii,1} = subtype(check);
        neuron_name = id2name(check);
        node_cell{ii,2} = neuron_name;
        
    elseif strcmp(style, 'class')
        data_type = readtable('230724_neurontype_mini.csv');
        subtype = data_type.Cook_type;
        node_cell{ii,1} = subtype(check);
        neuron_name = id2name_mini(check, 'class');
        node_cell{ii,2} = neuron_name;
        
    else
        fprintf("Please write down 'whole' or 'class', or type it manually\n")   ; 
    end

end



end








% ex. (data8_cell, dauer_cell, data8_list, dauer_list, 'mini', 0.2, neuron_id);

function [cytoscape_cell, node_cell] = cytoscape_cutoff_node(nondauer_circuit, dauer_circuit, nondauer_list_full, dauer_list_full, style, cutoff, neuron_id)




nondauer_list = nondauer_list_full(:, 1:2);
dauer_list = dauer_list_full(:, 1:2);
union_list = union(nondauer_list, dauer_list,'rows');



% remove self
self_list = [];
for ii = 1 : size(union_list, 1);
    pre = union_list(ii,1);
    post = union_list(ii,2);
    if pre == post
        self_list = [self_list; ii];
    end
    
    

% remove CEPsh
if strcmp(style, 'whole')
    if post == 220
       self_list = [self_list; ii];
    elseif post == 221
        self_list = [self_list; ii];
    elseif post == 222
        self_list = [self_list; ii];
    elseif post == 223
        self_list = [self_list; ii];
    end
elseif strcmp(style, 'class')
    if post == 85
       self_list = [self_list; ii];
    end
elseif strcmp(style, 'subclass')

    cd C:\Users\WORM-05\Documents\MATLAB\excel\
    key = readtable('segmentation_key_common_230926.csv');
    subclass_code = key.subclass_code;

    if post == subclass_code(220)
       self_list = [self_list; ii];
    end

    
else
    fprintf("Please write down 'whole' or 'class' or 'subclass\n'")   ; 
end


end
union_list(self_list, :) = [];

    
    
    
    
    
dauer_only_idx = [];
nondauer_only_idx = [];
common_idx = [];
for ii = 1 : size(union_list, 1);
    temp_row = union_list(ii, :);
    [~,~,i_dauer] = intersect(temp_row, dauer_list, 'rows');
    [~,~,i_nondauer] = intersect(temp_row, nondauer_list, 'rows');
    if i_dauer          
        if i_nondauer   
            common_idx = [common_idx; ii];
        else         
            dauer_only_idx = [dauer_only_idx; ii];
        end
    else              
        nondauer_only_idx = [nondauer_only_idx; ii];
    end
end
dauer_idx = union(dauer_only_idx, common_idx);
nondauer_idx = union(nondauer_only_idx, common_idx);

%
cytoscape_cell = cell(size(union_list, 1), 9);

for ii = 1 : size(union_list, 1);
    for jj = 1 : 9
        cytoscape_cell{ii,jj} = 0;
    end
end







for ii = 1 : size(union_list, 1);
    temp_row = union_list(ii, :);
    [~,~,i_dauer] = intersect(temp_row, dauer_list, 'rows');
    [~,~,i_nondauer] = intersect(temp_row, nondauer_list, 'rows');
    if i_dauer   
        pre = dauer_circuit{i_dauer,1};
        post = dauer_circuit{i_dauer,2};
        [~, type_name] = type_subclass(pre);
        cytoscape_cell{ii, 1} = type_name;        
        cytoscape_cell{ii, 2} = pre;
        cytoscape_cell{ii, 3} = post;               
    else         
        pre = nondauer_circuit{i_nondauer,1};
        post = nondauer_circuit{i_nondauer,2};
        [~, type_name] = type_subclass(pre);
        cytoscape_cell{ii, 1} = type_name;          
        cytoscape_cell{ii, 2} = pre;
        cytoscape_cell{ii, 3} = post;       
    end
end




for ii = 1 : size(union_list, 1);
    temp_row = union_list(ii, :);
    [~,~,i_dauer] = intersect(temp_row, dauer_list, 'rows');
    [~,~,i_nondauer] = intersect(temp_row, nondauer_list, 'rows');
    if i_dauer  
        pre = dauer_circuit{i_dauer,1};
        post = dauer_circuit{i_dauer,2};
        [~, type_name] = type_subclass(pre);
        cytoscape_cell{ii, 1} = type_name;        
        cytoscape_cell{ii, 2} = pre;
        cytoscape_cell{ii, 3} = post;               
    else         
        pre = nondauer_circuit{i_nondauer,1};
        post = nondauer_circuit{i_nondauer,2};
        [~, type_name] = type_subclass(pre);
        cytoscape_cell{ii, 1} = type_name;          
        cytoscape_cell{ii, 2} = pre;
        cytoscape_cell{ii, 3} = post;       
    end
end


for ii = 1 : size(union_list, 1)
    cytoscape_cell{ii, 4} = 'weight_dauer';
end
for ii = 1 : size(union_list, 1)
    cytoscape_cell{ii, 5} = 'weight_nondauer';
end


for ii = 1 : size(union_list, 1);
    temp_row = union_list(ii, :);
    [~,~,i_dauer] = intersect(temp_row, dauer_list, 'rows');
    [~,~,i_nondauer] = intersect(temp_row, nondauer_list, 'rows');

    if i_dauer      
        if i_nondauer   
            dauer_weight = dauer_list_full(i_dauer,3);
            nondauer_weight = nondauer_list_full(i_nondauer,3);
            cytoscape_cell{ii, 6} = dauer_weight;
            cytoscape_cell{ii, 7} = nondauer_weight;
            cytoscape_cell{ii, 8} = dauer_weight/nondauer_weight;
            cytoscape_cell{ii, 9} = nondauer_weight/dauer_weight;   
        else           
            dauer_weight = dauer_list_full(i_dauer,3);
            nondauer_weight = 0;
            cytoscape_cell{ii, 6} = dauer_weight;
            cytoscape_cell{ii, 8} = 100;     
            cytoscape_cell{ii, 9} = 0; 
        end
    else           
        dauer_weight = 0;
        nondauer_weight = nondauer_list_full(i_nondauer,3);        
        cytoscape_cell{ii, 7} = nondauer_weight;
        cytoscape_cell{ii, 8} = 0;  
        cytoscape_cell{ii, 9} = 100;  
    end

    
end








% cutoff

weight_dauer = [];
weight_nondauer = [];

for ii = 1 : size(cytoscape_cell,1);
    
    temp_dauer = cytoscape_cell{ii,6};
    temp_nondauer = cytoscape_cell{ii,7};
    if isempty(temp_dauer)
        temp_dauer = 0;
    end
    if isempty(temp_nondauer)
        temp_nondauer = 0;
    end   

    weight_dauer = [weight_dauer; temp_dauer];
    weight_nondauer = [weight_nondauer; temp_nondauer];

end

idx_dauer = find(weight_dauer < cutoff);
idx_nondauer = find(weight_nondauer < cutoff);
idx_remove = intersect(idx_dauer, idx_nondauer);
idx_survival = setdiff([1:size(cytoscape_cell,1)], idx_remove);

cytoscape_cell_new = cell(length(idx_survival), size(cytoscape_cell,2));

for ii = 1 : length(idx_survival)
    temp = idx_survival(ii);
    for jj = 1:size(cytoscape_cell,2)
        cytoscape_cell_new{ii,jj} = cytoscape_cell{temp,jj};
    end
end

cytoscape_cell = cytoscape_cell_new;





union_list(idx_remove, :) = [];




nodes = unique(union_list);

node_cell = cell(length(nodes),4);

num_node = length(nodes)-1;
radius = 200;  % variable
theta = 360/num_node;
timepoint = 0;

for ii = 1 : length(nodes)
    
    check = nodes(ii);

    if strcmp(style, 'whole')
        data_type = readtable('230724_neurontype_revised.csv');
        subtype = data_type.Cook_type;
        node_cell{ii,1} = subtype(check);
        neuron_name = id2name(check);
        node_cell{ii,2} = neuron_name;
        
    elseif strcmp(style, 'class')
        data_type = readtable('230724_neurontype_mini.csv');
        subtype = data_type.Cook_type;
        node_cell{ii,1} = subtype(check);
        neuron_name = id2name_mini(check, 'class');
        node_cell{ii,2} = neuron_name;

    elseif strcmp(style, 'subclass')
        data_type = readtable('230724_neurontype_mini.csv');
        subtype = data_type.Cook_subclass_type;
        node_cell{ii,1} = subtype(check);
        neuron_name = id2name_mini(check, 'subclass');
        node_cell{ii,2} = neuron_name;    
        
    else
        fprintf("Please write down 'whole' or 'class' or 'subclass'\n")   ; 
    end
    

    if check == neuron_id
        timepoint = ii;
        node_cell{ii,3} = 0;
        node_cell{ii,4} = 0;
    else
        if timepoint == 0
        temp_x = cosd(theta*(ii) - 90) * radius;
        temp_y =  sind(theta*(ii) - 90) * radius;
        node_cell{ii,3} = temp_x;
        node_cell{ii,4} = temp_y;
        else
        temp_x = cosd(theta*(ii-1) - 90) * radius;
        temp_y =  sind(theta*(ii-1) - 90) * radius;
        node_cell{ii,3} = temp_x;
        node_cell{ii,4} = temp_y;     
        end

    end
   
end









end










