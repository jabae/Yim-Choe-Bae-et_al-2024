
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




%% SD_cutoff

cutoff = 0.05;
[data1_sd_cutoff, sd1] = cutoff_adj(data1, cutoff);
[data2_sd_cutoff, sd2] = cutoff_adj(data2, cutoff);
[data3_sd_cutoff, sd3] = cutoff_adj(data3, cutoff);
[data4_sd_cutoff, sd4] = cutoff_adj(data4, cutoff);
[data5_sd_cutoff, sd5] = cutoff_adj(data5, cutoff);
[data6_sd_cutoff, sd6] = cutoff_adj(data6, cutoff);
[data7_sd_cutoff, sd7] = cutoff_adj(data7, cutoff);
[data8_sd_cutoff, sd8] = cutoff_adj(data8, cutoff);
[data_JSH_sd_cutoff, sd_JSH] = cutoff_adj(data_JSH, cutoff);
[data_N2U_sd_cutoff, sd_N2U] = cutoff_adj(data_N2U, cutoff);
[data_dauer_sd_cutoff, sd_dauer] = cutoff_adj(data_dauer, cutoff);


data1_list = adj2list_weight(data1_sd_cutoff);
data2_list = adj2list_weight(data2_sd_cutoff);
data3_list = adj2list_weight(data3_sd_cutoff);
data4_list = adj2list_weight(data4_sd_cutoff);
data5_list = adj2list_weight(data5_sd_cutoff);
data6_list = adj2list_weight(data6_sd_cutoff);
data7_list = adj2list_weight(data7_sd_cutoff);
data8_list = adj2list_weight(data8_sd_cutoff);
data_JSH_list = adj2list_weight(data_JSH_sd_cutoff);
data_N2U_list = adj2list_weight(data_N2U_sd_cutoff);
data_dauer_list = adj2list_weight(data_dauer_sd_cutoff);



%%

data1_neuron = data1_sd_cutoff(1:180, 1:180);
data2_neuron = data2_sd_cutoff(1:180, 1:180);
data3_neuron = data3_sd_cutoff(1:180, 1:180);
data4_neuron = data4_sd_cutoff(1:180, 1:180);
data5_neuron = data5_sd_cutoff(1:180, 1:180);
data6_neuron = data6_sd_cutoff(1:180, 1:180);
data7_neuron = data7_sd_cutoff(1:180, 1:180);
data8_neuron = data8_sd_cutoff(1:180, 1:180);
data_JSH_neuron = data_JSH_sd_cutoff(1:180, 1:180);
data_N2U_neuron = data_N2U_sd_cutoff(1:180, 1:180);
data_dauer_neuron = data_dauer_sd_cutoff(1:180, 1:180);


data1_neuron_list = adj2list_weight(data1_neuron);
data2_neuron_list = adj2list_weight(data2_neuron);
data3_neuron_list = adj2list_weight(data3_neuron);
data4_neuron_list = adj2list_weight(data4_neuron);
data5_neuron_list = adj2list_weight(data5_neuron);
data6_neuron_list = adj2list_weight(data6_neuron);
data7_neuron_list = adj2list_weight(data7_neuron);
data8_neuron_list = adj2list_weight(data8_neuron);
data_JSH_neuron_list = adj2list_weight(data_JSH_neuron);
data_N2U_neuron_list = adj2list_weight(data_N2U_neuron);
data_dauer_neuron_list = adj2list_weight(data_dauer_neuron);



%% mini matrix

data1_mini = mini_matrix(data1_sd_cutoff, 'class');
data2_mini = mini_matrix(data2_sd_cutoff, 'class');
data3_mini = mini_matrix(data3_sd_cutoff, 'class');
data4_mini = mini_matrix(data4_sd_cutoff, 'class');
data5_mini = mini_matrix(data5_sd_cutoff, 'class');
data6_mini = mini_matrix(data6_sd_cutoff, 'class');
data7_mini = mini_matrix(data7_sd_cutoff, 'class');
data8_mini = mini_matrix(data8_sd_cutoff, 'class');
data_JSH_mini = mini_matrix(data_JSH_sd_cutoff, 'class');
data_N2U_mini = mini_matrix(data_N2U_sd_cutoff, 'class');
data_dauer_mini = mini_matrix(data_dauer_sd_cutoff, 'class');



%% submini matrix

data1_submini = mini_matrix(data1_sd_cutoff, 'subclass');
data2_submini = mini_matrix(data2_sd_cutoff, 'subclass');
data3_submini = mini_matrix(data3_sd_cutoff, 'subclass');
data4_submini = mini_matrix(data4_sd_cutoff, 'subclass');
data5_submini = mini_matrix(data5_sd_cutoff, 'subclass');
data6_submini = mini_matrix(data6_sd_cutoff, 'subclass');
data7_submini = mini_matrix(data7_sd_cutoff, 'subclass');
data8_submini = mini_matrix(data8_sd_cutoff, 'subclass');
data_JSH_submini = mini_matrix(data_JSH_sd_cutoff, 'subclass');
data_N2U_submini = mini_matrix(data_N2U_sd_cutoff, 'subclass');
data_dauer_submini = mini_matrix(data_dauer_sd_cutoff, 'subclass');



%% syn (synapse expanded list) 

cd C:\Users\WORM-05\Documents\MATLAB\excel\

data1_syn = readmatrix('total_synapses_expanded_new.xlsx','Sheet','dataset1');
data2_syn = readmatrix('total_synapses_expanded_new.xlsx','Sheet','dataset2');
data3_syn = readmatrix('total_synapses_expanded_new.xlsx','Sheet','dataset3');
data4_syn = readmatrix('total_synapses_expanded_new.xlsx','Sheet','dataset4');
data5_syn = readmatrix('total_synapses_expanded_new.xlsx','Sheet','dataset5');
data6_syn = readmatrix('total_synapses_expanded_new.xlsx','Sheet','dataset6');
data7_syn = readmatrix('total_synapses_expanded_new.xlsx','Sheet','dataset7');
data8_syn = readmatrix('total_synapses_expanded_new.xlsx','Sheet','dataset8');
data_dauer_syn = readmatrix('total_synapses_expanded_new.xlsx','Sheet','dauer_221021');




%% ------------ Basic properties plotting ------------

%% No. of active zone

num_syn_1 = count_syn_cc(data1_syn);
num_syn_2 = count_syn_cc(data2_syn);
num_syn_3 = count_syn_cc(data3_syn);
num_syn_4 = count_syn_cc(data4_syn);
num_syn_5 = count_syn_cc(data5_syn);
num_syn_6 = count_syn_cc(data6_syn);
num_syn_7 = count_syn_cc(data7_syn);
num_syn_8 = count_syn_cc(data8_syn);
num_syn_dauer = count_syn_cc(data_dauer_syn);

x1 = [0 5 8 16 23 27 45 45];
y1 = [num_syn_1; num_syn_2; num_syn_3; num_syn_4; num_syn_5; num_syn_6; num_syn_7; num_syn_8];
y_lim1 = 4000;

figure(1)
subplot(1,3,[1,2]);
plot(x1, y1, '-o', 'color', 'black', 'LineWidth', 1, 'MarkerFaceColor','black')

xlabel('Hours after birth')
ylabel('No.of active zones')
xlim([0 50])
ylim([0 y_lim1])
set(gca,'xticklabel',[])
hold on;


plot([17 17],[0 y_lim1],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([24 24],[0 y_lim1],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([31 31],[0 y_lim1],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([38 38],[0 y_lim1],  'color', 'black', 'LineStyle' , ':')
hold on;



subplot(1,3,3);

x2 = 48;
y2 = num_syn_dauer;

plot(x2, y2, 'Marker', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'red')
xlim([24 72])
ylim([0 y_lim1])
set(gca,'xticklabel',[], 'yticklabel',[])





%% No. of total synapses

num_syn_expan_1 = count_syn_expan(data1_syn);
num_syn_expan_2 = count_syn_expan(data2_syn);
num_syn_expan_3 = count_syn_expan(data3_syn);
num_syn_expan_4 = count_syn_expan(data4_syn);
num_syn_expan_5 = count_syn_expan(data5_syn);
num_syn_expan_6 = count_syn_expan(data6_syn);
num_syn_expan_7 = count_syn_expan(data7_syn);
num_syn_expan_8 = count_syn_expan(data8_syn);
num_syn_expan_dauer = count_syn_expan(data_dauer_syn);


x3 = [0 5 8 16 23 27 45 45];
y3 = [num_syn_expan_1; num_syn_expan_2; num_syn_expan_3; num_syn_expan_4; num_syn_expan_5; num_syn_expan_6; num_syn_expan_7; num_syn_expan_8];
y_lim2 = 10000;

figure(2)
subplot(1,3,[1,2]);
plot(x3, y3, '-o', 'color', 'black', 'LineWidth', 1, 'MarkerFaceColor','black')

xlabel('Hours after birth')
ylabel('No.of synapses')
xlim([0 50])
ylim([0 y_lim2])
set(gca,'xticklabel',[])
hold on;


plot([17 17],[0 y_lim2],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([24 24],[0 y_lim2],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([31 31],[0 y_lim2],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([38 38],[0 y_lim2],  'color', 'black', 'LineStyle' , ':')
hold on;


subplot(1,3,3);

x4 = 48;
y4 = num_syn_expan_dauer;

plot(x4, num_syn_expan_dauer, 'Marker', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'red')
xlim([24 72])
ylim([0 y_lim2])
set(gca,'xticklabel',[], 'yticklabel',[])







%% No. of connections

num_conn_1 = num_connection(data1);
num_conn_2 = num_connection(data2);
num_conn_3 = num_connection(data3);
num_conn_4 = num_connection(data4);
num_conn_5 = num_connection(data5);
num_conn_6 = num_connection(data6);
num_conn_7 = num_connection(data7);
num_conn_8 = num_connection(data8);
num_conn_dauer = num_connection(data_dauer);

x5 = [0 5 8 16 23 27 45 45];
y5 = [num_conn_1; num_conn_2; num_conn_3; num_conn_4; num_conn_5; num_conn_6; num_conn_7; num_conn_8];
y_lim3 = 2500;

figure(3)
subplot(1,3,[1,2]);
plot(x5, y5, '-o', 'color', 'black', 'LineWidth', 1, 'MarkerFaceColor','black')

xlabel('Hours after birth')
ylabel('No.of connections')
xlim([0 50])
ylim([0 y_lim3])
set(gca,'xticklabel',[])
hold on;


plot([17 17],[0 y_lim3],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([24 24],[0 y_lim3],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([31 31],[0 y_lim3],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([38 38],[0 y_lim3],  'color', 'black', 'LineStyle' , ':')
hold on;


subplot(1,3,3);

x6 = 48;
y6 = num_conn_dauer;

plot(x6, y6, 'Marker', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'red')
xlim([24 72])
ylim([0 y_lim3])
set(gca,'xticklabel',[], 'yticklabel',[])




%% No. of connections - only neurons (include JSH, N2U) 

num_conn_neuron_1 = num_connection(data1_neuron);
num_conn_neuron_2 = num_connection(data2_neuron);
num_conn_neuron_3 = num_connection(data3_neuron);
num_conn_neuron_4 = num_connection(data4_neuron);
num_conn_neuron_5 = num_connection(data5_neuron);
num_conn_neuron_6 = num_connection(data6_neuron);
num_conn_neuron_7 = num_connection(data7_neuron);
num_conn_neuron_8 = num_connection(data8_neuron);
num_conn_neuron_JSH = num_connection(data_JSH_neuron);
num_conn_neuron_N2U = num_connection(data_N2U_neuron);
num_conn_neuron_dauer = num_connection(data_dauer_neuron);

x5 = [0 5 8 16 23 27 36 45 45 48];   
y5 = [num_conn_neuron_1; num_conn_neuron_2; num_conn_neuron_3; num_conn_neuron_4; num_conn_neuron_5; num_conn_neuron_6; num_conn_neuron_JSH; num_conn_neuron_7; num_conn_neuron_8; num_conn_neuron_N2U];
y_lim4 = 2500;

figure(4)
subplot(1,3,[1,2]);
plot(x5, y5, '-o', 'color', 'black', 'LineWidth', 1, 'MarkerFaceColor','black')

xlabel('Hours after birth')
ylabel('No.of connections')
xlim([0 50])
ylim([0 y_lim4])
set(gca,'xticklabel',[])
hold on;


plot([17 17],[0 y_lim4],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([24 24],[0 y_lim4],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([31 31],[0 y_lim4],  'color', 'black', 'LineStyle' , ':')
hold on;
plot([38 38],[0 y_lim4],  'color', 'black', 'LineStyle' , ':')
hold on;


subplot(1,3,3);

x6 = 48;
y6 = num_conn_dauer;

plot(x6, y6, 'Marker', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'red')
xlim([24 72])
ylim([0 y_lim4])
set(gca,'xticklabel',[], 'yticklabel',[])



%% nr. of synapses per active zone

data = data_dauer_syn(:,1);

tbl = tabulate(data);  
B = (tbl(:,3));
idx = find(B == 0);
B(idx) = [];
C = B*length(data)/100;

mean(C)



%% nr. synapses per connections

synperconn1 = num_syn_expan_1/num_conn_1
synperconn2 = num_syn_expan_2/num_conn_2
synperconn3 = num_syn_expan_3/num_conn_3
synperconn4 = num_syn_expan_4/num_conn_4
synperconn5 = num_syn_expan_5/num_conn_5
synperconn6 = num_syn_expan_6/num_conn_6
synperconn7 = num_syn_expan_7/num_conn_7
synperconn8 = num_syn_expan_8/num_conn_8
synperconn_dauer = num_syn_expan_dauer/num_conn_dauer




%% ------------ dataset_only connection ------------




% dauer only

ultra_common = data1_neuron+data2_neuron+data3_neuron+data4_neuron+data5_neuron+data6_neuron+data7_neuron+data8_neuron+data_JSH_neuron+data_N2U_neuron;

[only_1, only_2, common_du] = diff_mat(data_dauer_neuron,ultra_common);
data_dauer_neuron_only_list = adj2list(only_1);


% dataset1 only 

ultra_common = data2_neuron+data3_neuron+data4_neuron+data5_neuron+data6_neuron+data7_neuron+data8_neuron+data_JSH_neuron+data_N2U_neuron+data_dauer_neuron;

[only_1, only_2, common_du] = diff_mat(data1_neuron,ultra_common);
data1_neuron_only_list = adj2list(only_1);  % no dataset1 only

% dataset2 only

ultra_common = data1_neuron+data3_neuron+data4_neuron+data5_neuron+data6_neuron+data7_neuron+data8_neuron+data_JSH_neuron+data_N2U_neuron+data_dauer_neuron;

[only_1, only_2, common_du] = diff_mat(data2_neuron,ultra_common);
data2_neuron_only_list = adj2list(only_1);

% dataset3 only

ultra_common = data1_neuron+data2_neuron+data4_neuron+data5_neuron+data6_neuron+data7_neuron+data8_neuron+data_JSH_neuron+data_N2U_neuron+data_dauer_neuron;

[only_1, only_2, common_du] = diff_mat(data3_neuron,ultra_common);
data3_neuron_only_list = adj2list(only_1);

% dataset4 only

ultra_common = data1_neuron+data2_neuron+data3_neuron+data5_neuron+data6_neuron+data7_neuron+data8_neuron+data_JSH_neuron+data_N2U_neuron+data_dauer_neuron;

[only_1, only_2, common_du] = diff_mat(data4_neuron,ultra_common);
data4_neuron_only_list = adj2list(only_1);

% dataset5 only

ultra_common = data1_neuron+data2_neuron+data3_neuron+data4_neuron+data6_neuron+data7_neuron+data8_neuron+data_JSH_neuron+data_N2U_neuron+data_dauer_neuron;

[only_1, only_2, common_du] = diff_mat(data5_neuron,ultra_common);
data5_neuron_only_list = adj2list(only_1);

% dataset6 only

ultra_common = data1_neuron+data2_neuron+data3_neuron+data4_neuron+data5_neuron+data7_neuron+data8_neuron+data_JSH_neuron+data_N2U_neuron+data_dauer_neuron;

[only_1, only_2, common_du] = diff_mat(data6_neuron,ultra_common);
data6_neuron_only_list = adj2list(only_1);


% dataset7 only

ultra_common = data1_neuron+data2_neuron+data3_neuron+data4_neuron+data5_neuron+data6_neuron+data8_neuron+data_JSH_neuron+data_N2U_neuron+data_dauer_neuron;

[only_1, ~, common_du] = diff_mat(data7_neuron,ultra_common);
data7_neuron_only_list = adj2list(only_1);

% dataset8 only

ultra_common = data1_neuron+data2_neuron+data3_neuron+data4_neuron+data5_neuron+data6_neuron+data7_neuron+data_JSH_neuron+data_N2U_neuron+data_dauer_neuron;

[only_1, only_2, common_du] = diff_mat(data8_neuron,ultra_common);
data8_neuron_only_list = adj2list(only_1);

% JSH only

ultra_common = data1_neuron+data2_neuron+data3_neuron+data4_neuron+data5_neuron+data6_neuron+data7_neuron+data8_neuron+data_N2U_neuron+data_dauer_neuron;

[only_1, only_2, common_du] = diff_mat(data_JSH_neuron,ultra_common);
data_JSH_neuron_only_list = adj2list(only_1);

% N2U only

ultra_common = data1_neuron+data2_neuron+data3_neuron+data4_neuron+data5_neuron+data6_neuron+data7_neuron+data8_neuron+data_JSH_neuron+data_dauer_neuron;

[only_1, only_2, common_du] = diff_mat(data_N2U_neuron,ultra_common);
data_N2U_neuron_only_list = adj2list(only_1);



%% plot

conn_1 = sum(logical(data1_neuron), 'all');
conn_2 = sum(logical(data2_neuron), 'all');
conn_3 = sum(logical(data3_neuron), 'all');
conn_4 = sum(logical(data4_neuron), 'all');
conn_5 = sum(logical(data5_neuron), 'all');
conn_6 = sum(logical(data6_neuron), 'all');
conn_7 = sum(logical(data7_neuron), 'all');
conn_8 = sum(logical(data8_neuron), 'all');
conn_JSH = sum(logical(data_JSH_neuron), 'all');
conn_N2U = sum(logical(data_N2U_neuron), 'all');
conn_dauer = sum(logical(data_dauer_neuron), 'all');

y1 = length(data1_neuron_only_list)  / conn_1;
y2 = length(data2_neuron_only_list) / conn_2;
y3 = length(data3_neuron_only_list) / conn_3;
y4 = length(data4_neuron_only_list) / conn_4;
y5 = length(data5_neuron_only_list) / conn_5;
y6 = length(data6_neuron_only_list) / conn_6;
y7 = length(data7_neuron_only_list) / conn_7;
y8 = length(data8_neuron_only_list) / conn_8;
y_JSH = length(data_JSH_neuron_only_list) / conn_JSH;
y_N2U = length(data_N2U_neuron_only_list) / conn_N2U;
y_dauer = length(data_dauer_neuron_only_list) / conn_dauer;

x = [3 5 8 16 23 27 36 45 45 48];  % x = [3 5 8 16 23 27 36 45 45 48];  
y = [y1 y2 y3 y4 y5 y6 y_JSH y7 y8 y_N2U];



% Linear regression

x = x'
y = y'
X = [ones(length(x),1) x];
fit = X\y;
b = fit(1);  
a = fit(2); 

y_fitted = a*x + b;




figure(1)
subplot(1,3,[1,2]);
plot(x, y, '-o', 'color', 'none', 'LineWidth', 1, 'MarkerFaceColor','black')

xlim([0 50])
ylim([0 0.3])

hold on;
p2 = plot(x,y_fitted)

fitting = sprintf('y = %0.5fx %0.5f', a, b);
legend(p2, fitting, 'Location', 'northwest');


x2 = 47.5;
y2 = y_dauer;
subplot(1,3,3);
plot(x2, y2, 'Marker', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'red')
xlim([47 48])
ylim([0 0.3])

average = sum(y)/10;







%% stage-conserved connections

clear common_12;
[~, ~, common_12] = diff_mat(data1_neuron,data2_neuron);
[~, ~, common_12] = diff_mat(common_12,data3_neuron);
[~, ~, common_12] = diff_mat(common_12,data4_neuron);
[~, ~, common_12] = diff_mat(common_12,data5_neuron);
[~, ~, common_12] = diff_mat(common_12,data6_neuron);
[~, ~, common_12] = diff_mat(common_12,data7_neuron);
[~, ~, common_12] = diff_mat(common_12,data8_neuron);
[~, ~, common_12] = diff_mat(common_12,data_JSH_neuron);
[~, ~, common_12] = diff_mat(common_12,data_N2U_neuron);

[nondauer_only, ~, dauer_common_12] = diff_mat(common_12, data_dauer_neuron);
common_list = adj2list(dauer_common_12);
dauer_loss_list = adj2list(nondauer_only);

common_cell = id2name(common_list); 







%% ---------------------- L / R homolog plot -------------------

adj = data_dauer_sd_cutoff;



% make pair

dauer_list = adj2list_weight(adj);
dauer_list_3 = list2pair(dauer_list);  
[dauer_pair_left, solo_list] = onlypair(dauer_list_3);  
[paired_list, summary]=  pair2paired_list(dauer_pair_left); 

num_paired = summary(1)*2;    % nr. paired connection
num_solo = summary(2);        % nr. solo connection 
num_not_paired = length(dauer_list) - num_paired - num_solo;   


% weight(paired)
weight_list_pair = conn2weight_onlypair_and_solo(paired_list, 'dauer_221021_norm', 2); 


% weight(solo)
soloed_list = pair2paired_list(solo_list); 
soloed_list(:, 3:4) = [];


only_weight = [];
for ii = 1 : size(soloed_list)
    pre = soloed_list(ii,1);
    post = soloed_list(ii,2);
    only_weight = [only_weight; adj(pre,post)];
end

weight_list_solo = zeros(size(soloed_list,1), 2);
for ii = 1 : size(soloed_list,1)
    pre_code = solo_list(ii,4);
    if pre_code == 1;
        weight_list_solo(ii,1) = only_weight(ii);
        weight_list_solo(ii,2) = 0;
    else
        weight_list_solo(ii,1) = 0;
        weight_list_solo(ii,2) = only_weight(ii);
    end  
end

weight_list = [weight_list_pair; weight_list_solo]; 










% plotting


right = weight_list(:,2);
left = weight_list(:,1);

% Linear regression2
fit = fitlm(right,left);
fit_table = fit.Coefficients;
% y = ax + b
a = fit_table.Estimate(2);
b = fit_table.Estimate(1);


figure(1)
p1 = scatter(right,left, 10,'black','filled');
%p1= loglog(right, left, 'o','MarkerFaceColor','black', 'MarkerEdgeColor', 'None');
%title('Left versus right (weight)')
xlabel('Right homolog')
ylabel('Left homolog')
xlim([0 33])
ylim([0 33])

hold on;
p3 = plot([0 33], [0 33],':', 'Color', 'black')


hold on;
%plot([0 50000000000], [0 50000000000], 'Color', 'black');
x_fit = [0 33];
y_fit = [b 33*a+b];
p2 = plot(x_fit, y_fit, 'Color', 'red');  

fit_legend = sprintf('y = %0.2fx + %0.2f', a, b);
legend([p3 p2], 'y = x', fit_legend, 'Location', 'northwest', 'FontSize', 10); 








%% ------------------------ functions ------------------------


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



% list > weighted list

function weight_list = list2weight(edge_list, sheet_name)


adj = readmatrix('total_synapse_arranged.xlsx','Sheet', sheet_name,'Range','C3:IH242');
weight_list = zeros(size(edge_list, 1), 3);

for ii = 1 : size(edge_list, 1);
    
    pre = edge_list(ii,1);
    post = edge_list(ii,2);
    if ~isnan(pre)
        weight = adj(pre,post);
        weight_list(ii,:) = [pre post weight];
    end
    

end

end


% weighted list > rearrange 

function weight_arranged = weight_arrange(edge_list)

weight = edge_list(:,3);
[~, idx] = sort(weight, 'descend');
weight_arranged = edge_list(idx, :);


end




% cutoff matrix 

function [new_adj, sd] = cutoff_adj(adj, cutoff)

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





%% function - id transformation


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

idx = isnan(unique_code);
unique_code(idx) = [];

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








% connection type

function type_list = conn_type(list)

data_type = readtable('230724_neurontype_revised.csv');
subtype = data_type.Cook_subtypecode;

type_list = zeros(size(list,1), 3);

for ii = 1 : size(list,1)
    pre = list(ii, 1); 
    post = list(ii, 2);
    type_list(ii, 1:2) = [subtype(pre) subtype(post)];
    
    temp_pre = num2str(subtype(pre));
    temp_post = num2str(subtype(post));
    temp =  str2double([temp_pre temp_post]);
    type_list(ii, 3) = temp;

end

end





function type_list = conn_type_simple(list)

cd C:\Users\WORM-05\Documents\MATLAB\excel
data_type = readtable('230724_neurontype_revised.csv');
subtype = data_type.Cook_typecode2;

type_list = zeros(size(list,1), 3);

for ii = 1 : size(list,1)
    pre = list(ii, 1); 
    post = list(ii, 2);
    type_list(ii, 1:2) = [subtype(pre) subtype(post)];
    
    temp_pre = num2str(subtype(pre));
    temp_post = num2str(subtype(post));
    temp =  str2double([temp_pre temp_post]);
    type_list(ii, 3) = temp;

end

end







%% function - Basic properties plotting


% nr.synapse

function num_syn = count_syn_cc(list)

syn_id = list(:,1);
syn_id = unique(syn_id);
num_syn = length(syn_id);

end


% nr.synapse (expanded)

function num_syn = count_syn_expan(list)

num_syn = size(list, 1);

end


% nr.connection

function num_conn = num_connection(adj)   

adj = logical(adj);
adj_tf = double(adj); 

num_conn = sum(adj_tf, 'all');


end





% Comparing connections

function [only_1, only_2, common_12] = diff_mat(adj1,adj2)

num_node = size(adj1, 1);

adj1 = logical(adj1);
adj2 = logical(adj2);

only_1 = zeros(num_node);
only_2 = zeros(num_node);
common_12 = zeros(num_node);

for ii = 1 : num_node^2
    
    if adj1(ii) - adj2(ii) > 0
       only_1(ii) = adj1(ii);
    elseif adj1(ii) - adj2(ii) < 0
       only_2(ii) = adj2(ii);
    elseif adj1(ii)+adj2(ii) == 2
        common_12(ii) = adj1(ii);
    end

end

conn_1 = sum(logical(adj1), 'all');  % adj1 connection nr
conn_only_1 = sum(only_1, 'all');    % adj1 only connection nr
conn_2 = sum(logical(adj2), 'all');  % adj2 connection nr
conn_only_2 = sum(only_2, 'all');    % adj2 only connection nr
comm_12 = sum(common_12, 'all');     % adj1, adj2 common connection nr

fprintf('adj1 : total connection-%d, only-%d\nadj2 : total connection-%d, only-%d\ncommon connection : %d\n', conn_1, conn_only_1, conn_2, conn_only_2, comm_12)

end







%% pairing functions



function new_list = list2pair(test_list)


cd C:\Users\WORM-05\Documents\MATLAB\excel\
key = readtable('segmentation_key_common_230926.csv');

common_ID = key.Common_ID;
pair_code = key.pair_code;
LR_code = key.LR_code;

key_list = [common_ID pair_code LR_code];




% test list (n * 2) > pair_code, LR_code, pair_code, LR_code (n * 4)


new_list = zeros(size(test_list, 1), 4);


for ii = 1 : size(test_list, 1)
    
    pre = test_list(ii,1);
    post = test_list(ii,2);
    
    idx_pre = pre==common_ID;
    pair_pre = pair_code(idx_pre);
    LR_pre = LR_code(idx_pre);
    
    idx_post = post==common_ID;
    pair_post = pair_code(idx_post);
    LR_post = LR_code(idx_post);    
    
    new_list(ii,:) = [pair_pre LR_pre pair_post LR_post];
    
end


end






% new_list(n*4) > only 'pair' left

function [paired_list, solo_list] = onlypair(new_list)


paired_list = [];
solo_list = [];

for ii = 1 : size(new_list, 1)
    
    temp = new_list(ii, :);      
    
    
    if temp(2) == 1
        if temp(4) == 1            
            target = [temp(1) 2 temp(3) 2];            

        elseif temp(4) == 2
            target = [temp(1) 2 temp(3) 1];            
            
        else % temp(4) == 3
            target = [temp(1) 2 temp(3) 3];    
            
        end

        
    elseif temp(2) == 2
        if temp(4) == 1
            target = [temp(1) 1 temp(3) 2];                
            
        elseif temp(4) == 2
            target = [temp(1) 1 temp(3) 1];    
            
        else % temp(4) == 3
            target = [temp(1) 1 temp(3) 3];    
            
        end
        
    else % temp(2) == 3
        if temp(4) == 1
            target = [temp(1) 3 temp(3) 2];    
            
        elseif temp(4) == 2
            target = [temp(1) 3 temp(3) 1];    
            
        else  % 3-3 case
            target = [temp(1) 3 temp(3) 3];

        end
        
    end
    
    idx = ismember(new_list, target,'rows');  
    if sum(idx)
        paired_list = [paired_list; temp];
    else
        solo_list = [solo_list; temp];
    end
    
    
    
end




end






% paired_list (n *4) > paired edge_list, how much pair, how much solo

function [new_edge_list, summary] =  pair2paired_list(paired_list)



cd C:\Users\WORM-05\Documents\MATLAB\excel\
key = readtable('segmentation_key_common_230926.csv');

common_ID = key.Common_ID;
pair_code = key.pair_code;
LR_code = key.LR_code;
pair_LR = [pair_code LR_code];





new_edge_list = [];
copy_paired_list = paired_list;
solo_count = 0;

for ii = 1 : size(paired_list, 1)
    
    temp = paired_list(ii, :);
    pre = temp(1:2);
    pre1 = temp(1);
    pre2 = temp(2);
    post = temp(3:4);
    post1 = temp(3);
    post2 = temp(4);
    
    idx = ismember(copy_paired_list, temp,'rows');
    
    if sum(idx)

        if pre2 == 1
            if post2 == 1
                pair_pre = [pre1 2];
                pair_post = [post1 2];
            elseif post2 == 2
                pair_pre = [pre1 2];
                pair_post = [post1 1];
            else
                pair_pre = [pre1 2];
                pair_post = [post1 3];                
            end
                    
        elseif pre2 == 2
            if post2 == 1
                pair_pre = [pre1 1];
                pair_post = [post1 2];
            elseif post2 == 2
                pair_pre = [pre1 1];
                pair_post = [post1 1];
            else
                pair_pre = [pre1 1];
                pair_post = [post1 3];         
            end
            
        else  % pre2 == 3
            if post2 == 1
                pair_pre = [pre1 3];
                pair_post = [post1 2];                
            elseif post2 == 2
                pair_pre = [pre1 3];
                pair_post = [post1 1];                   
            else  
                pair_pre = [0 0];
                pair_post = [0 0];
            end
        end
    

        check = [pair_pre pair_post];
        check_idx = ismember(copy_paired_list, check,'rows');

    
    
        if sum(idx + check_idx)
            copy_paired_list(logical(idx + check_idx), :) = [];
        end

    
    
        % common_id
    
        idx_pre = ismember(pair_LR, pre,'rows');    
        idx_post = ismember(pair_LR, post,'rows');
        pre_comm = common_ID(idx_pre);
        post_comm = common_ID(idx_post);

        idx_check_pre = ismember(pair_LR, pair_pre,'rows');    
        idx_check_post = ismember(pair_LR, pair_post,'rows');
        pre_check_comm = common_ID(idx_check_pre);
        post_check_comm = common_ID(idx_check_post);   
    
        pair_row = [pre_comm post_comm pre_check_comm post_check_comm];
    
        if length(pair_row) == 2
            pair_row = [pair_row NaN NaN];
            solo_count = solo_count + 1;
        end
    
        new_edge_list = [new_edge_list; pair_row];
        
        
        
    else  

    end    
  
    
    
        
end


pair_count = size(new_edge_list,1) - solo_count;
summary = [pair_count, solo_count];


end   






function weight_list = conn2weight_onlypair(edge_list, sheet_name) 



cd C:\Users\WORM-05\Documents\MATLAB\excel\

adj = readmatrix('total_connection_arranged_new.xlsx','Sheet', sheet_name ,'Range','C3:IH242');






weight_list = zeros(size(edge_list,1), size(edge_list,2)/2 );
solo_list = [];

for ii = 1 : size(edge_list, 1)
    
    temp = edge_list(ii, :);
    
    if length(temp) == 2
        
        pre = temp(1);
        post = temp(2);   
        weight = adj(pre,post);
        weight_list(ii) = weight;
           
    else % length(temp) >= 4
        
        for jj = 1 : size(edge_list,2)/2
            pre = temp(jj*2-1);
            post = temp(jj*2);   
            if ~isnan(pre)
                weight = adj(pre,post);
                weight_list(ii, jj) = weight;   
            else
                weight = 0;
                weight_list(ii, jj) = weight;   
                solo_list = [solo_list; ii];
            end
        end
    end
    
    
end

weight_list(solo_list, :) = [];


end






function weight_list = conn2weight_onlypair_and_solo(edge_list, sheet_name, code) 


cd C:\Users\WORM-05\Documents\MATLAB\excel\
adj = readmatrix('total_connection_arranged_new.xlsx','Sheet', sheet_name ,'Range','C3:IH242');

if code == 1
    list = adj2list(adj(1:222,1:222));
elseif code == 2
    list = adj2list(adj(1:180, 1:180));
else
print('all adj = 1, only neuron = 2' )
end


% solo pair 

weight_list = [];

for ii = 1 : size(edge_list, 1)
    
    temp = edge_list(ii, :);

        for jj = 1 : 2
            pre = temp(jj*2-1);
            post = temp(jj*2);   
            if ~isnan(pre)
                weight = adj(pre,post);
                weight_list(ii, jj) = weight;   
            else
                weight = weight_list(ii, 1);
                weight_list(ii, jj) = weight;   
            end
        end
    
end



end




