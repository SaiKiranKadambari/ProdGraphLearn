clc;clear all;close all

% This code is for generaring the Figure 5 in the paper. 



load('ResultsAQI_12months.mat') 
% rng('default');
% Plot India map
I = imread('IndiaMap_2.jpg');
figure(1);clf;imagesc(I);
axis off
hold on;

% Plot the AQI data Space graph

load('PixelCoords.mat')
coords = NameCoords;
C = coords;

Wp_i(abs(Wp_i)<0.05) = 0;
Gp = gsp_graph(Wp_i);
Gp.coords = coords;
figure(1);
%gsp_plot_graph(Gp)
%gsp_plot_signal(Gp, X_data(:,24))
[Xg,Yg] = gplot(Gp.L, coords);
plot(Xg,Yg,'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',8,...
    'Marker','o',...
    'Color',[0.494117647409439 0.494117647409439 0.494117647409439]);


names{8} = '';
names{12} = '';
names{24} = '';
names{13} = '';
names{22} = '';
names{16} = '';
names{7}= '';
names{8} = '';
names{28} = '';
% names{23} = '';
names{20} = '';
names{1} = '';
names{2} = '';
names{6} = '';
names{19} = '';
names{5} = '';
names{4} = '';
names{14} = '';
Txtcoords = coords;
Txtcoords(17,2) = coords(17,2) + 20;  % offset for overlapping text labels
Txtcoords(21,1) = coords(21,1) + 20;  % offset for overlapping text labels
Txtcoords(3,1) = coords(3,1) + 10;  % offset for overlapping text labels
Txtcoords(30,1) = coords(30,1) + 10;  % offset for overlapping text labels
Txtcoords(27,1) = coords(27,1) + 15;  % offset for overlapping text labels
text(Txtcoords(:,1), Txtcoords(:,2), names, 'VerticalAlignment','bottom',...
    'HorizontalAlignment','center','Position',[0.1,0.1,0], 'FontSize',11.5)

set(gca,'ydir','reverse')


print('-depsc', '.\Figures\AQI_Data_SpaceGraph_up.eps');



G2 = graph(Wq_i);
figure(2);fq = plot(G2,'Layout','force')

% % Months = {'Jan-18','Feb-18','Mar-18','Apr-18','May-18','Jun-18','Jul-18','Aug-18'...
%     ,'Sept-18','Oct-18','Nov-18','Dec-18','Jan-19','Feb-19','Mar-19','Apr-19','May-19'...
%     ,'Jun-19','July','Aug-19','Sept-19','Oct-19','Nov-19','Dec-19'};
Months = {'January','February','March','April','May','June','July','August','September','October','November','December'};
fq.NodeFontSize = 12;
fq.NodeLabel = Months;
fq.MarkerSize = 7;
% G2.Nodes.value = X_data(8,:)';
% G2.Nodes.NodeColors = G2.Nodes.value;
% fq.NodeCData = G2.Nodes.NodeColors;
set(gca,'Visible','off')

print('-depsc', '.\Figures\AQI_Data_TimeGraph_up.eps');



G2 = graph(Wq_i);
figure(3);fq = plot(G2,'Layout','force')

% % Months = {'Jan-18','Feb-18','Mar-18','Apr-18','May-18','Jun-18','Jul-18','Aug-18'...
%     ,'Sept-18','Oct-18','Nov-18','Dec-18','Jan-19','Feb-19','Mar-19','Apr-19','May-19'...
%     ,'Jun-19','July','Aug-19','Sept-19','Oct-19','Nov-19','Dec-19'};
Months = {'January','February','March','April','May','June','July','August','September','October','November','December'};
fq.NodeFontSize = 12;
fq.NodeLabel = Months;
fq.MarkerSize = 7;
G2.Nodes.value = X_data(8,:)';
G2.Nodes.NodeColors = G2.Nodes.value;
fq.NodeCData = G2.Nodes.NodeColors;
colorbar;
set(gca,'Visible','off')
print('-depsc', '.\Figures\AQI_Data_TimeGraph_up_withNodevalues.eps');