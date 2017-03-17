% fid = fopen('~/Documents/mpra/data/mpra_data.txt');
% 
% A = textscan(fid)

A = readtable('~/Documents/mpra/data/mpra_data.txt','Delimiter','\t');

%% Distributions of input DNA
figure
scatter(A{:,'DNAInput_ETotal'},A{:,'DNAInput_PTotal'})
xlabel('DNAInput_ETotal','Interpreter','None')
ylabel('DNAInput_PTotal','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');

figure
scatter(log2(A{:,'DNAInput_ETotal'}),log2(A{:,'DNAInput_PTotal'}))
xlabel('log2(DNAInput_ETotal)','Interpreter','None')
ylabel('log2(DNAInput_PTotal)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');

%% Replicate Correlation
figure
scatter(A{:,'Rep1_ERatio'},A{:,'Rep2_ERatio'})
figure
scatter(A{:,'Rep1_PRatio'},A{:,'Rep2_PRatio'})

%% E and P Activity
figure
scatter(A{:,'Rep1_ERatio'},A{:,'Rep1_PRatio'})
figure
scatter(A{:,'Rep2_ERatio'},A{:,'Rep2_PRatio'})
