% fid = fopen('~/Documents/mpra/data/mpra_data.txt');
% 
% A = textscan(fid)

mpra_data = readtable('~/Documents/mpra/data/mpra_data.txt','Delimiter','\t');
construct_annot = readtable('~/Documents/mpra/annot/mpra_data_and_construct_annot.txt','Delimiter','\t');
combined_data = readtable('/Users/jnasser/Documents/mpra/annot/mpra_data_and_construct_annot.txt','Delimiter','\t');
%% Distributions of input DNA
figure
scatter(mpra_data{:,'DNAInput_ETotal'},mpra_data{:,'DNAInput_PTotal'})
xlabel('DNAInput_ETotal','Interpreter','None')
ylabel('DNAInput_PTotal','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');

figure
scatter(log2(mpra_data{:,'DNAInput_ETotal'}),log2(mpra_data{:,'DNAInput_PTotal'}))
xlabel('log2(DNAInput_ETotal)','Interpreter','None')
ylabel('log2(DNAInput_PTotal)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');

%% Replicate Correlation
figure
scatter(mpra_data{:,'Rep1_ERatio'},mpra_data{:,'Rep2_ERatio'})
figure
scatter(mpra_data{:,'Rep1_PRatio'},mpra_data{:,'Rep2_PRatio'})

%% E and P Activity
figure
scatter(mpra_data{:,'Rep1_ERatio'},mpra_data{:,'Rep1_PRatio'})
figure
scatter(mpra_data{:,'Rep2_ERatio'},mpra_data{:,'Rep2_PRatio'})
