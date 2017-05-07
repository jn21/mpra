mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');

%Remove non-finite entries
idx = isfinite(mpra_data{:,'enhancer_activity'}) & isfinite(mpra_data{:,'promoter_activity'});
mpra_data = mpra_data(idx,:);

figdir = '~/Documents/mpra/fig/technical_quality';

%% Distributions of input DNA
figure
data1 = log2(mpra_data{:,'DNAInput_ATotal'});
data2 = log2(mpra_data{:,'DNAInput_BTotal'});
scatter(data1,data2);
xlabel('log2(Input "A" DNA Barcodes)','Interpreter','None')
ylabel('log2(Input "B" DNA Barcodes)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');
ax.XLim = [0 14];
ax.YLim = [0 14];

input_corr = corr(data1,data2);
title_str = sprintf('Input DNA \n Pearson Correlation = %.3f',input_corr);
title(title_str)

%saveas(gcf,fullfile(figdir,'input_scatter'),'png')

%% Replicate Reproducibility - A Barcode Counts
figure
data1 = log2(mpra_data{:,'Rep1_ATotal'});
data2 = log2(mpra_data{:,'Rep2_ATotal'});
scatter(data1,data2);
xlabel('log2("A" RNA Barcode Counts - Replicate 1)','Interpreter','None')
ylabel('log2("A" RNA Barcode Counts - Replicate 2)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');
ax.XLim = [0 14];
ax.YLim = [0 14];

input_corr = corr(data1,data2);
title_str = sprintf('"A" Barcodes \n Pearson Correlation = %.3f',input_corr);
title(title_str)

%saveas(gcf,fullfile(figdir,'construct_barcode_replicates'),'png')

%% Replicate Reproducibility - Barcode Counts - GFP Barcodes

figure
data1 = log2(mpra_data{:,'Rep1_BTotal'});
data2 = log2(mpra_data{:,'Rep2_BTotal'});
scatter(data1,data2);
xlabel('log2("B" RNA Barcode Counts - Replicate 1)','Interpreter','None')
ylabel('log2("B" RNA Barcode Counts - Replicate 2)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');
ax.XLim = [0 14];
ax.YLim = [0 14];

input_corr = corr(data1,data2);
title_str = sprintf('GFP Barcodes \n Pearson Correlation = %.3f',input_corr);
title(title_str)

%saveas(gcf,fullfile(figdir,'gfp_barcode_replicates'),'png')

%% Replicate Reproducibility - Promoter Activity
figure
data1 = mpra_data{:,'rep1_A_ratio'};
data2 = mpra_data{:,'rep2_A_ratio'};
scatter(data1,data2);
xlabel('"A" Barcode Ratio - Replicate 1','Interpreter','None')
ylabel('"A" Barcode Ratio - Replicate 2','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
min_lim = min([ax.YLim ax.XLim]);
hold on
plot([min_lim max_lim],[min_lim max_lim],'r:');
ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];

input_corr = corr(data1,data2);
title_str = sprintf('"A" Barcode Ratio (Promoter Activity) \n Pearson Correlation = %.3f',input_corr);
title(title_str)

%saveas(gcf,fullfile(figdir,'promoter_activity_replicates'),'png')

%% Replicate Reproducibility - Enhancer Activity
figure
data1 = mpra_data{:,'rep1_B_ratio'};
data2 = mpra_data{:,'rep2_B_ratio'};
scatter(data1,data2);
xlabel('"B" Barcode Ratio - Replicate 1','Interpreter','None')
ylabel('"B" Barcode Ratio - Replicate 2','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
min_lim = min([ax.YLim ax.XLim]);
hold on
plot([min_lim max_lim],[min_lim max_lim],'r:');
ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];

input_corr = corr(data1,data2);
title_str = sprintf('"B" Barcode Ratio (Enhancer Activity) \n Pearson Correlation = %.3f',input_corr);
title(title_str)

%saveas(gcf,fullfile(figdir,'enhancer_activity_replicates'),'png')
