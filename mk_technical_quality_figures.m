mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');

%Remove non-finite entries
idx = isfinite(mpra_data{:,'P_ratio_avg_rep'}) & isfinite(mpra_data{:,'E_ratio_avg_rep'});
mpra_data = mpra_data(idx,:);

figdir = '~/Documents/mpra/fig/technical_quality';

%% Distributions of input DNA
figure
data1 = log2(mpra_data{:,'DNAInput_ETotal'});
data2 = log2(mpra_data{:,'DNAInput_PTotal'});
scatter(data1,data2);
xlabel('log2(DNAInput_ETotal)','Interpreter','None')
ylabel('log2(DNAInput_PTotal)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');
ax.XLim = [0 14];
ax.YLim = [0 14];

input_corr = corr(data1,data2);
title_str = sprintf('Input DNA \n Pearson Correlation = %.2f',input_corr);
title(title_str)

saveas(gcf,fullfile(figdir,'input_scatter'),'png')

%% Replicate Reproducibility - Barcode Counts - Construct Barcodes
figure
data1 = log2(mpra_data{:,'Rep1_ETotal'});
data2 = log2(mpra_data{:,'Rep2_ETotal'});
scatter(data1,data2);
xlabel('log2(Rep1 Barcode Counts)','Interpreter','None')
ylabel('log2(Rep2 Barcode Counts)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');
ax.XLim = [0 14];
ax.YLim = [0 14];

input_corr = corr(data1,data2);
title_str = sprintf('Construct Barcodes \n Pearson Correlation = %.2f',input_corr);
title(title_str)

saveas(gcf,fullfile(figdir,'construct_barcode_replicates'),'png')

%% Replicate Reproducibility - Barcode Counts - GFP Barcodes

figure
data1 = log2(mpra_data{:,'Rep1_PTotal'});
data2 = log2(mpra_data{:,'Rep2_PTotal'});
scatter(data1,data2);
xlabel('log2(Rep1 Barcode Counts)','Interpreter','None')
ylabel('log2(Rep2 Barcode Counts)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');
ax.XLim = [0 14];
ax.YLim = [0 14];

input_corr = corr(data1,data2);
title_str = sprintf('GFP Barcodes \n Pearson Correlation = %.2f',input_corr);
title(title_str)

saveas(gcf,fullfile(figdir,'gfp_barcode_replicates'),'png')

%% Replicate Reproducibility - Promoter Activity
figure
data1 = mpra_data{:,'rep1_E_ratio'};
data2 = mpra_data{:,'rep2_E_ratio'};
scatter(data1,data2);
xlabel('log2(Rep1 Construct Ratio)','Interpreter','None')
ylabel('log2(Rep2 Construct Ratio)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
min_lim = min([ax.YLim ax.XLim]);
hold on
plot([min_lim max_lim],[min_lim max_lim],'r:');
ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];

input_corr = corr(data1,data2);
title_str = sprintf('Construct Barcode Ratio (Promoter Activity) \n Pearson Correlation = %.2f',input_corr);
title(title_str)

saveas(gcf,fullfile(figdir,'promoter_activity_replicates'),'png')

%% Replicate Reproducibility - Enhancer Activity
figure
data1 = mpra_data{:,'rep1_P_ratio'};
data2 = mpra_data{:,'rep2_P_ratio'};
scatter(data1,data2);
xlabel('log2(Rep1 GFP Ratio)','Interpreter','None')
ylabel('log2(Rep2 GFP Ratio)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
min_lim = min([ax.YLim ax.XLim]);
hold on
plot([min_lim max_lim],[min_lim max_lim],'r:');
ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];

input_corr = corr(data1,data2);
title_str = sprintf('GFP Barcodes Ratio (Enhancer Activity) \n Pearson Correlation = %.2f',input_corr);
title(title_str)

saveas(gcf,fullfile(figdir,'enhancer_activity_replicates'),'png')
