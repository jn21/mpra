mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot_normalized.txt','Delimiter','\t');

%Remove non-finite entries
idx = isfinite(mpra_data{:,'enhancer_activity'}) & isfinite(mpra_data{:,'promoter_activity'});
mpra_data = mpra_data(idx,:);

%Remove noisy points
% THRESHOLD = 2^7;
% idx = (mpra_data{:,'Rep1_ATotal'} > THRESHOLD & ...
%     mpra_data{:,'Rep2_ATotal'} > THRESHOLD & ...
%     mpra_data{:,'Rep1_BTotal'} > THRESHOLD & ...
%     mpra_data{:,'Rep2_BTotal'} > THRESHOLD);
% 
% mpra_data = mpra_data(idx,:);

mpra_data_intact = subset_table(mpra_data,'is_intact_sequence',1);

figdir = '~/Documents/mpra/fig/technical_quality';

%% Distributions of input DNA
figure
data1 = log2(mpra_data{:,'DNAInput_ATotal'});
data2 = log2(mpra_data{:,'DNAInput_BTotal'});
scatter(data1,data2);
xlabel('Input "A" DNA Barcodes  (log_2)')
ylabel('Input "B" DNA Barcodes  (log_2)')
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

%saveas(gcf,fullfile(figdir,'input_dna_scatter'),'png')

%% Replicate Reproducibility - A Barcode Counts
figure
data1 = log2(mpra_data{:,'Rep1_ATotal'});
data2 = log2(mpra_data{:,'Rep2_ATotal'});
scatter(data1,data2);
xlabel('"A" RNA Barcode Counts - Replicate 1  (log_2)')
ylabel('"A" RNA Barcode Counts - Replicate 2  (log_2)')
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

%saveas(gcf,fullfile(figdir,'A_RNA_barcode_replicates'),'png')

%% Replicate Reproducibility - B Barcode Counts 

figure
data1 = log2(mpra_data{:,'Rep1_BTotal'});
data2 = log2(mpra_data{:,'Rep2_BTotal'});
scatter(data1,data2);
xlabel('"B" RNA Barcode Counts - Replicate 1  (log_2)')
ylabel('"B" RNA Barcode Counts - Replicate 2  (log_2)')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');
ax.XLim = [0 14];
ax.YLim = [0 14];

input_corr = corr(data1,data2);
title_str = sprintf('"B" Barcodes \n Pearson Correlation = %.3f',input_corr);
title(title_str)

%saveas(gcf,fullfile(figdir,'B_RNA_barcode_replicates'),'png')

%% Replicate Reproducibility - A RNA Ratio - Promoter Activity
figure
data1 = mpra_data{:,'rep1_A_ratio'};
data2 = mpra_data{:,'rep2_A_ratio'};
scatter(data1,data2);
xlabel('"A" Barcode Ratio - Replicate 1')
ylabel('"A" Barcode Ratio - Replicate 2')
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

%saveas(gcf,fullfile(figdir,'A_ratio_promoter_activity_replicates'),'png')

%% Replicate Reproducibility - B RNA Ratio - Enhancer Activity
figure
data1 = mpra_data{:,'rep1_B_ratio'};
data2 = mpra_data{:,'rep2_B_ratio'};
scatter(data1,data2);
xlabel('"B" Barcode Ratio - Replicate 1')
ylabel('"B" Barcode Ratio - Replicate 2')
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

%saveas(gcf,fullfile(figdir,'B_ratio_enhancer_activity_replicates'),'png')

%% Replicate Reproducibility - A RNA Ratio - Promoter Activity - Intact sequences
figure
data1 = mpra_data_intact{:,'rep1_A_ratio'};
data2 = mpra_data_intact{:,'rep2_A_ratio'};
scatter(data1,data2);
xlabel('"A" Barcode Ratio - Replicate 1')
ylabel('"A" Barcode Ratio - Replicate 2')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
min_lim = min([ax.YLim ax.XLim]);
hold on
plot([min_lim max_lim],[min_lim max_lim],'r:');
ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];

input_corr = corr(data1,data2);
title_str = sprintf('"A" Barcode Ratio (Promoter Activity) \n Intact Sequences \n Pearson Correlation = %.3f',input_corr);
title(title_str)

%saveas(gcf,fullfile(figdir,'A_ratio_promoter_activity_replicates_intact'),'png')

%% Replicate Reproducibility - B RNA Ratio - Enhancer Activity - Intact sequences
figure
data1 = mpra_data_intact{:,'rep1_B_ratio'};
data2 = mpra_data_intact{:,'rep2_B_ratio'};
scatter(data1,data2);
xlabel('"B" Barcode Ratio - Replicate 1')
ylabel('"B" Barcode Ratio - Replicate 2')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
min_lim = min([ax.YLim ax.XLim]);
hold on
plot([min_lim max_lim],[min_lim max_lim],'r:');
ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];

input_corr = corr(data1,data2);
title_str = sprintf('"B" Barcode Ratio (Enhancer Activity) \n Intact Sequences \n Pearson Correlation = %.3f',input_corr);
title(title_str)

%saveas(gcf,fullfile(figdir,'B_ratio_enhancer_activity_replicates_intact'),'png')
