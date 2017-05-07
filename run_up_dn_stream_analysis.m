mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');

%% Run ID analysis for promoter and enhancer activity
save_fig = false;
up_gfp = examine_activity_by_up_or_dn_id('up','enhancer',save_fig);
dn_gfp = examine_activity_by_up_or_dn_id('dn','gfp',save_fig);
up_construct = examine_activity_by_up_or_dn_id('up','construct',save_fig);
dn_construct = examine_activity_by_up_or_dn_id('dn','construct',save_fig);

%% Scatter E and P activity of all constructs
bad_idx = mpra_data{:,'dnstream_is_reverse'} |...
    mpra_data{:,'upstream_is_reverse'} |...
    mpra_data{:,'dnstream_is_modified'} |...
    ~isfinite(mpra_data{:,'enhancer_activity'}) |...
    ~isfinite(mpra_data{:,'promoter_activity'});
good_mpra_data = mpra_data(~bad_idx,:);

scatter(good_mpra_data{:,'enhancer_activity'},good_mpra_data{:,'promoter_activity'});
scatter_corr = corr(good_mpra_data{:,'enhancer_activity'},good_mpra_data{:,'promoter_activity'});

xlabel('Enhancer Activity')
ylabel('Promoter Activity')
title_str = sprintf('All non-Modified Constructs \n correlation = %.2f',scatter_corr);
title(title_str)
%saveas(gcf,'~/Documents/mpra/fig/up_dn_analysis/global_ep_scatter','png')

%% Scatter Enhancer and Promoter Activity of Median UP and DN stream sequences
up_gfp_no_rev = subset_table(up_gfp,'is_reverse',0);
dn_gfp_no_rev = subset_table(dn_gfp,'is_reverse',0);
up_construct_no_rev = subset_table(up_construct,'is_reverse',0);
dn_construct_no_rev = subset_table(dn_construct,'is_reverse',0);

figure
scatter(up_gfp_no_rev{:,'median_ratio'},up_construct_no_rev{:,'median_ratio'})
scatter_corr = corr(up_gfp_no_rev{:,'median_ratio'},up_construct_no_rev{:,'median_ratio'});
% hold on
% ax = gca;
% plot([up_gfp{1,'neg_control_median_ratio'} up_gfp{1,'neg_control_median_ratio'}], ax.YLim, 'r:')
% hold on
% plot(ax.XLim, [up_construct{1,'neg_control_median_ratio'} up_construct{1,'neg_control_median_ratio'}], 'r:')
xlabel('Median Enhancer Activity')
ylabel('Median Promoter Activity')
title_str = sprintf('Upstream Sequences - Medians \n correlation = %.2f',scatter_corr);
title(title_str)
%saveas(gcf,'~/Documents/mpra/fig/up_dn_analysis/up_ep_median_scatter','png')

figure
scatter(dn_gfp_no_rev{:,'median_ratio'},dn_construct_no_rev{:,'median_ratio'})
scatter_corr = corr(dn_gfp_no_rev{:,'median_ratio'},dn_construct_no_rev{:,'median_ratio'});
% hold on
% ax = gca;
% plot([dn_gfp{1,'neg_control_median_ratio'} dn_gfp{1,'neg_control_median_ratio'}], ax.YLim, 'r:')
% hold on
% plot(ax.XLim, [dn_construct{1,'neg_control_median_ratio'} dn_construct{1,'neg_control_median_ratio'}], 'r:')
xlabel('Median Enhancer Activity')
ylabel('Median Promoter Activity')
title_str = sprintf('Downstream Sequences - Medians \n correlation = %.2f',scatter_corr);
title(title_str)
%saveas(gcf,'~/Documents/mpra/fig/up_dn_analysis/dn_ep_median_scatter','png')

%% FDR 
fdr_cutoff = .1;
nnz(up_gfp{:,'qval'} < fdr_cutoff)
nnz(dn_gfp{:,'qval'} < fdr_cutoff)
nnz(up_construct{:,'qval'} < fdr_cutoff)
nnz(dn_construct{:,'qval'} < fdr_cutoff)

%%
% figure
% gscatter([up_gfp{:,'median_ratio'}; dn_gfp{:,'median_ratio'}],...
%     [up_construct{:,'median_ratio'}; dn_construct{:,'median_ratio'}],...
%     [repmat('upstream',height(up_gfp),1); repmat('dnstream',height(dn_gfp),1)])
% hold on
% % ax = gca;
% % plot([up_gfp{1,'neg_control_median_ratio'} up_gfp{1,'neg_control_median_ratio'}], ax.YLim, 'r:')
% % hold on
% % plot(ax.XLim, [up_construct{1,'neg_control_median_ratio'} up_construct{1,'neg_control_median_ratio'}], 'r:')
% xlabel('Median Enhancer Activity')
% ylabel('Median Promoter Activity')
% %title('Upstream Sequences')

%% Promoter Activty / Enhancer Activity
data = up_construct{:,'median_ratio'} ./ up_gfp{:,'median_ratio'};
boxplot(data,up_construct{:,'median_ratio'});

%% Enhancer Activity by region ID
up_gfp_no_reverse = subset_table(up_gfp,'is_reverse',0);
dn_gfp_no_reverse = subset_table(dn_gfp,'is_reverse',0);
common_regions = intersect([up_gfp_no_reverse.region_id],[dn_gfp_no_reverse.region_id]);
up_gfp_no_reverse = subset_table(up_gfp,'region_id',common_regions);
dn_gfp_no_reverse = subset_table(dn_gfp,'region_id',common_regions);
figure
assert(isequal(up_gfp_no_reverse{:,'region_id'},dn_gfp_no_reverse{:,'region_id'}),'regions not in same order')
gscatter(up_gfp_no_reverse{:, 'median_ratio'}, dn_gfp_no_reverse{:, 'median_ratio'},up_gfp_no_reverse{:, 'relation_to_tss'})
xlabel('Upstream Sequence')
ylabel('Downstream Sequence')
title('Enhancer Activity (GFP Barcode Ratio)')
ax = gca;
max_lim = max([ax.XLim ax.YLim]);
min_lim = min([ax.XLim ax.YLim]);

ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];
hold on
plot([min_lim max_lim],[min_lim max_lim],'r:')
l = legend;
set(l, 'Interpreter', 'None')
grid on
saveas(gcf,'~/Documents/mpra/fig/up_dn_analysis/enhancer_activity_scatter_by_region','png')

%% Promoter Activity by region ID
up_construct_no_reverse = subset_table(up_construct,'is_reverse',0);
dn_construct_no_reverse = subset_table(dn_construct,'is_reverse',0);
common_regions = intersect([up_construct_no_reverse.region_id],[dn_construct_no_reverse.region_id]);
up_construct_no_reverse = subset_table(up_construct,'region_id',common_regions);
dn_construct_no_reverse = subset_table(dn_construct,'region_id',common_regions);
figure
assert(isequal(up_construct_no_reverse{:,'region_id'},dn_construct_no_reverse{:,'region_id'}),'regions not in same order')
gscatter(up_construct_no_reverse{:,'median_ratio'},dn_construct_no_reverse{:,'median_ratio'},up_construct_no_reverse{:,'relation_to_tss'})
xlabel('Upstream Sequence')
ylabel('Downstream Sequence')
title('Promoter Activity (Construct Barcode Ratio)')
ax = gca;
max_lim = max([ax.XLim ax.YLim]);
min_lim = min([ax.XLim ax.YLim]);

ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];
hold on
plot([min_lim max_lim],[min_lim max_lim],'r:')
l = legend;
set(l, 'Interpreter', 'None')
grid on
saveas(gcf,'~/Documents/mpra/fig/up_dn_analysis/promoter_activity_scatter_by_region','png')