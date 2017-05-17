mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot_normalized.txt','Delimiter','\t');

%% Run ID analysis for promoter and enhancer activity
save_fig = true;
enhancer_activity_by_up = examine_activity_by_up_or_dn_id('up','enhancer_activity',save_fig);
enhancer_activity_by_dn = examine_activity_by_up_or_dn_id('dn','enhancer_activity',save_fig);
promoter_activity_by_up = examine_activity_by_up_or_dn_id('up','promoter_activity',save_fig);
promoter_activity_by_dn = examine_activity_by_up_or_dn_id('dn','promoter_activity',save_fig);

%% Scatter E and P activity of Native sequences
bad_idx = mpra_data{:,'dnstream_is_reverse'} |...
    mpra_data{:,'upstream_is_reverse'} |...
    mpra_data{:,'dnstream_is_modified'} |...
    ~isfinite(mpra_data{:,'enhancer_activity'}) |...
    ~isfinite(mpra_data{:,'promoter_activity'});
good_mpra_data = mpra_data(~bad_idx,:);

intact_data = subset_table(good_mpra_data,'is_intact_sequence',1);
[intact_corr, intact_pval] = corr(intact_data{:,'enhancer_activity'},intact_data{:,'promoter_activity'});
scatter_corr = corr(good_mpra_data{:,'enhancer_activity'},good_mpra_data{:,'promoter_activity'});

figure;
scatter(intact_data{:,'enhancer_activity'},...
    intact_data{:,'promoter_activity'},...
    'DisplayName','200bp Native Sequences');
% gscatter(intact_data{:,'enhancer_activity'},...
%     intact_data{:,'promoter_activity'},...
%     intact_data{:,'upstream_sequence_derivation'});
xlabel('Enhancer Activity')
ylabel('Promoter Activity')
title_str = sprintf('Promoter vs Enhancer Activity \n Pearson correlation = %.2f',...
    intact_corr);
title(title_str)
grid on
ax = gca;
ax.XLim = [-2 1];
ax.YLim = [-8 0];
l = legend;
set(l,'Interpreter','None')
%saveas(gcf,'~/Documents/mpra/fig/ep_scatter/native_sequences','png')

%% Scatter E and P activity of all constructs
figure;
scatter(good_mpra_data{:,'enhancer_activity'},...
    good_mpra_data{:,'promoter_activity'},...
    'DisplayName','All Constructs');
hold on
scatter(intact_data{:,'enhancer_activity'},...
    intact_data{:,'promoter_activity'},...
    'filled',...
    'MarkerEdgeColor','flat',...
    'DisplayName','200bp Native Sequences');
xlabel('Enhancer Activity')
ylabel('Promoter Activity')
title_str = sprintf(['Promoter vs Enhancer Activity \n' ...
    'Native Sequence Construct - correlation = %.2f \n'...
    'All Non-modified Constructs - correlation = %.2f'],...
    scatter_corr,intact_corr);
title(title_str)
legend show
saveas(gcf,'~/Documents/mpra/fig/ep_scatter/all_non_modified_sequences','png')

%% Scatter Enhancer and Promoter Activity of Median UP and DN stream sequences
up_gfp_no_rev = subset_table(enhancer_activity_by_up,'is_reverse',0);
dn_gfp_no_rev = subset_table(enhancer_activity_by_dn,'is_reverse',0);
up_construct_no_rev = subset_table(promoter_activity_by_up,'is_reverse',0);
dn_construct_no_rev = subset_table(promoter_activity_by_dn,'is_reverse',0);

figure
scatter(up_gfp_no_rev{:,'median_ratio'},up_construct_no_rev{:,'median_ratio'})
% gscatter(up_gfp_no_rev{:,'median_ratio'},...
%     up_construct_no_rev{:,'median_ratio'},...
%     up_construct_no_rev{:,'derivation'})
scatter_corr = corr(up_gfp_no_rev{:,'median_ratio'},up_construct_no_rev{:,'median_ratio'});
% hold on
% ax = gca;
% plot([up_gfp{1,'neg_control_median_ratio'} up_gfp{1,'neg_control_median_ratio'}], ax.YLim, 'r:')
% hold on
% plot(ax.XLim, [up_construct{1,'neg_control_median_ratio'} up_construct{1,'neg_control_median_ratio'}], 'r:')
xlabel('Median Enhancer Activity')
ylabel('Median Promoter Activity')
title_str = sprintf('Upstream Sequences - Medians \n Correlation = %.2f',scatter_corr);
title(title_str)
%saveas(gcf,'~/Documents/mpra/fig/ep_scatter/up_ep_median_scatter','png')
l = legend;
set(l,'Interpreter','None')

figure
scatter(dn_gfp_no_rev{:,'median_ratio'},dn_construct_no_rev{:,'median_ratio'})
% gscatter(dn_gfp_no_rev{:,'median_ratio'},...
%     dn_construct_no_rev{:,'median_ratio'},...
%     dn_construct_no_rev{:,'derivation'})
scatter_corr = corr(dn_gfp_no_rev{:,'median_ratio'},dn_construct_no_rev{:,'median_ratio'});
% hold on
% ax = gca;
% plot([dn_gfp{1,'neg_control_median_ratio'} dn_gfp{1,'neg_control_median_ratio'}], ax.YLim, 'r:')
% hold on
% plot(ax.XLim, [dn_construct{1,'neg_control_median_ratio'} dn_construct{1,'neg_control_median_ratio'}], 'r:')
xlabel('Median Enhancer Activity')
ylabel('Median Promoter Activity')
title_str = sprintf('Downstream Sequences - Medians \n Correlation = %.2f',scatter_corr);
title(title_str)
%saveas(gcf,'~/Documents/mpra/fig/ep_scatter/dn_ep_median_scatter','png')
l = legend;
set(l,'Interpreter','None')

%% FDR 
fdr_cutoff = .1;
nnz(enhancer_activity_by_up{:,'qval'} < fdr_cutoff)
nnz(enhancer_activity_by_dn{:,'qval'} < fdr_cutoff)
nnz(promoter_activity_by_up{:,'qval'} < fdr_cutoff)
nnz(promoter_activity_by_dn{:,'qval'} < fdr_cutoff)

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
data = promoter_activity_by_up{:,'median_ratio'} ./ enhancer_activity_by_up{:,'median_ratio'};
boxplot(data,promoter_activity_by_up{:,'median_ratio'});

%% Enhancer Activity by region ID
up_enhancer_no_reverse = subset_table(enhancer_activity_by_up,'is_reverse',0);
dn_enhancer_no_reverse = subset_table(enhancer_activity_by_dn,'is_reverse',0);
common_regions = intersect([up_enhancer_no_reverse.region_id],[dn_enhancer_no_reverse.region_id]);
up_enhancer_no_reverse = subset_table(enhancer_activity_by_up,'region_id',common_regions);
dn_enhancer_no_reverse = subset_table(enhancer_activity_by_dn,'region_id',common_regions);
figure
assert(isequal(up_enhancer_no_reverse{:,'region_id'},dn_enhancer_no_reverse{:,'region_id'}),'regions not in same order')
gscatter(up_enhancer_no_reverse{:, 'median_ratio'}, dn_enhancer_no_reverse{:, 'median_ratio'},up_enhancer_no_reverse{:, 'derivation'})
%scatter(up_enhancer_no_reverse{:, 'median_ratio'}, dn_enhancer_no_reverse{:, 'median_ratio'})

xlabel('Upstream Sequence')
ylabel('Downstream Sequence')
title('Enhancer Activity')
ax = gca;
max_lim = max([ax.XLim ax.YLim]);
min_lim = min([ax.XLim ax.YLim]);

ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];
% hold on
% plot([min_lim max_lim],[min_lim max_lim],'r:')
l = legend;
set(l, 'Interpreter', 'None')
grid on
%saveas(gcf,'~/Documents/mpra/fig/ep_scatter/enhancer_activity_scatter_by_region','png')

%% Promoter Activity by region ID
up_promoter_no_reverse = subset_table(promoter_activity_by_up,'is_reverse',0);
dn_promoter_no_reverse = subset_table(promoter_activity_by_dn,'is_reverse',0);
common_regions = intersect([up_promoter_no_reverse.region_id],[dn_promoter_no_reverse.region_id]);
up_promoter_no_reverse = subset_table(promoter_activity_by_up,'region_id',common_regions);
dn_promoter_no_reverse = subset_table(promoter_activity_by_dn,'region_id',common_regions);
figure
assert(isequal(up_promoter_no_reverse{:,'region_id'},dn_promoter_no_reverse{:,'region_id'}),'regions not in same order')
%gscatter(up_promoter_no_reverse{:,'median_ratio'},dn_promoter_no_reverse{:,'median_ratio'},up_promoter_no_reverse{:,'derivation'})
scatter(up_promoter_no_reverse{:,'median_ratio'},dn_promoter_no_reverse{:,'median_ratio'})
xlabel('Upstream Sequence')
ylabel('Downstream Sequence')
title('Promoter Activity')
ax = gca;
max_lim = max([ax.XLim ax.YLim]);
min_lim = min([ax.XLim ax.YLim]);

ax.XLim = [min_lim max_lim];
ax.YLim = [min_lim max_lim];
l = legend;
set(l, 'Interpreter', 'None')
% hold on
% plot([min_lim max_lim],[min_lim max_lim],'r:')

grid on
saveas(gcf,'~/Documents/mpra/fig/ep_scatter/promoter_activity_scatter_by_region','png')