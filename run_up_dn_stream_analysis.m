
%% Run ID analysis for promoter and enhancer activity
up_gfp = examine_activity_by_up_or_dn_id('up','gfp');
dn_gfp = examine_activity_by_up_or_dn_id('dn','gfp');
up_construct = examine_activity_by_up_or_dn_id('up','construct');
dn_construct = examine_activity_by_up_or_dn_id('dn','construct');

%% Scatter Enhancer and Promoter Activity of UP and DN stream sequences
figure
scatter(up_gfp{:,'median_ratio'},up_construct{:,'median_ratio'})
hold on
ax = gca;
plot([up_gfp{1,'neg_control_median_ratio'} up_gfp{1,'neg_control_median_ratio'}], ax.YLim, 'r:')
hold on
plot(ax.XLim, [up_construct{1,'neg_control_median_ratio'} up_construct{1,'neg_control_median_ratio'}], 'r:')
xlabel('Median Enhancer Activity')
ylabel('Median Promoter Activity')
title('Upstream Sequences')

figure
scatter(dn_gfp{:,'median_ratio'},dn_construct{:,'median_ratio'})
hold on
ax = gca;
plot([dn_gfp{1,'neg_control_median_ratio'} dn_gfp{1,'neg_control_median_ratio'}], ax.YLim, 'r:')
hold on
plot(ax.XLim, [dn_construct{1,'neg_control_median_ratio'} dn_construct{1,'neg_control_median_ratio'}], 'r:')
xlabel('Median Enhancer Activity')
ylabel('Median Promoter Activity')
title('Downstream Sequences')

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