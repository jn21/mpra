% fid = fopen('~/Documents/mpra/data/mpra_data.txt');
% 
% A = textscan(fid)

% mpra_data = readtable('~/Documents/mpra/data/mpra_data.txt','Delimiter','\t');
% construct_annot = readtable('~/Documents/mpra/annot/mpra_data_and_construct_annot.txt','Delimiter','\t');
% combined_data = readtable('/Users/jnasser/Documents/mpra/annot/mpra_data_and_construct_annot.txt','Delimiter','\t');

mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');

%% Distributions of input DNA
% figure
% scatter(mpra_data{:,'DNAInput_ETotal'},mpra_data{:,'DNAInput_PTotal'})
% xlabel('DNAInput_ETotal','Interpreter','None')
% ylabel('DNAInput_PTotal','Interpreter','None')
% grid on
% ax = gca;
% max_lim = max([ax.YLim ax.XLim]);
% hold on
% plot([0 max_lim],[0 max_lim],'r:');

figure
scatter(log2(mpra_data{:,'DNAInput_ETotal'}),log2(mpra_data{:,'DNAInput_PTotal'}))
xlabel('log2(DNAInput_ETotal)','Interpreter','None')
ylabel('log2(DNAInput_PTotal)','Interpreter','None')
grid on
ax = gca;
max_lim = max([ax.YLim ax.XLim]);
hold on
plot([0 max_lim],[0 max_lim],'r:');

%% Make a quick annotation table describing if each prefix is from an enhancer or promoter
% up_prefixes = mpra_data{:,'upstream_prefix'};
% up_is_within_100nt_of_tss = mpra_data{:,'upstream_sequence_wt_100nt_of_tss'};
% dn_prefixes = mpra_data{:,'dnstream_prefix'};
% dn_is_within_100nt_of_tss = mpra_data{:,'dnstream_sequence_wt_100nt_of_tss'};
% annot_mat = unique([vertcat(up_prefixes,dn_prefixes) vertcat(up_is_within_100nt_of_tss,dn_is_within_100nt_of_tss)],'rows');
% annot_table = array2table(annot_mat,'VariableNames',{'prefix','is_within_100nt_of_tss'});
% writetable(annot_table,'~/Documents/mpra/data/prefix_is_from_promoter_annot.txt','Delimiter','\t')

%% Examine individual distributions 
no_mods = subset_table(mpra_data,'dnstream_is_modified',0);
unique_upstream_ids = unique(no_mods{:,'upstream_full_id'});

for ii = 101:120
    figure
    this_table = subset_table(no_mods, 'upstream_full_id', unique_upstream_ids{ii});
    histogram(this_table{:,'E_ratio_avg_rep'},...
        'BinWidth',.1)
end


%% Get some numbers

length(unique(mpra_data{:,'dnstream_id'}))
length(unique(mpra_data{:,'upstream_id'}))

%% Dnstream signals
signals = {'addStrongPAS','addPAS','delPAS','addU1'};
ratios = {'E_ratio_avg_rep','P_ratio_avg_rep'};

for ii = 1:length(signals)
    for jj = 1:length(ratios)
        explore_dnstream_sequence_signal(signals{ii},ratios{jj},true,true)
    end
end

%% Intact pairs
intact_idx = mpra_data{:,'upstream_region_id'} == mpra_data{:,'dnstream_region_id'};
normal_idx = mpra_data{:,'upstream_is_reverse'} == 0 & ...
    mpra_data{:,'dnstream_is_reverse'} == 0 & ...
    mpra_data{:,'dnstream_addPAS'} == 0 & ...
    mpra_data{:,'dnstream_delPAS'} == 0 & ...
    mpra_data{:,'dnstream_addStrongPAS'} == 0 & ...
    mpra_data{:,'dnstream_num_delU1'} == 0 & ...
    mpra_data{:,'dnstream_addU1'} == 0;

finite_idx = isfinite(mpra_data{:,'E_ratio_avg_rep'}) & isfinite(mpra_data{:,'P_ratio_avg_rep'});
to_use_idx = intact_idx & normal_idx & finite_idx;

% wt_100bp_idx = strcmp(mpra_data{:,'upstream_sequence_relation_to_tss'},'wt_100nt');
% 
% boxplot(horzcat(mpra_data{to_use_idx & wt_100bp_idx,'P_ratio_avg_rep'},...
%     mpra_data{to_use_idx & ~wt_100bp_idx,'P_ratio_avg_rep'}))

scatter(mpra_data{:,'E_ratio_avg_rep'},mpra_data{:,'P_ratio_avg_rep'})
hold on
scatter(mpra_data{to_use_idx,'E_ratio_avg_rep'},mpra_data{to_use_idx,'P_ratio_avg_rep'})

corr(mpra_data{finite_idx & normal_idx,'E_ratio_avg_rep'},mpra_data{finite_idx & normal_idx,'P_ratio_avg_rep'})
corr(mpra_data{to_use_idx,'E_ratio_avg_rep'},mpra_data{to_use_idx,'P_ratio_avg_rep'})
            
%% SignRank 
DIMENSION = 10000;
TRIALS = 1000;

for ii = 1:TRIALS
    ii
    x = rand(DIMENSION,1);
    y = rand(DIMENSION,1);
    
    pval(ii) = signrank(x,y);
end

histogram(pval)
% %% Replicate Correlation
% figure
% scatter(mpra_data{:,'Rep1_ERatio'},mpra_data{:,'Rep2_ERatio'})
% figure
% scatter(mpra_data{:,'Rep1_PRatio'},mpra_data{:,'Rep2_PRatio'})
% 
% %% E and P Activity
% figure
% scatter(mpra_data{:,'Rep1_ERatio'},mpra_data{:,'Rep1_PRatio'})
% figure
% scatter(mpra_data{:,'Rep2_ERatio'},mpra_data{:,'Rep2_PRatio'})

% %% Remove some low 'E' count entries
% CUTOFF = 20;
% rep1_idx = (mpra_data{:,'Rep1_ETotal'} < CUTOFF);
% rep2_idx = (mpra_data{:,'Rep2_ETotal'} < CUTOFF);
% mpra_data = mpra_data(~(rep1_idx | rep2_idx),:);
% 
% %% Linear regression
% lm_table = convert_data_table_to_lm_table(mpra_data);
% 
% % this_table = lm_table;
% pratio_table = lm_table; 
% pratio_table.E_ratio = [];
% 
% eratio_table = lm_table; 
% eratio_table.P_ratio = [];
% 
% pratio_lm = fitlm(pratio_table,'ResponseVar','P_ratio');
% eratio_lm = fitlm(eratio_table,'ResponseVar','E_ratio');
% 
% %% Weird shit
% A = table2array(pratio_table);
% [m,n] = size(A);
% A(:,n) = [];
% 
% A_pseudo = pinv(A);
% A_pseudo*lm_table{:,'P_ratio'};
% A_pseudo*lm_table{:,'E_ratio'};