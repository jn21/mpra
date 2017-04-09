mpra_data = readtable('~/Documents/mpra/data/mpra_processed_data_with_annot.txt','Delimiter','\t');

%% Number of regions from Enhancers and from Promoters
unique_up_region = unique(mpra_data{:,'upstream_region_id'});
unique_dn_region = unique(mpra_data{:,'dnstream_region_id'});

num_regions = length(unique(union(unique_up_region,unique_dn_region)))

T = readtable('/Users/jnasser/Documents/mpra/data/prefix_is_from_promoter_annot.txt','Delimiter','\t');

num_promoter_derived = nnz(T{:,'is_within_100nt_of_tss'})
num_enhancer_derived = nnz(~T{:,'is_within_100nt_of_tss'})

%% Number of unique upstream id's

up_rev = subset_table(mpra_data,'upstream_is_reverse',1);
dn_rev = subset_table(mpra_data,'dnstream_is_reverse',1);
num_unique_upstream_regions_with_reverse = length(unique(up_rev{:,'upstream_region_id'}))
num_unique_dnstream_regions_with_reverse = length(unique(dn_rev{:,'dnstream_region_id'}))

num_unique_upstream_sequences = length(unique(mpra_data{:,'upstream_full_id'}))
num_unique_dnstream_sequences = length(unique(mpra_data{:,'dnstream_full_id'}))


