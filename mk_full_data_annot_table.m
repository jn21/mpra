
mpra_data = readtable('~/Documents/mpra/data/mpra_data.txt','Delimiter','\t');

construct_names = mpra_data{:,'construct'};

construct_names_parsed = parse_construct_name(construct_names);

T = join(mpra_data,construct_names_parsed,'Keys','construct');