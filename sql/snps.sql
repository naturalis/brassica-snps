-- all files we're loading will be tab separated
.separator "\t"

create table if not exists chromosomes (
        chromosome_id integer constraint chromosome_pk primary key asc autoincrement,
	chromosome_name text,
        centromere_feat_start integer, -- index
        centromere_feat_end integer, -- index
        centromere_hta_start integer, -- index
        centromere_hta_end integer --index
);
.import chromosomes.tsv chromosomes

create table if not exists features (
	feature_id integer constraint feature_pk primary key asc autoincrement,
	chromosome_id integer constraint chromosome_fk references chromosomes (chromosome_id),
	source_name text,
	feature_type text,
	feat_start integer, --index
	feat_end integer, --index
	score text,
	strand text,
	phase text,
	attributes text
);
.import features.tsv features

create table if not exists linkages (
	linkage_id integer constraint linkage_pk primary key asc autoincrement,
	chromosome_id integer constraint chromosome_fk references chromosomes (chromosome_id),
	position integer,
	centi_morgan real,
	marker_name text,
	marker_type text,
	fw_primer text,
	rev_primer text,
	feat_start integer, --index
	feat_end integer --index
);
.import linkages.tsv linkages

create table if not exists qtl_regions (
	qtl_region_id integer constraint qtl_region_pk primary key asc autoincrement,
        chromosome_id integer constraint chromosome_fk references chromosomes (chromosome_id),
	qtl integer,
	start integer,
	end integer,
	length integer,
	n_snps integer,
	avg_snps_mb real,
	peak_delta_snp real,
	avg_delta_snp real,
	max_g_prime real,
	mean_g_prime real,
	sd_g_prime real,
	auc_a_t real,
	mean_p_val real,
	mean_q_val real,
	bsa_contrast text
);
.import qtl_regions.tsv qtl_regions

create table if not exists snps ( 
	snp_id integer constraint snp_pk primary key asc autoincrement, -- 1
	chromosome_id integer constraint chromosome_fk references chromosomes (chromosome_id), -- 3
	position integer, -- 13496034
	ref text, -- A
	alt text, -- T
	ad_ref_low integer, -- 35
	ad_alt_low integer, -- 37
	dp_low integer, -- 72
	pl_low text, -- "1253,0,1225"
	gq_low integer, -- 99
	snp_index_low real, -- 0.513888888888889
	ad_ref_high integer, -- 36
	ad_alt_high integer, -- 33
	dp_high integer, -- 69
	pl_high text, -- "1208,0,1404"
	gq_high integer, -- 99
	snp_index_high real, -- 0.478260869565217
	ref_freq real, -- 0.50354609929078
	delta_snp real, -- -0.0356280193236714
	n_snps integer, -- 10626
	tricube_delta_snp real, -- 0.0412010241428342
	g real, -- 0.178946352691828
	g_prime real, -- 2.50042210787922
	p_value real, -- 0.0903783142874682
	neg_log10_p_value real, -- 1.04393576327565
	q_value real -- 0.55772879456407
);
.import snps.tsv snps
