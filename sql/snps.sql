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
