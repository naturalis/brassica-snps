use utf8;
package My::Brassica::Result::QtlRegion;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

My::Brassica::Result::QtlRegion

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<qtl_regions>

=cut

__PACKAGE__->table("qtl_regions");

=head1 ACCESSORS

=head2 qtl_region_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 chromosome_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 1

=head2 qtl

  data_type: 'integer'
  is_nullable: 1

=head2 start

  data_type: 'integer'
  is_nullable: 1

=head2 end

  data_type: 'integer'
  is_nullable: 1

=head2 length

  data_type: 'integer'
  is_nullable: 1

=head2 n_snps

  data_type: 'integer'
  is_nullable: 1

=head2 avg_snps_mb

  data_type: 'real'
  is_nullable: 1

=head2 peak_delta_snp

  data_type: 'real'
  is_nullable: 1

=head2 avg_delta_snp

  data_type: 'real'
  is_nullable: 1

=head2 max_g_prime

  data_type: 'real'
  is_nullable: 1

=head2 mean_g_prime

  data_type: 'real'
  is_nullable: 1

=head2 sd_g_prime

  data_type: 'real'
  is_nullable: 1

=head2 auc_a_t

  data_type: 'real'
  is_nullable: 1

=head2 mean_p_val

  data_type: 'real'
  is_nullable: 1

=head2 mean_q_val

  data_type: 'real'
  is_nullable: 1

=head2 bsa_contrast

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "qtl_region_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "chromosome_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
  "qtl",
  { data_type => "integer", is_nullable => 1 },
  "start",
  { data_type => "integer", is_nullable => 1 },
  "end",
  { data_type => "integer", is_nullable => 1 },
  "length",
  { data_type => "integer", is_nullable => 1 },
  "n_snps",
  { data_type => "integer", is_nullable => 1 },
  "avg_snps_mb",
  { data_type => "real", is_nullable => 1 },
  "peak_delta_snp",
  { data_type => "real", is_nullable => 1 },
  "avg_delta_snp",
  { data_type => "real", is_nullable => 1 },
  "max_g_prime",
  { data_type => "real", is_nullable => 1 },
  "mean_g_prime",
  { data_type => "real", is_nullable => 1 },
  "sd_g_prime",
  { data_type => "real", is_nullable => 1 },
  "auc_a_t",
  { data_type => "real", is_nullable => 1 },
  "mean_p_val",
  { data_type => "real", is_nullable => 1 },
  "mean_q_val",
  { data_type => "real", is_nullable => 1 },
  "bsa_contrast",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</qtl_region_id>

=back

=cut

__PACKAGE__->set_primary_key("qtl_region_id");

=head1 RELATIONS

=head2 chromosome

Type: belongs_to

Related object: L<My::Brassica::Result::Chromosome>

=cut

__PACKAGE__->belongs_to(
  "chromosome",
  "My::Brassica::Result::Chromosome",
  { chromosome_id => "chromosome_id" },
  {
    is_deferrable => 0,
    join_type     => "LEFT",
    on_delete     => "NO ACTION",
    on_update     => "NO ACTION",
  },
);


# Created by DBIx::Class::Schema::Loader v0.07049 @ 2018-08-30 12:25:46
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:dk5GpWTUJ88H7z6H6eF+HQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
