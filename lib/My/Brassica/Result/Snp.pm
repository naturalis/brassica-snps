use utf8;
package My::Brassica::Result::Snp;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

My::Brassica::Result::Snp

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<snps>

=cut

__PACKAGE__->table("snps");

=head1 ACCESSORS

=head2 snp_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 chromosome_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 1

=head2 position

  data_type: 'integer'
  is_nullable: 1

=head2 ref

  data_type: 'text'
  is_nullable: 1

=head2 alt

  data_type: 'text'
  is_nullable: 1

=head2 ad_ref_low

  data_type: 'integer'
  is_nullable: 1

=head2 ad_alt_low

  data_type: 'integer'
  is_nullable: 1

=head2 dp_low

  data_type: 'integer'
  is_nullable: 1

=head2 gq_low

  data_type: 'integer'
  is_nullable: 1

=head2 snp_index_low

  data_type: 'real'
  is_nullable: 1

=head2 ad_ref_high

  data_type: 'integer'
  is_nullable: 1

=head2 ad_alt_high

  data_type: 'integer'
  is_nullable: 1

=head2 dp_high

  data_type: 'integer'
  is_nullable: 1

=head2 gq_high

  data_type: 'integer'
  is_nullable: 1

=head2 snp_index_high

  data_type: 'real'
  is_nullable: 1

=head2 ref_freq

  data_type: 'real'
  is_nullable: 1

=head2 delta_snp

  data_type: 'real'
  is_nullable: 1

=head2 n_snps

  data_type: 'integer'
  is_nullable: 1

=head2 tricube_delta_snp

  data_type: 'real'
  is_nullable: 1

=head2 g

  data_type: 'real'
  is_nullable: 1

=head2 g_prime

  data_type: 'real'
  is_nullable: 1

=head2 p_value

  data_type: 'real'
  is_nullable: 1

=head2 neg_log10_p_value

  data_type: 'real'
  is_nullable: 1

=head2 q_value

  data_type: 'real'
  is_nullable: 1

=head2 contrast

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "snp_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "chromosome_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
  "position",
  { data_type => "integer", is_nullable => 1 },
  "ref",
  { data_type => "text", is_nullable => 1 },
  "alt",
  { data_type => "text", is_nullable => 1 },
  "ad_ref_low",
  { data_type => "integer", is_nullable => 1 },
  "ad_alt_low",
  { data_type => "integer", is_nullable => 1 },
  "dp_low",
  { data_type => "integer", is_nullable => 1 },
  "gq_low",
  { data_type => "integer", is_nullable => 1 },
  "snp_index_low",
  { data_type => "real", is_nullable => 1 },
  "ad_ref_high",
  { data_type => "integer", is_nullable => 1 },
  "ad_alt_high",
  { data_type => "integer", is_nullable => 1 },
  "dp_high",
  { data_type => "integer", is_nullable => 1 },
  "gq_high",
  { data_type => "integer", is_nullable => 1 },
  "snp_index_high",
  { data_type => "real", is_nullable => 1 },
  "ref_freq",
  { data_type => "real", is_nullable => 1 },
  "delta_snp",
  { data_type => "real", is_nullable => 1 },
  "n_snps",
  { data_type => "integer", is_nullable => 1 },
  "tricube_delta_snp",
  { data_type => "real", is_nullable => 1 },
  "g",
  { data_type => "real", is_nullable => 1 },
  "g_prime",
  { data_type => "real", is_nullable => 1 },
  "p_value",
  { data_type => "real", is_nullable => 1 },
  "neg_log10_p_value",
  { data_type => "real", is_nullable => 1 },
  "q_value",
  { data_type => "real", is_nullable => 1 },
  "contrast",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</snp_id>

=back

=cut

__PACKAGE__->set_primary_key("snp_id");

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


# Created by DBIx::Class::Schema::Loader v0.07049 @ 2018-09-07 14:08:56
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:/3rWLPEtowduxQAbBGv8Cw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
