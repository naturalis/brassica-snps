use utf8;
package My::Brassica::Result::Chromosome;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

My::Brassica::Result::Chromosome

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<chromosomes>

=cut

__PACKAGE__->table("chromosomes");

=head1 ACCESSORS

=head2 chromosome_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 chromosome_name

  data_type: 'text'
  is_nullable: 1

=head2 centromere_feat_start

  data_type: 'integer'
  is_nullable: 1

=head2 centromere_feat_end

  data_type: 'integer'
  is_nullable: 1

=head2 centromere_hta_start

  data_type: 'integer'
  is_nullable: 1

=head2 centromere_hta_end

  data_type: 'integer'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "chromosome_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "chromosome_name",
  { data_type => "text", is_nullable => 1 },
  "centromere_feat_start",
  { data_type => "integer", is_nullable => 1 },
  "centromere_feat_end",
  { data_type => "integer", is_nullable => 1 },
  "centromere_hta_start",
  { data_type => "integer", is_nullable => 1 },
  "centromere_hta_end",
  { data_type => "integer", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</chromosome_id>

=back

=cut

__PACKAGE__->set_primary_key("chromosome_id");

=head1 RELATIONS

=head2 features

Type: has_many

Related object: L<My::Brassica::Result::Feature>

=cut

__PACKAGE__->has_many(
  "features",
  "My::Brassica::Result::Feature",
  { "foreign.chromosome_id" => "self.chromosome_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 kale_snps

Type: has_many

Related object: L<My::Brassica::Result::KaleSnp>

=cut

__PACKAGE__->has_many(
  "kale_snps",
  "My::Brassica::Result::KaleSnp",
  { "foreign.chromosome_id" => "self.chromosome_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 linkages

Type: has_many

Related object: L<My::Brassica::Result::Linkage>

=cut

__PACKAGE__->has_many(
  "linkages",
  "My::Brassica::Result::Linkage",
  { "foreign.chromosome_id" => "self.chromosome_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 qtl_regions

Type: has_many

Related object: L<My::Brassica::Result::QtlRegion>

=cut

__PACKAGE__->has_many(
  "qtl_regions",
  "My::Brassica::Result::QtlRegion",
  { "foreign.chromosome_id" => "self.chromosome_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 snps

Type: has_many

Related object: L<My::Brassica::Result::Snp>

=cut

__PACKAGE__->has_many(
  "snps",
  "My::Brassica::Result::Snp",
  { "foreign.chromosome_id" => "self.chromosome_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07049 @ 2018-11-09 14:13:59
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:hhw649aCsNPIiiusLzB9VQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
