use utf8;
package My::Brassica::Result::KaleSnp;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

My::Brassica::Result::KaleSnp

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<kale_snps>

=cut

__PACKAGE__->table("kale_snps");

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


# Created by DBIx::Class::Schema::Loader v0.07049 @ 2018-11-09 14:13:59
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:fGbLg4SNHd+DElRvSDEQ2Q


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
