use utf8;
package My::Brassica::Result::Linkage;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

My::Brassica::Result::Linkage

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<linkages>

=cut

__PACKAGE__->table("linkages");

=head1 ACCESSORS

=head2 linkage_id

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

=head2 centi_morgan

  data_type: 'real'
  is_nullable: 1

=head2 marker_name

  data_type: 'text'
  is_nullable: 1

=head2 marker_type

  data_type: 'text'
  is_nullable: 1

=head2 fw_primer

  data_type: 'text'
  is_nullable: 1

=head2 rev_primer

  data_type: 'text'
  is_nullable: 1

=head2 feat_start

  data_type: 'integer'
  is_nullable: 1

=head2 feat_end

  data_type: 'integer'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "linkage_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "chromosome_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
  "position",
  { data_type => "integer", is_nullable => 1 },
  "centi_morgan",
  { data_type => "real", is_nullable => 1 },
  "marker_name",
  { data_type => "text", is_nullable => 1 },
  "marker_type",
  { data_type => "text", is_nullable => 1 },
  "fw_primer",
  { data_type => "text", is_nullable => 1 },
  "rev_primer",
  { data_type => "text", is_nullable => 1 },
  "feat_start",
  { data_type => "integer", is_nullable => 1 },
  "feat_end",
  { data_type => "integer", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</linkage_id>

=back

=cut

__PACKAGE__->set_primary_key("linkage_id");

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
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:1GJo0bt7YwG2PnlnwGCcKg


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
