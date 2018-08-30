use utf8;
package My::Brassica::Result::Feature;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

My::Brassica::Result::Feature

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<features>

=cut

__PACKAGE__->table("features");

=head1 ACCESSORS

=head2 feature_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 chromosome_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 1

=head2 source_name

  accessor: undef
  data_type: 'text'
  is_nullable: 1

=head2 feature_type

  data_type: 'text'
  is_nullable: 1

=head2 feat_start

  data_type: 'integer'
  is_nullable: 1

=head2 feat_end

  data_type: 'integer'
  is_nullable: 1

=head2 score

  data_type: 'text'
  is_nullable: 1

=head2 strand

  data_type: 'text'
  is_nullable: 1

=head2 phase

  data_type: 'text'
  is_nullable: 1

=head2 attributes

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "feature_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "chromosome_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
  "source_name",
  { accessor => undef, data_type => "text", is_nullable => 1 },
  "feature_type",
  { data_type => "text", is_nullable => 1 },
  "feat_start",
  { data_type => "integer", is_nullable => 1 },
  "feat_end",
  { data_type => "integer", is_nullable => 1 },
  "score",
  { data_type => "text", is_nullable => 1 },
  "strand",
  { data_type => "text", is_nullable => 1 },
  "phase",
  { data_type => "text", is_nullable => 1 },
  "attributes",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</feature_id>

=back

=cut

__PACKAGE__->set_primary_key("feature_id");

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
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:E2imNZwVW5L3JRjypkuLug


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
