package My::GeneDB::Record;
use strict;
use warnings;
our $AUTOLOAD;

sub new {
    my ( $package, %record ) = @_;
    return bless \%record, $package;
}

sub link_name {
    my ( $self, $field ) = @_;
    my $value = $self->$field;
    if ( $value =~ /^\[(.+?)\]/ ) {
        my $name = $1;
        return $name;
    }
    die $field;
}

sub to_markdown {
    my ( $self, @columns ) = @_;
    my @values = map { $self->$_ } @columns;
    no warnings 'uninitialized';
    return '|' . join( '|', @values ) . '|';
}

sub AUTOLOAD {
    my ( $self, $arg ) = @_;
    my $method = $AUTOLOAD;
    $method =~ s/.+://;
    if ( exists $self->{$method} ) {
        if ( defined $arg ) {
            $self->{$method} = $arg;
            return $self;
        }
        else {
            return $self->{$method};
        }
    }
    else {
        if ( $method !~ /^[A-Z]+$/ ) {
            die "No method $method";
        }
    }
}

package My::GeneDB;
use strict;
use warnings;
our $AUTOLOAD;

sub new {
    my ( $package, $file ) = @_;
    my $records = $package->read($file);
    return bless $records, $package;
}

sub read {
    my ( $package, $file ) = @_;
    my ( @records, @header );
    open my $in, '<', $file or die $!;
    LINE: while(<$in>) {
        chomp;
        my @line = split /\|/, $_;
        shift @line;
        #pop @line;
        if ( not @header ) {
            @header = @line;
            next LINE;
        }
        elsif ( /^[\|\-]+$/ ) {
            next LINE;
        }
        else {
            my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
            push @records, My::GeneDB::Record->new(%record);
        }
    }

    # ->columns and ->records are the methods
    return {
        columns => \@header,
        records => \@records,
    };
}

sub to_markdown {
    my $self = shift;
    my @lines;
    push @lines, '|' . join( '|', $self->columns ) . '|';
    push @lines, '|' . join( '|', map { '-' x length($_) } $self->columns ) . '|';
    push @lines, join( "\n", map { $_->to_markdown($self->columns) } $self->records );
    return join "\n", @lines;
}

sub AUTOLOAD {
    my ( $self, @arg ) = @_;
    my $method = $AUTOLOAD;
    $method =~ s/.+://;
    if ( exists $self->{$method} ) {
        if ( @arg ) {
            $self->{$method} = [ @arg ];
            return $self;
        }
        else {
            return @{ $self->{$method} };
        }
    }
    else {
        if ( $method !~ /^[A-Z]+$/ ) {
            die "No method $method";
        }
    }
}

1;