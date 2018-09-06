my $pk = 1;
while(<>) {
	my @fields = split /\t/;
	$fields[0] =~ s/C//;
	unshift @fields, $pk++;
	print join "\t", @fields;
}
