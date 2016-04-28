#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;

my @files;
my $mapfile;
my $outfile;

GetOptions(		'files|f=s{1,}' => \@files,
				'mapfile|m=s' => \$mapfile,
				'outfile|o=s' => \$outfile,

	);


open(MAP, "<", $mapfile) or die $!;

my %lgs = read_map($mapfile);


my %groups;
foreach my $lg (keys %lgs) {
	foreach my $locus (keys $lgs{$lg}) {
		$groups{$locus} = $lg;
	}
}

my %blocks;

foreach my $file (@files) {
	open(IN, "<", $file) or die $!;
	<IN>;
	my @chunks = split(/\./, $file);
	my $genome = $chunks[0];
	while(<IN>) {
		last if $_ =~ /^\s/;
		my ($locus, $lg, $left, $right, $left_pos, $right_pos) = split;
		$blocks{$genome}{$locus} = [$left, $right];
	}
	close IN;
}

open(DUMP, ">", 'dumper.out') or die $!;
#print DUMP Dumper(\%blocks);

my %comb_blocks;
foreach my $genome (keys %blocks) {
	foreach my $locus (keys $blocks{$genome}) {
		if (! $comb_blocks{$locus}) {
			$comb_blocks{$locus} = [ [$blocks{$genome}{$locus}[0], $blocks{$genome}{$locus}[1]], [$genome] ];
		} elsif ($comb_blocks{$locus}) {
			# compare the two blocks to find the smaller interval

			# the existing block:

			my $e_left = $comb_blocks{$locus}[0][0];
			my $e_right = $comb_blocks{$locus}[0][1];
			my $e_int = lookup_interval($e_left, $e_right);

			# the new proposed block

			my $n_left = $blocks{$genome}{$locus}[0];
			my $n_right = $blocks{$genome}{$locus}[1];
			my $n_int = lookup_interval($n_left, $n_right);

			if ($e_int <= $n_int) {
				push @{$comb_blocks{$locus}[1]}, $genome;
			} else {
				$comb_blocks{$locus}[0][0] = $n_left;
				$comb_blocks{$locus}[0][1] = $n_right;
				push @{$comb_blocks{$locus}[1]}, $genome;
			}

		}
	}
}

print DUMP Dumper(\%comb_blocks);

open(OUT, ">", $outfile) or die $!;
print OUT join("\t", 'LOCUS', 'LEFT_LOCUS', 'RIGHT_LOCUS', 'LEFT_POS', 'RIGHT_POS', 'SUPPORT', 'INTERVAL_SIZE', 'LG'), "\n";
foreach my $locus (keys %comb_blocks) {

	my $map_lg;
	if ($groups{$comb_blocks{$locus}[0][0]} eq $groups{$comb_blocks{$locus}[0][1]}) {
		$map_lg = $groups{$comb_blocks{$locus}[0][0]};
	} else {
		print "Error mapping locus $locus - flanking loci on different linkage groups\n";
	}
	print OUT join("\t", $locus, $comb_blocks{$locus}[0][0], $comb_blocks{$locus}[0][1], $lgs{$map_lg}{$comb_blocks{$locus}[0][0]}, $lgs{$map_lg}{$comb_blocks{$locus}[0][1]}, join(',', @{$comb_blocks{$locus}[1]}), lookup_interval($comb_blocks{$locus}[0][0], $comb_blocks{$locus}[0][1]), $map_lg), "\n";
}

sub lookup_interval {
	my $left = shift;
	my $right = shift;

	my $lpos;
	my $rpos;
	foreach my $group (keys %lgs) {
		if ($lgs{$group}{$left}) {
			$lpos = $lgs{$group}{$left};
		}
		if ($lgs{$group}{$right}) {
			$rpos = $lgs{$group}{$right};
		}
	}
	my $interval = abs($lpos - $rpos);
	return $interval;
}

sub read_map {
	my $mapfile = shift;
	open(MAP, "<", $mapfile) or die $!;
	my %lgs;
	my $group;
	<MAP>;
	while(<MAP>) {
		next if $_ =~ /^\s/;
		chomp;
		my ($locus, $lg, $pos) = split;
		$lgs{$lg}{$locus} = $pos;

	}
	close MAP;

	return %lgs;
}
