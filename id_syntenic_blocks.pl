#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(firstidx);
use Statistics::Descriptive;

pod2usage(-verbose => 1) if @ARGV == 0;

my $infile;
my $mapfile;
my $genome;

GetOptions(		'infile|i=s' => \$infile,
				'mapfile|m=s' => \$mapfile,
				'genome|g=s' => \$genome,

		);
	
open(CSV, "<", $infile) or die $!;
open(DUMP, ">", 'dumper.out');
open(SUM, ">", $genome . '.synteny.sum');
open(LOG, ">", $genome . '.synteny.log');
open(OUT, ">", $genome . '.synteny.blocks');
open(STAT, ">", $genome . '.synteny.stats');
open(TAB, ">", $genome . '.synteny.tab');


chomp(my $headers = <CSV>);
my @fields = split(',', $headers);

my %map_info;
my %comp_info;

my $chr_field;
my $prop_field;
my $pos_field;
for (my $i = 0; $i < scalar(@fields); $i++) {
	if ($fields[$i] =~ /($genome)_chr/) {
		$chr_field = $i;
	}
	if ($fields[$i] =~ /($genome)_prop/) {
		$prop_field = $i;
	}
	if ($fields[$i] =~ /($genome)_pos/) {
		$pos_field = $i;
	}
}


while(<CSV>) {
	chomp;
	my @fields = split(',', $_);
	$map_info{$fields[0]} = [$fields[1], $fields[2], $fields[3]];
	$comp_info{$fields[0]} = [$fields[$chr_field], $fields[$prop_field], $fields[$pos_field]];
}
close CSV;

my %lgs = read_map($mapfile);



my %map_lgs;
my @map_lgs;
my %pos_by_lg;
foreach my $locus (sort { ($map_info{$a}[0] <=> $map_info{$b}[0]) || ($map_info{$a}[1] <=> $map_info{$b}[1])} keys %comp_info) {
	next unless $comp_info{$locus}[0] && $map_info{$locus}[0] =~ /\d/;
	push @map_lgs, $map_info{$locus}[0] unless $map_lgs{$map_info{$locus}[0]};
	$map_lgs{$map_info{$locus}[0]}++;
	if (! $pos_by_lg{$map_info{$locus}[0]}) {
		$pos_by_lg{$map_info{$locus}[0]} = [ [ $locus, $map_info{$locus}[1], $comp_info{$locus}[0], $comp_info{$locus}[1] ] ];
	} else {
		push @{$pos_by_lg{$map_info{$locus}[0]}}, [ $locus, $map_info{$locus}[1], $comp_info{$locus}[0], $comp_info{$locus}[1] ];
	}
	#print OUT join("\t", $comp_info{$locus}[0], $locus, $map_info{$locus}[0], $map_info{$locus}[1]), "\n";
}
#print join("\n", @comp_lgs), "\n";

#print LOG LOGer(\%pos_by_lg);


my %blocks;
foreach my $lg (@map_lgs) {
	print LOG "-----Current comp LG: $lg-----\n";
	$blocks{$lg} = [];
	my @loci;
	my $pos;
	my $prev_pos;
	my $prev_dir;
	my $map_lg;
	my $prev_map_lg;
	for (my $i = 0; $i < scalar(@{$pos_by_lg{$lg}}); $i++) {
		my $locus = $pos_by_lg{$lg}[$i][0];
		print LOG "Current locus: $locus\n";
		if (scalar(@loci) == 0) { # first locus
			print LOG "Array is empty, starting a new group\n";
			push @loci, $locus;
			$prev_pos = $pos_by_lg{$lg}[$i][3];
			$prev_map_lg = $pos_by_lg{$lg}[$i][2];
			#print join("\t", $locus, $lg, $pos_by_lg{$lg}[$i][3], $pos_by_lg{$lg}[$i][2]), "\n";
			
		} elsif (scalar(@loci) == 1) {
			print LOG "Group has one locus, checking to see if the current locus is syntenic\n";
			my $same_lg = same_lg($pos_by_lg{$lg}[$i][2], $pos_by_lg{$lg}[$i-1][2]);
			my $adjacent = adjacent($i, $pos_by_lg{$lg}[$i][3], $pos_by_lg{$lg}[$i-1][3], $lg);
			if ($same_lg && $adjacent) {
				print LOG "Same LG:$same_lg, Adjacent:$adjacent\n";
				print LOG "Locus is syntenic - adding to group\n";
				push @loci, $locus;
				$prev_dir = get_dir($pos_by_lg{$lg}[$i][3], $pos_by_lg{$lg}[$i-1][3]);
				print LOG "Establishing a direction: $prev_dir\n";
			} else {
				print LOG "Same LG:$same_lg, Adjacent:$adjacent\n";
				print LOG "Locus is not syntenic - starting a new group\n";
				push $blocks{$lg}, [ @loci ];
				@loci = ($locus);
			}
			
		} elsif (scalar(@loci) > 1) {
			print LOG "Group has more than one locus, checking to see if the current locus is syntenic\n";
			my $same_lg = same_lg($pos_by_lg{$lg}[$i][2], $pos_by_lg{$lg}[$i-1][2]);
			my $adjacent = adjacent($i, $pos_by_lg{$lg}[$i][3], $pos_by_lg{$lg}[$i-1][3], $lg);
			my $dir = get_dir($pos_by_lg{$lg}[$i][3], $pos_by_lg{$lg}[$i-1][3]);
			my $same_dir = same_dir($dir, $prev_dir);
			if ($same_lg && $adjacent && $same_dir) {
				print LOG "Same LG:$same_lg, Adjacent:$adjacent, Same dir:$same_dir ($dir, $prev_dir)\n";
				print LOG "Locus is syntenic - adding to group\n";
				push @loci, $locus;
				$prev_dir = $dir;
			} else {
				print LOG "Same LG:$same_lg, Adjacent:$adjacent, Same dir:$same_dir\n";
				if ((abs($pos_by_lg{$lg}[$i][3] - $pos_by_lg{$lg}[$i][3]) <= 0.05) && $same_lg && $adjacent) {
					print LOG "SPECIAL CASE: small difference\n";
				}
				print LOG "Locus is not syntenic - starting a new group\n";
				push $blocks{$lg}, [ @loci ];
				@loci = ($locus);
			}
		}
		
		print LOG "Status of group array after considering locus:\n";
		print LOG "\t", join("\t", @loci), "\n";
		
		if ($i == scalar(@{$pos_by_lg{$lg}} - 1)) { # The last locus on the comparison lg
			print LOG "End of map LG\n";
			print LOG "Adding current group to blocks\n";
			push $blocks{$lg}, [ @loci ];
			print LOG "Resetting the block array\n";
			@loci = ();
		}
	}
}

#print OUT Dumper(\%blocks);

print OUT join("\t", 'Locus', 'MAP_LG', 'MAP_POS', 'COMP_CHR', 'COMP_POS'), "\n";
foreach my $lg (keys %blocks) {
	foreach my $block (@{$blocks{$lg}}) {
		foreach my $locus (@{$block}) {
			print OUT join("\t", $locus, $map_info{$locus}[0], $map_info{$locus}[1], $comp_info{$locus}[0], $comp_info{$locus}[1]), "\n";
		}
		print OUT "---------------------------------------\n"
	}
	#print "\n";
}

# Produce some summary statistics for the run

print SUM join("\t", 'LOCI', 'NUM_LOCI', 'COMP_CHR', 'COMP_START', 'COMP_END', 'COMP_SIZE', 'MAP_LG', 'MAP_START', 'MAP_END', 'MAP_SIZE'), "\n";
print TAB join("\t", 'START_LOCUS', 'END_LOCUS', 'NUM_LOCI', 'COMP_CHR', 'COMP_START', 'COMP_END', 'COMP_SIZE', 'MAP_LG', 'MAP_START', 'MAP_END', 'MAP_SIZE'), "\n";

my %block_stats;
my @block_size;
my @comp_block_size;
my @num_loci_per_block;
foreach my $lg (keys %blocks) {
	foreach my $block (@{$blocks{$lg}}) {
		next if scalar(@{$block}) == 1;
		my $size = calc_block_size_map($block);
		my $comp_size = calc_block_size_comp($block);
		$block_stats{'num_blocks'}++;
		push @num_loci_per_block, scalar(@{$block});
		push @block_size, $size;
		push @comp_block_size, $comp_size;
		my $map_lg;
		foreach my $group (keys %lgs) {
			if (defined $lgs{$group}{$block->[0]}) {
				$map_lg = $group;
				$map_lg =~ s/LG//;
				last;
			}
		}
		
		print SUM join("\t", join(',', @{$block}), scalar(@{$block}), $comp_info{$block->[0]}[0], $comp_info{$block->[0]}[2], $comp_info{$block->[-1]}[2], $comp_size, $map_lg, $map_info{$block->[0]}[2], $map_info{$block->[-1]}[2], sprintf("%.2f", $size)), "\n";
		print TAB join("\t", $block->[0], $block->[-1], scalar(@{$block}), $comp_info{$block->[0]}[0], $comp_info{$block->[0]}[2], $comp_info{$block->[-1]}[2], $comp_size, $map_lg, $map_info{$block->[0]}[2], $map_info{$block->[-1]}[2], sprintf("%.2f", $size)), "\n";
		#print SUM "Block size: " . sprintf("%.2f", $size) . " cM\n";
		#print SUM "Comp block size: " . sprintf("%.2f", $comp_size) . " bp\n";
	}
	#print "\n";
}

my $size_stat = Statistics::Descriptive::Full->new();
$size_stat->add_data(@block_size);
#my $mean = $size_stat->mean();
#print $mean, "\n";

my $comp_size_stat = Statistics::Descriptive::Full->new();
$comp_size_stat->add_data(@comp_block_size);
#$comp_size_stat->mean();
#print $mean, "\n";

my $num_per_block_stat = Statistics::Descriptive::Full->new();
$num_per_block_stat->add_data(@num_loci_per_block);
#$mean = $num_per_block_stat->mean();
#print $mean, "\n";


print STAT "Total Loci in Blocks:\t", $num_per_block_stat->sum(), "\n";
print STAT "Total Blocks:\t", scalar(@num_loci_per_block), "\n";
print STAT "Mean loci per block:\t", $num_per_block_stat->mean(), "\n";
print STAT "Min loci per block:\t", $num_per_block_stat->min(), "\n";
print STAT "Max loci per block:\t", $num_per_block_stat->max(), "\n";
print STAT "Total Comp Block Size (bp):\t", $comp_size_stat->sum(), "\n";
print STAT "Mean Comp Block Size (bp):\t", $comp_size_stat->mean(), "\n";
print STAT "Min Comp Block Size (bp):\t", $comp_size_stat->min(), "\n";
print STAT "Max Comp Block Size (bp):\t", $comp_size_stat->max(), "\n";
print STAT "Total Map Block Size (cM):\t", $size_stat->sum(), "\n";
print STAT "Mean Map Block Size (cM):\t", $size_stat->mean(), "\n";
print STAT "Min Map Block Size (cM):\t", $size_stat->min(), "\n";
print STAT "Max Map Block Size (cM):\t", $size_stat->max(), "\n";


### Subroutines ###

sub get_dir {
	my $pos = shift;
	my $prev_pos = shift;
	my $diff = abs($pos - $prev_pos);
	if ($diff > 0.05) {
		if ($pos > $prev_pos) {
			return 1;
		} elsif ($pos < $prev_pos) {
			return -1;
		}
	} else {
		return 0;
	}
}

# sub get_dir {
	# my $pos = shift;
	# my $prev_pos = shift;
	# if ($pos > $prev_pos) {
		# return '+';
	# } elsif ($pos < $prev_pos) {
		# return '-';
	# } elsif ($pos == $prev_pos) {
		# return '=';
	# }
# }

sub same_dir {
	my $dir = shift;
	my $prev_dir = shift;
	if ($dir == $prev_dir || $prev_dir == 0 || $dir == 0) {
		return '1';
	} else {
	
		return '0';
	}
}

sub same_lg {
	my $map_lg = shift;
	my $prev_map_lg = shift;
	if ($map_lg eq $prev_map_lg) {
		return '1';
	} else {
		return '0';
	}
}

sub adjacent {
	my $num = shift;
	my $pos = shift;
	my $prev_pos = shift;
	my $lg = shift;
	for (my $z; $z < scalar(@{$pos_by_lg{$lg}}); $z++) {
		next if ($z == $num || $z == $num - 1); # don't consider the current locus or previous locus
		next if $pos_by_lg{$lg}[$z][2] ne $pos_by_lg{$lg}[$num][2]; # don't consider loci not on the same map linkage group
		if ($pos < $prev_pos) {
			if ($pos_by_lg{$lg}[$z][3] > $pos && $pos_by_lg{$lg}[$z][3] < $prev_pos) {
				if (abs($pos_by_lg{$lg}[$z][3] - $pos) <= 0.05 || abs($pos_by_lg{$lg}[$z][3] - $prev_pos) <= 0.05) {
					next;
				}
				return 0;
			}
		} elsif ($pos > $prev_pos) {
			if ($pos_by_lg{$lg}[$z][3] < $pos && $pos_by_lg{$lg}[$z][3] > $prev_pos) {
				if (abs($pos_by_lg{$lg}[$z][3] - $pos) <= 0.05 || abs($pos_by_lg{$lg}[$z][3] - $prev_pos) <= 0.05) {
					next;
				}
				return 0;
			}
		}
		
	}
	return 1;
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

sub calc_block_size_map {
	my @group = @{$_[0]};
	my $small = 500;
	my $large = 0;
	my $left;
	my $right;
	foreach my $locus (@group) {
		my $pos = $map_info{$locus}[2];
		if ($pos < $small) {
			$small = $pos;
			$left = $locus;
		}
		if ($pos > $large) {
			$large = $pos;
			$right = $locus;
		}
	}
	my $size = $large - $small;
	print OUT "Size: $large - $small, Left: $left ($small), Right: $right ($large)\n";
	return $size;
		
}

sub calc_block_size_comp {
	my @group = @{$_[0]};
	my $small = 100000000000;
	my $large = 0;
	my $left;
	my $right;
	foreach my $locus (@group) {
		my $pos = $comp_info{$locus}[2];
		if ($pos < $small) {
			$small = $pos;
			$left = $locus;
		}
		if ($pos > $large) {
			$large = $pos;
			$right = $locus;
		}
	}
	my $size = $large - $small;
	print OUT "Size: $large - $small, Left: $left ($small), Right: $right ($large)\n";
	return $size;
		
}

__END__

=head1 NAME

id_syntenic_blocks.pl

=head1 SYNOPSIS

id_syntenic_blocks.pl

Options:
     -i     infile
     -m     mapfile
     -g		genome

=head1 OPTIONS

=over 8

=item B<-i>

'.csv' formatted infile containing synteny info

=item B<-m>

'.map' file of linkage map

=item B<-g>

Genome of the comparison species - needs to be in the infile. Currently supports 'onil', 'gacu', 'tnig', 'tsub'

=back

=head1 DESCRIPTION

B<id_syntenic_blocks.pl> uses synteny data to identify blocks of shared synteny between red drum and a comparison species

=cut
