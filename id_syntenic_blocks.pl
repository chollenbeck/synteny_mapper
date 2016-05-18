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
my $error_level = 0.025;

GetOptions(		'infile|i=s' => \$infile,
				'mapfile|m=s' => \$mapfile,
				'genome|g=s' => \$genome,
				'error_level|e=s' => \$error_level,

		);


if ($error_level > 1 || $error_level < 0) {
	die "Error level needs to be between 0 and 1\n";
}

open(CSV, "<", $infile) or die $!;
open(DUMP, ">", 'dumper.out');
open(SUM, ">", $genome . '.synteny.sum');
open(LOG, ">", $genome . '.synteny.log');
open(OUT, ">", $genome . '.synteny.blocks');
open(STAT, ">", $genome . '.synteny.stats');
open(TAB, ">", $genome . '.synteny.tab');

open(DUMP, ">", 'dumper.out') or die $!;

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
	$map_info{$fields[0]} = [$fields[1], sprintf('%.2f', $fields[2]),  $fields[3]];
	$comp_info{$fields[0]} = [$fields[$chr_field], sprintf('%.2f', $fields[$prop_field]), $fields[$pos_field]];
}
close CSV;

my %lgs = read_map($mapfile);

my %map_lgs;
my @map_lgs;
my %pos_by_lg;

foreach my $locus (sort { ($map_info{$a}[0] <=> $map_info{$b}[0]) || ($map_info{$a}[1] <=> $map_info{$b}[1]) || ($comp_info{$a}[1] <=> $comp_info{$b}[1])} keys %comp_info) {
	next unless $comp_info{$locus}[0] && $map_info{$locus}[0] =~ /\d/;
	push @map_lgs, $map_info{$locus}[0] unless $map_lgs{$map_info{$locus}[0]};
	$map_lgs{$map_info{$locus}[0]}++;
	if (! $pos_by_lg{$map_info{$locus}[0]}) {
		$pos_by_lg{$map_info{$locus}[0]} = [ [ $locus, $map_info{$locus}[1], $comp_info{$locus}[0], $comp_info{$locus}[1] ] ];
	} else {
		push @{$pos_by_lg{$map_info{$locus}[0]}}, [ $locus, $map_info{$locus}[1], $comp_info{$locus}[0], $comp_info{$locus}[1] ];
	}
}

my %comp_chroms;
my @comp_chroms;
my %pos_by_chrom;

foreach my $locus (sort { ($comp_info{$a}[0] <=> $comp_info{$b}[0]) || ($comp_info{$a}[1] <=> $comp_info{$b}[1])} keys %map_info) {
	next unless $comp_info{$locus}[0] && $map_info{$locus}[0] =~ /\d/;
	push @comp_chroms, $comp_info{$locus}[0] unless $comp_chroms{$comp_info{$locus}[0]};
	$comp_chroms{$comp_info{$locus}[0]}++;
	if (! $pos_by_chrom{$comp_info{$locus}[0]}) {
		$pos_by_chrom{$comp_info{$locus}[0]} = [ [ $locus, $comp_info{$locus}[1], $map_info{$locus}[0], $map_info{$locus}[1] ] ];
	} else {
		push @{$pos_by_chrom{$comp_info{$locus}[0]}}, [ $locus, $comp_info{$locus}[1], $map_info{$locus}[0], $map_info{$locus}[1] ];
	}
}
#print join("\n", @comp_lgs), "\n";

#print LOG LOGer(\%pos_by_lg);


print "\tForming initial blocks\n";

# Step 1: Form initial blocks
my %blocks;
my %removed;
foreach my $lg (@map_lgs) {
	print LOG "-----Current map LG: $lg-----\n";
	$blocks{$lg} = [];
	my @loci;
	my $pos;
	my $prev_dir;
	my $map_lg;
	for (my $i = 0; $i < scalar(@{$pos_by_lg{$lg}}); $i++) {
		my $locus = $pos_by_lg{$lg}[$i][0];
		print LOG "Current locus: $locus\n";
		if (scalar(@loci) == 0) { # first locus
			print LOG "Array is empty, starting a new group\n";
			push @loci, $locus;
			#print join("\t", $locus, $lg, $pos_by_lg{$lg}[$i][3], $pos_by_lg{$lg}[$i][2]), "\n";

		} elsif (scalar(@loci) == 1) {
			print LOG "Group has one locus, checking to see if the current locus is syntenic\n";
			my $prev_locus = $pos_by_lg{$lg}[$i-1][0];
			my $same_lg = same_lg($pos_by_lg{$lg}[$i][2], $pos_by_lg{$lg}[$i-1][2]); # Checking to see if locus is on the same comparison chrom
			my $adjacent = adj($locus, $prev_locus, $pos_by_lg{$lg}[$i][3], $pos_by_lg{$lg}[$i-1][3], $pos_by_lg{$lg}[$i][2], \%removed, 0);
			if ($same_lg && $adjacent) {
				print LOG "Same LG:$same_lg, Adjacent:$adjacent\n";
				print LOG "Locus is syntenic - adding to group\n";
				push @loci, $locus;
				$prev_dir = get_dir($pos_by_lg{$lg}[$i][3], $pos_by_lg{$lg}[$i-1][3], 0);
				print LOG "Establishing a direction: $prev_dir\n";
			} else {
				print LOG "Same LG:$same_lg, Adjacent:$adjacent\n";
				print LOG "Locus is not syntenic - starting a new group\n";
				push $blocks{$lg}, [ @loci ];
				@loci = ($locus);
			}

		} elsif (scalar(@loci) > 1) {
			print LOG "Group has more than one locus, checking to see if the current locus is syntenic\n";
			my $prev_locus = $pos_by_lg{$lg}[$i-1][0];
			my $same_lg = same_lg($pos_by_lg{$lg}[$i][2], $pos_by_lg{$lg}[$i-1][2]);
			my $adjacent = adj($locus, $prev_locus, $pos_by_lg{$lg}[$i][3], $pos_by_lg{$lg}[$i-1][3], $pos_by_lg{$lg}[$i][2], \%removed, 0);
			my $dir = get_dir($pos_by_lg{$lg}[$i][3], $pos_by_lg{$lg}[$i-1][3], 0);
			my $same_dir = same_dir($dir, $prev_dir);
			if ($same_lg && $adjacent && $same_dir) {
				print LOG "Same LG:$same_lg, Adjacent:$adjacent, Same dir:$same_dir ($dir, $prev_dir)\n";
				print LOG "Locus is syntenic - adding to group\n";
				push @loci, $locus;
				$prev_dir = $dir;
			} else {
				print LOG "Same LG:$same_lg, Adjacent:$adjacent, Same dir:$same_dir\n";
				if ((abs($pos_by_lg{$lg}[$i][3] - $pos_by_lg{$lg}[$i][3]) <= $error_level) && $same_lg && $adjacent) {
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

print "\tInitial round: ", count_blocks(\%blocks), "\n";



# Step 2: Remove single loci that map to different chromosomes (errors/mulitcopy/small translocations)
print "\tRemoving single locus translocations", "\n";
my ($new_blocks, $removed) = remove_single_breakers(\%blocks, \%comp_info);
%blocks = %{$new_blocks};
%removed = %{$removed};




# Step 3: Iteratively merge blocks that can now be combined
print "\tMerging blocks", "\n";
my $prev_blocks;
my $i;
for ($i = 1; $i < 100; $i++) {
	my %new_blocks = merge(\%blocks, \%removed);
	print "\tRound $i: ", count_blocks(\%new_blocks), "\n";
	%blocks = %new_blocks;
	last if $prev_blocks == count_blocks(\%new_blocks);
	$prev_blocks = count_blocks(\%new_blocks);
}



# Step 4: Swap positions (in the array) of single blocks with the same map position
print "\tSwapping single locus blocks with same map position (error rate = 0)", "\n";
%blocks = swap_singles(\%blocks, 0);

#print DUMP Dumper(\%blocks);

# Step 5: Iteratively merge blocks that can now be combined
print "\tMerging blocks", "\n";
my $j;
for ($j = $i + 1; $j < 100; $j++) {
	my %new_blocks = merge(\%blocks, \%removed);
	print "\tRound $j: ", count_blocks(\%new_blocks), "\n";
	%blocks = %new_blocks;
	last if $prev_blocks == count_blocks(\%new_blocks);
	$prev_blocks = count_blocks(\%new_blocks);
}

# Step 6: Swap positions (in the array) of single blocks with the same map position (with some error rate)
print "\tSwapping single locus blocks with same map position (error rate = $error_level)", "\n";
%blocks = swap_singles(\%blocks, $error_level);

# Step 7: Iteratively merge blocks that can now be combined
print "\tMerging blocks", "\n";
my $k;
for ($k = $j + 1; $k < 100; $k++) {
	my %new_blocks = merge(\%blocks, \%removed);
	print "\tRound $j: ", count_blocks(\%new_blocks), "\n";
	%blocks = %new_blocks;
	last if $prev_blocks == count_blocks(\%new_blocks);
	$prev_blocks = count_blocks(\%new_blocks);
}



#print DUMP Dumper(\%blocks);

#my %overmerged = test_blocks(\%blocks, \%comp_info, \%map_info, $error_level);

#print Dumper(\%overmerged);

# Re-sort the comp positions relative to map positions and rebuild


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
	my $pos = $_[0];
	my $prev_pos = $_[1];
	my $error_rate = $_[2];

	my $diff = abs($pos - $prev_pos);
	if ($diff > $error_rate) {
		if ($pos > $prev_pos) {
			return 1;
		} elsif ($pos < $prev_pos) {
			return -1;
		}
	} else {
		return 0;
	}
}

sub same_dir {
	my $dir = shift;
	my $prev_dir = shift;
	if ($dir == $prev_dir || $prev_dir == 0 || $dir == 0) {
		return 1;
	} else {

		return 0;
	}
}

sub same_lg {
	my $map_lg = shift;
	my $prev_map_lg = shift;
	if ($map_lg eq $prev_map_lg) {
		return 1;
	} else {
		return 0;
	}
}

sub adj {
	my $locus = $_[0];
	my $prev_locus = $_[1];
	my $pos = $_[2];
	my $prev_pos = $_[3];
	my $chrom = $_[4];
	my %removed = %{$_[5]};
	my $error_level = $_[6];

	for (my $z; $z < scalar(@{$pos_by_chrom{$chrom}}); $z++) {

		if ($pos_by_chrom{$chrom}[$z][0] eq $locus || $pos_by_chrom{$chrom}[$z][0] eq $prev_locus) { # don't consider the current locus or previous locus
			next;
		}
		if ($removed{$pos_by_chrom{$chrom}[$z][0]}) { # skip the locus if it has been removed from consideration
			next;
		}
		if ($pos < $prev_pos) {
			if ($pos_by_chrom{$chrom}[$z][1] > $pos && $pos_by_chrom{$chrom}[$z][1] < $prev_pos) {
				#if (abs($pos_by_chrom{$chrom}[$z][1] - $pos) <= 0.02 || abs($pos_by_chrom{$chrom}[$z][1] - $prev_pos) <= 0.02) {
				#	print "Too small, skipping\n" if $locus eq 'Contig_1597';
				#	next;
				#}
				#print "This is what happened\n" if $locus eq "Soc1128";
				return 0;
			}
		} elsif ($pos > $prev_pos) {
			if ($pos_by_chrom{$chrom}[$z][1] < $pos && $pos_by_chrom{$chrom}[$z][1] > $prev_pos) {
				#if (abs($pos_by_chrom{$chrom}[$z][1] - $pos) <= 0.02 || abs($pos_by_chrom{$chrom}[$z][1] - $prev_pos) <= 0.02) {
				#	print "Too small, skipping\n" if $locus eq 'Contig_1597';
				#	next;
				#}
				return 0;
			}
		}

	}
	return 1;
}

sub remove_single_breakers {
	my %blocks = %{$_[0]};
	my %comp = %{$_[1]};

	my %removed;
	my %new_blocks;
	foreach my $lg (keys %blocks) {

		# Go through each block and check to see if it is a single locus
		my %remove;
		for (my $i = 1; $i < scalar(@{$blocks{$lg}} - 1); $i++) { # indexes block (skips the first and last locus

			# Skip the block if it has more than one locus
			if (scalar(@{$blocks{$lg}[$i]}) != 1) {
				next;
			}
			my $locus = $blocks{$lg}[$i][0];
			my $prev_locus = $blocks{$lg}[$i - 1][0];
			my $next_locus = $blocks{$lg}[$i + 1][0];

			my $chrom = $comp{$locus}[0];
			my $prev_chrom = $comp{$prev_locus}[0];
			my $next_chrom = $comp{$next_locus}[0];

			if ($chrom ne $prev_chrom && $chrom ne $next_chrom) {
				# The block is not on the same chromosome as the ones before and after; remove it
				$remove{$lg}{$i} = 1;
			}


		}

		$new_blocks{$lg} = [ ];
		for (my $i = 0; $i < scalar(@{$blocks{$lg}}); $i++) {
			if ($remove{$lg}{$i}) {
				$removed{$blocks{$lg}[$i][0]} = 1;
				next;
			}
			push @{$new_blocks{$lg}}, $blocks{$lg}[$i];
		}


	}

	return (\%new_blocks, \%removed);
}

sub merge {

	my %blocks = %{$_[0]};
	my %removed = %{$_[1]};

	my %new_blocks;
	foreach my $lg (keys %blocks) {
		my @merge_codes = ();
		my $code = 0;
		# Go through each block and check if they are mergable
		for (my $i = 0; $i < scalar(@{$blocks{$lg}} - 1); $i++) { # indexes block
			my $curr_block = $blocks{$lg}[$i];
			my $next_block = $blocks{$lg}[$i + 1];

			# Check to see if the two blocks are mergable

			my $merge = mergable($curr_block, $next_block, \%comp_info, \%map_info, \%pos_by_chrom, \%pos_by_lg, $lg, \%removed);

			if ($merge == 1) { # Mergable
				$merge_codes[$i] = $code;
				$merge_codes[$i + 1] = $code;
				#die;
			} else { # Not mergable
				$merge_codes[$i] = $code;
				$merge_codes[$i + 1] = $code + 1;
				$code++;
			}

		}

		# Merge blocks that are able to be merged

		for (my $i = 0; $i < scalar(@{$blocks{$lg}}); $i++) {
			if ( ! $new_blocks{$lg}[$merge_codes[$i]] ) {
				$new_blocks{$lg}[$merge_codes[$i]] = [ @{$blocks{$lg}[$i]} ];
			} else {
				push @{$new_blocks{$lg}[$merge_codes[$i]]}, @{$blocks{$lg}[$i]};
			}
		}


	}

	return %new_blocks;
}

sub mergable { # Get it?
	my @block1 = @{$_[0]};
	my @block2 = @{$_[1]};
	my %comp = %{$_[2]};
	my %map = %{$_[3]};
	my %pos_by_chrom = %{$_[4]};
	my %pos_by_lg = %{$_[5]};
	my $lg = $_[6];
	my %removed = %{$_[7]};

	my $error_level = 0.025;

	 # Groups will be mergable if:
	 # 1) there are no other loci in the region represented by both groups

	 # Make sure the blocks are on the same comparison chromosome:



	if ($comp{$block1[0]}[0] ne $comp{$block2[0]}[0]) {
		 return 0;
	}

	 print LOG "Checking to see if two blocks are mergable\n";

	my $chrom = $comp{$block1[0]}[0];

	# Combine the groups
	my @comb = (@block1, @block2);

	print LOG "Combined block:\n";
	print LOG Dumper(\@comb);
	my %combined;
	foreach my $locus (@comb) {
		$combined{$locus} = [ $map{$locus}[0], $map{$locus}[1], $comp{$locus}[0], $comp{$locus}[1] ]; # Map LG, Map Prop, Comp LG, Comp Prop
	}

	# Calculate the distance between groups


	# Sort by comparision position and get the range of the blocks in the comparison species
	my @merged;
	my %merged_loci;
	foreach my $locus (sort { $combined{$a}[3] <=> $combined{$b}[3] } keys %combined ) {
		push @merged, $combined{$locus};
		$merged_loci{$locus} = 1;
	}
	#print Dumper(\@merged);

	my $high = $merged[-1][3];
	my $low = $merged[0][3];
	#print join("\t", $low, $high), "\n";

	# Check to see if any other loci occupy the same range

	# Check to see if there are blocks that exist within the combined block
	for (my $z; $z < scalar(@{$pos_by_chrom{$chrom}}); $z++) {
		my $locus = $pos_by_chrom{$chrom}[$z][0];
		my $pos = $pos_by_chrom{$chrom}[$z][1];
		next if $merged_loci{$locus};
		next if $removed{$locus};
		if ($pos > $low && $pos < $high) {
			print LOG "Not mergable: Other blocks within combined block\n";
			return 0;
		}

	}
	print LOG "No other blocks within combined block\n";


	# Check to see if the combined block has other problems
	my $prev_dir = 0;
	for (my $i = 1; $i < scalar(@comb); $i++) {

		my $locus = $comb[$i];
		my $prev_locus = $comb[$i - 1];
		my $next_locus = $comb[$i + 1] if $i < scalar(@comb);
		my $pos = $combined{$locus}[3];
		my $prev_pos = $combined{$prev_locus}[3];
		my $next_pos = $combined{$next_locus}[3] if $i < scalar(@comb);

		my $chrom = $combined{$locus}[2];
		my $adjacent = adj($locus, $prev_locus, $pos, $prev_pos, $chrom, \%removed, $error_level);
		my $dir = get_dir($pos, $prev_pos, $error_level);
		my $same_dir = same_dir($dir, $prev_dir);

		if ($adjacent && $same_dir) {
			$prev_dir = $dir;
			next;
		} elsif (! $adjacent || ! $same_dir) {
			if (abs($pos - $prev_pos) <= $error_level) {
				$prev_dir = $dir;
				next;
			} elsif ($next_locus) { # Try to see if the problem can be resolved by swapping the locus with the next locus
				my $adj_2 = adj($prev_locus, $next_locus, $prev_pos, $next_pos, $chrom, \%removed, $error_level);
				print LOG "Not adjacent. Trying to swap locus with next locus...\n";
				if ($adj_2) {
					if ((abs($pos - $next_pos) <= $error_level) && (abs($combined{$locus}[1] - $combined{$next_locus}[1]) <= $error_level) ) {
						print LOG "Adjacent and within error rate. Swap success\n";
						$prev_dir = 0;
						$i += 1;
						next;
					} else {
						print LOG "Adjacent but outside of error rate. Swap fail - not mergable\n";
						return 0;
					}
				} else {
					print LOG "Next locus not adjacent. Swap fail - not mergable\n";
					return 0;
				}
			} else {
				print LOG "No locus to swap - not mergable\n";
				return 0;
			}
		} else {
			print LOG "Something strange happened.\n";
			return 0;
		}



		#if (! $adjacent ||
	}


	return 1;
}

sub test_blocks {
	my %blocks = %{$_[0]};
	my %comp = %{$_[1]};
	my %map = %{$_[2]};
	my $error_level = $_[3];

	my %new_blocks;
	my %overmerged;
	foreach my $lg (keys %blocks) {

		for (my $i = 0; $i < scalar(@{$blocks{$lg}}); $i++) { # indexes block

			# Skip single locus blocks
			if (scalar(@{$blocks{$lg}[$i]}) == 1) {
				next;
			}

			my %info;
			for (my $j = 0; $j < scalar(@{$blocks{$lg}[$i]}); $j++) {
				my $locus = $blocks{$lg}[$i][$j];
				$info{$locus} = [ $map{$locus}[0], $map{$locus}[1], $comp{$locus}[0], $comp{$locus}[1] ];
			}

			my $num_loci = scalar(keys %info);

			#print Dumper(\%info);


			# Iterate through each locus, checking for overmerged blocks
			my $over1;
			my $over2;

			#print "Going forward:\n";
			my $prev_dir = 0;
			for (my $j = 1; $j < scalar(@{$blocks{$lg}[$i]}); $j++) {
				my $locus = $blocks{$lg}[$i][$j];
				my $prev_locus = $blocks{$lg}[$i][$j - 1];
				my $next_locus = $blocks{$lg}[$i][$j + 1] if $j < $num_loci - 1;


				my $pos = $info{$locus}[3];
				my $prev_pos = $info{$prev_locus}[3];
				my $next_pos = $info{$next_locus}[3] if $j < $num_loci - 1;

				my $chrom = $info{$locus}[2];

				#print "Considering locus: $locus\n";
				#print join('--', $locus, $prev_locus, $next_locus, $pos, $prev_pos, $next_pos, $chrom), "\n";




				my $adjacent = adj($locus, $prev_locus, $pos, $prev_pos, $chrom, \%removed, $error_level);
				my $dir = get_dir($pos, $prev_pos, $error_level);
				my $same_dir = same_dir($dir, $prev_dir);

				if (! $adjacent || ! $same_dir) {
					if (abs($pos - $prev_pos) <= $error_level) {
						next;
					} else {
						$over1 = 1;
						last;
					}
				}
			}

			#print Dumper(\%info);


			# Try the other direction
			#print "Going backwards:\n";
			my $prev_dir = 0;
			for (my $j = -2; $j >= (-1 * $num_loci); $j--) {
				my $locus = $blocks{$lg}[$i][$j];
				my $prev_locus = $blocks{$lg}[$i][$j + 1];
				my $next_locus = $blocks{$lg}[$i][$j - 1] if $j < -1 * $num_loci;


				my $pos = $info{$locus}[3];
				my $prev_pos = $info{$prev_locus}[3];
				my $next_pos = $info{$next_locus}[3] if $j < -1 * $num_loci;
				my $chrom = $info{$locus}[2];

				#print "Considering locus: $locus\n";
				#print Dumper(\%info);

				my $adjacent = adj($locus, $prev_locus, $pos, $prev_pos, $chrom, \%removed, $error_level);
				#my $adjacent = 1;
				my $dir = get_dir($pos, $prev_pos, 0);
				my $same_dir = same_dir($dir, $prev_dir);

				if (! $adjacent || ! $same_dir) {
					if (abs($pos - $prev_pos) <= $error_level) {
						next;
					} else {
						$over2 = 1;
						last;
					}
				}
			}
			if ($over1 && $over2) {
				$overmerged{$blocks{$lg}[$i][0]} = 1;
			}

		}

	}

	return %overmerged;
}

sub swap_singles {
	my %blocks = %{$_[0]};
	my $error_level = $_[1];


	my %new_blocks;
	foreach my $lg (keys %blocks) {
		# Look for single locus blocks that have the same map position
		my %positions;
		for (my $i = 0; $i < scalar(@{$blocks{$lg}}); $i++) { # indexes block
			if (scalar(@{$blocks{$lg}[$i]}) != 1) {
				$positions{$lg}{$i} = $i;
				next;
			}

			#print $i, "\t", scalar(@{$blocks{11}}), "\n" if $lg == 11;


			my $locus = $blocks{$lg}[$i][0];
			my $num_matches = 0;
			for (my $j = 1; $j < scalar(@{$blocks{$lg}}); $j++) {
				last unless $blocks{$lg}[$i + $j];
				my $locus2 = $blocks{$lg}[$i + $j][0];
				if (scalar(@{$blocks{$lg}[$i + $j]}) == 1 && (abs($map_info{$locus}[1] - $map_info{$locus2}[1]) <= $error_level) && ($comp_info{$locus}[0] eq $comp_info{$locus2}[0])) {
					$num_matches++;
				} else {
					last;
				}
			}



			if ($num_matches > 0) {
				# Put them in reverse order
				my $num = $i + $num_matches;
				for (my $j = $i; $j <= $i + $num_matches; $j++) {
					$positions{$lg}{$j} = $num;
					$num--;
				}
				$i += $num_matches;

			} else { # no matches
				$positions{$lg}{$i} = $i;
			}
		}

		#Swap the block positions
		#$new_blocks{$lg} = [ ];
		for (my $i = 0; $i < scalar(@{$blocks{$lg}}); $i++) { # indexes block
			$new_blocks{$lg}[$positions{$lg}{$i}] = $blocks{$lg}[$i];
		}



	}


	return %new_blocks;
}

sub count_blocks {
	my %blocks = %{$_[0]};

	my $num_blocks = 0;
	foreach my $lg (keys %blocks) {
		foreach my $block (@{$blocks{$lg}}) {
			next unless scalar(@$block) == 1;
			$num_blocks++;
		}
	}

	return $num_blocks;

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
