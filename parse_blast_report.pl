#!/usr/bin/perl
use strict;
use List::MoreUtils qw/ uniq /;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

pod2usage(-verbose => 1) if @ARGV == 0;

my @files;
my $outfile = 'test.out';
my $mapfile;
my $locfile;


GetOptions(		'files|f=s{1,}' => \@files,
				'outfile|o=s' => \$outfile,
				'mapfile|m=s' => \$mapfile,
				'locfile|l=s' => \$locfile,

	);



open (DUMP, ">", 'dumper.out');


# Read in the locus data

my %loci;
open(LOC, "<", $locfile) or die $!;
while(<LOC>) {
	next if $_ =~ /^\s/;
	chomp;
	$loci{$_} = 1;
}

# Read in the map data

my %map = read_map($mapfile);
my %map_mod = convert_lg_pos(\%map);


my %map_loci;
foreach my $lg (keys %map) {
	foreach my $locus (keys $map{$lg}) {
		$map_loci{$locus}++;
	}
}

my @chr_lengths = read_chr_lengths();

my %data;
my %stats;
foreach my $file (@files) {

	my $genome;
	if ($file =~ /^(\w+?)\./) {
		$genome = $1;
	}

	$stats{'total_map_loci'} = scalar(keys %map_loci);

	open(IN, "<", $file);

	my %chroms;
	my %positions;
	while(<IN>) {
		# Get the locus name from the line
		my $locus;
		if ($_ =~ /^gi\|(\d+)\|gb/) {
			$locus = $1;
		} elsif ($_ =~ /^(\S+)\s/) {
			$locus = $1;
		}


		if (! $chroms{$locus}) {
			$chroms{$locus} = [];
		}
		my $chrom;
		my $pos;
		if ($genome eq 'drer' || $genome eq 'tnig') {
			if ($_ =~ /^\S+\s(\S+)/) {
				$chrom = $1;
			}
			if ($_ =~ /(\d+)\s\d+$/) {
				$pos = $1;
			}
		}
		if ($genome eq 'gacu') {
			if ($_ =~ /^\S+\s(\S+)/) {
				$chrom = $1;
			}
			if ($_ =~ /(\d+)\s\d+$/) {
				$pos = $1;
			}
			next if $chrom eq 'chrUn';
			$chrom = convert_gacu($chrom);

		}
		if ($genome eq 'onil') {
			if ($_ =~ /LG(\S+)?,/) {
				$chrom = $1;
				$chrom =~ s/-\d\d//;
			}
			if ($_ =~ /(\d+)\s\d+$/) {
				$pos = $1;
			}
		}
		if ($genome eq 'trub') {
			if ($_ =~ /Chr_(\w+)/) {
				$chrom = $1;
			}
			if ($_ =~ /(\d+)\s\d+$/) {
				$pos = $1;
			}
		}
		if ($genome eq 'lcal') {
			if ($_ =~ /LG(\w+)/) {
				$chrom = $1;
			}

			if ($chrom == '16_LG22') {
				$chrom = '16-22';
			}
			if ($_ =~ /(\d+)\s\d+$/) {
				$pos = $1;
			}

		}
		if ($genome eq 'dlab') {
			next if $_ =~ /UN/;
			if ($_ =~ /LG(\w+)/) {
				$chrom = $1;
			}
			if ($chrom == '22') {
				$chrom = '22-25';
			}
			if ($chrom == '18') {
				$chrom = '18-21';
			}
			if ($_ =~ /(\d+)\s\d+$/) {
				$pos = $1;
			}
		}
		my $prop = convert_pos($chrom, $pos, $genome);
		push @{$chroms{$locus}}, $chrom;
		if (! $data{$genome}{$locus}) {
			$data{$genome}{$locus} = [[$chrom, $prop, $pos]];
		} else {
			push @{$data{$genome}{$locus}}, [$chrom, $prop, $pos];
		}

	}
	close IN;

	foreach my $lg (keys %map) {
		foreach my $locus (keys %{$map{$lg}}) {

			next unless $map_loci{$locus} == 1;
			if (! $data{$genome}{$locus}) { # if there is no hit for a locus...
				$stats{$genome}{'map_no_hit'}++;
				next;
			} else { # there are hits
				my @unique = uniq(@{$chroms{$locus}});
				if (scalar(@unique == 1)) { # if the locus hit a single chromosome
					print "Single\n" if $locus eq 'Contig_2132';
					$stats{$genome}{'map_hit'}++;
					$stats{$genome}{'map_single_hit'}++;
					# use the information from the first hit in the file
					$data{$genome}{$locus} = $data{$genome}{$locus}[0];
				} else { # the locus hit multiple chromosomes
					print "Mult\n" if $locus eq 'Contig_2132';
					$stats{$genome}{'map_hit'}++;
					$stats{$genome}{'map_multiple_hit'}++;

					# code the locus as no hit
					$data{$genome}{$locus} = ['', '', ''];
				}
			}
		}
	}
	#write_summary_file(\%stats, $genome);
}

#print DUMP Dumper(\%map_loci);
#print DUMP Dumper(\%data);
#print DUMP Dumper(\%stats);


# Open the outfile and write the combined data, separated by commas

open(OUT, ">", $outfile);

print OUT 'locus,map_chr,map_prop,map_pos,';
foreach my $genome (keys %data) {
	print OUT join(',', $genome . '_chr', $genome . '_prop', $genome . '_pos'), ',';
}
print OUT "\n";

foreach my $locus (keys %map_loci) {
	next unless $loci{$locus};
	print OUT $locus, ',';
	foreach my $lg (keys %map) {
		if (defined $map{$lg}{$locus}) {
			print OUT join(',', $lg, $map_mod{$lg}{$locus}, $map{$lg}{$locus}), ',';
			last;
		}
	}
	foreach my $genome (keys %data) {
		if (! $data{$genome}{$locus}) { # If there is no match for the species
			print OUT ',,,'; # leave the columns blank
		} else { # Otherwise, if there is a match for the species
			print OUT join(',', @{$data{$genome}{$locus}}), ','; # Record the info
		}
	}
	print OUT "\n";
}


sub convert_gacu {
	my $chrom = shift;
	$chrom =~ s/chr//;
	my %gacu_chrom = (
	'I' => 1,
	'II' => 2,
	'III' => 3,
	'IV' => 4,
	'V' => 5,
	'VI' => 6,
	'VII' => 7,
	'VIII' => 8,
	'IX' => 9,
	'X' => 10,
	'XI' => 11,
	'XII' => 12,
	'XIII' => 13,
	'XIV' => 14,
	'XV' => 15,
	'XVI' => 16,
	'XVII' => 17,
	'XVIII' => 18,
	'XIX' => 19,
	'XX' => 20,
	'XXI' => 21,
	);

	my $new_chrom = $gacu_chrom{$chrom};
	return $new_chrom;
}

sub read_chr_lengths {

	my @genomes = ();
	my $counter = 0;
	while(<DATA>) {
		if ($_ =~ /^>(\w+)/) {
			push @genomes, { 'name' => $1 };
			next;
		}
		if ($_ =~ /^\s/) {
			$counter++;
			next;
		}
		my ($chrom, $length) = split;
		if (${$genomes[$counter]}{'name'} eq 'gacu') {
			$chrom = convert_gacu($chrom);
		}
		${$genomes[$counter]}{$chrom} = $length;
	}
	#print DUMP Dumper(\@genomes);
	return @genomes;
}

sub convert_pos {
	my $chrom = shift;
	my $pos = shift;
	my $genome = shift;
	my %lengths;
	foreach my $org (@chr_lengths) {
		if ($$org{'name'} eq $genome) {
			%lengths = %{$org};
		}
	}
	my $total = $lengths{$chrom};
	my $prop = $pos/$total;
	# print $prop, "\n";
	return sprintf("%.3f", $prop);
}

sub write_summary_file {
	my %stats = %{$_[0]};
	my $genome = $_[1];
	open(SUM, ">", $genome . '_blast_summary.out');
	print SUM Dumper(\%stats);

	print SUM "Summary of $genome BLAST:\n";
	print SUM "Total_loci: ", $stats{'total_loci'}, "\n";
	print SUM "Total Map Loci: ", $stats{'total_map_loci'}, " ($stats{'total_map_neutral'} neutral, $stats{'total_map_est'} EST)", "\n";
	print SUM "\tMap Non-Hits: ", $stats{'map_no_hit'}, "\n";
	print SUM "\tMap Hits: ", $stats{'map_hit'}, "\n";
	print SUM "\t\tMap Hits (ESTs): ", $stats{'map_hit_est'}, "\n";
	print SUM "\t\t\tMultiple hits (ESTs): ", $stats{'map_multiple_hit_est'}, "\n";
	print SUM "\t\tMap Hits (neutral): ", $stats{'map_hit_neutral'}, "\n";
	print SUM "\t\t\tMultiple hits (neutral): ", $stats{'map_multiple_hit_neutral'}, "\n";
	print SUM "Total comparisons: ", $stats{'map_hit'} - $stats{'map_multiple_hit_neutral'} - $stats{'map_multiple_hit_est'}, "\n";

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

sub convert_lg_pos {
	my %lgs= %{$_[0]};
	my %lg_mod;
	foreach my $lg (keys %lgs) {
		my @positions;
		foreach my $locus (sort { ($lgs{$lg}{$a} <=> $lgs{$lg}{$b})} keys $lgs{$lg}) {
			push @positions, $lgs{$lg}{$locus};
		}
		my $length = $positions[-1];
		foreach my $locus (sort { ($lgs{$lg}{$a} <=> $lgs{$lg}{$b})} keys $lgs{$lg}) {
			my $new_pos = sprintf("%.3f", $lgs{$lg}{$locus} / $length);
			$lg_mod{$lg}{$locus} = $new_pos;
		}
	}
	return %lg_mod;
}


__END__
>drer
10	46591166
11	46661319
12	50697278
13	54093808
14	53733891
15	47442429
16	58780683
17	53984731
18	49877488
19	50254551
1	60348388
20	55952140
21	44544065
22	42261000
23	46386876
24	43947580
25	38499472
2	60300536
3	63268876
4	62094675
5	75682077
6	59938731
7	77276063
8	56184765
9	58232459
MT	16596

>gacu
I	28185914
II	23295652
III	16798506
IV	32632948
IX	20249479
M	15742
Un	62550211
V	12251397
VI	17083675
VII	27937443
VIII	19368704
X	15657440
XI	16706052
XII	18401067
XIII	20083130
XIV	15246461
XIX	20240660
XV	16198764
XVI	18115788
XVII	14603141
XVIII	16282716
XX	19732071
XXI	11717487

>onil
10	17092887
11	33447472
12	34679706
13	32787261
14	34191023
15	26684556
16	34890008
17	31749960
18	26198306
19	27159252
1	31194787
20	31470686
22	26410405
23	20779993
2	25048291
3	19325363
4	28679955
5	37389089
6	36725243
7	51042256
8	29447820
9	20956653

>tnig
10	13272281
11	11954808
12	12622881
13	13302670
14	10246949
15	7320470
16	9031048
17	12136232
18	11077504
19	7272499
1	22981688
20	3798727
21	5834722
2	21591555
3	15489435
4	9874776
5	13390619
6	7024381
7	11693588
8	10512681
9	10554956
MT	16462

>trub
10	9126174
11	11238946
12	10661439
13	16952093
14	12328026
15	11516971
16	10310320
17	11184067
18	8053132
19	14185242
1	23260604
20	13411851
21	14900632
22	11609757
2	11407286
3	13093193
4	13818948
5	12184708
6	9412079
7	14304132
8	14690302
9	13922460

>dlab
10     24053896
11     26215200
12     23234908
13     27622843
14     28395245
15     25566315
16     25821574
17     22893611
18-21  16453912
19     23202889
1A     29036788
1B     17995764
2      26218737
20     28398722
22-25  26439989
24     13918872
3      13451259
4      27569954
5      32612477
6      28185029
7      28480941
8      23067096
9      22362490
x      17779090

>lcal
23	18168282
4	25538952
12	27842965
24	19811778
8	25919959
3	23499962
17	27673719
1	25703306
21	28676982
5	28963731
10	27937307
11	23293155
19	24524913
13	27244013
7_1	23258384
7_2	13910880
2	30394535
6	27924252
15	30776907
18	19193443
9	22990584
20	23753645
16-22	25848596
14	14073782


=head1 NAME

parse_map_blast_report.pl

=head1 SYNOPSIS

parse_map_blast_report.pl

Options:
     -f     files
     -m     mapfile
     -l		locusfile
     -o		outfile

=head1 OPTIONS

=over 8

=item B<-f>

BLAST reports to be parsed and combined (may be multiple files)

=item B<-m>

'.map' file of linkage map

=item B<-l>

file of loci to be listed in the outfile

=item B<-o>

output file written in 'csv' format

=back

=head1 DESCRIPTION

B<parse_map_blast_report.pl> parses BLAST reports to generate a file containing BLAST hit information across multiple species

=cut
