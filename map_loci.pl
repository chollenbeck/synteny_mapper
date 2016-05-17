#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;
use List::MoreUtils qw/ uniq /;

pod2usage(-verbose => 1) if @ARGV == 0;

my $infile; # FASTA format
my $synfile;
my $blockfile;
my $genome;
my $blast_db;
my $db_dir = '/home/blast/db/';
my $min_length = 50;
my $argfile;


GetOptions(		'infile|i=s' => \$infile,
				'synfile|s=s' => \$synfile,
				'blockfile|l=s' => \$blockfile,
				'genome|g=s' => \$genome,
				'blast_db|b=s' => \$blast_db,
				'db_dir|z=s' => \$db_dir,
				'min_length|m=s' => \$min_length,
				'argfile|a=s' => \$argfile,

	);

open(DUMP, ">", 'dumper.out') or die $!;

# BLAST loci

# Open the file with BLAST arguments and read in the parameters

open(ARG, "<", $argfile) or die $!;
my @args;
while(<ARG>) {
	next if $_ =~ /^#/;  # Skip comment lines
	chomp;
	my ($flag, $option) = split;
	push @args, "-$flag";
	push @args, $option;
}

# Read in the FASTA file
my $seq_in = Bio::SeqIO->new(-file => $infile,
							-format => 'fasta');

my @seqs;
while (my $seq = $seq_in->next_seq) {
	push @seqs, $seq;
}

# Create the BLAST factory object
my $fac = Bio::Tools::Run::StandAloneBlastPlus->new( -db_name => $blast_db, -db_dir => $db_dir );

open(TAB, ">", $genome . '.synteny.hits.tab') or die $!;

my %hits;
foreach my $seq (@seqs) {
	my $result = $fac->blastn( -query => $seq,
                        -method_args => \@args );

		while (my $hit = $result->next_hit ) {
			while (my $hsp = $hit->next_hsp ) {
				if ($hsp->length('total') >= $min_length ) {
					print TAB join("\t", $result->query_name, $hit->name, $hit->description, $hit->significance, $hsp->start('hit'), $hsp->end('hit')), "\n";
				}
			}
		}
}
$fac->cleanup;

close TAB;


# Record hit data
open(HIT, "<", $genome . '.synteny.hits.tab') or die $!;

my @chr_lengths = read_chr_lengths();

my %chroms;
my %data;
while(<HIT>) {
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
	if (! $data{$locus}) {
		$data{$locus} = [[$chrom, $prop, $pos]];
	} else {
		push @{$data{$locus}}, [$chrom, $prop, $pos];
	}

}
close HIT;

print DUMP Dumper(\%data);

my %single_hits;
my %stats;
# Check to see if the locus has single or multiple hits and record some stats
foreach my $locus (keys %data) {

	if (! $data{$locus}) { # if there is no hit for a locus...
		$stats{'no_hit'}++;
		next;
	} else { # there are hits
		my @unique = uniq(@{$chroms{$locus}});
		if (scalar(@unique == 1)) { # if the locus hit a single chromosome
			$stats{'hit'}++;
			$stats{'single_hit'}++;
			# use the information from the first hit in the file
			$single_hits{$locus} = $data{$locus}[0];
		} else { # the locus hit multiple chromosomes
			$stats{'hit'}++;
			$stats{'multiple_hit'}++;
		}
	}
}

print DUMP Dumper(\%single_hits);

# Identify if loci are in syntenic block

# Read in the synteny information
open(CSV, "<", $synfile) or die $!;

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

# Read in the syntenic block info

open(BLOCK, "<", $blockfile);
my @blocks;
my @block_coords;
<BLOCK>; # skip the header line
while(<BLOCK>) {
	last if $_ =~ /^\s/;
	my %block = {};
	my @fields = split;
	my @block = split(',', $fields[0]);
	my $chr = $fields[2];
	$block{'chr'} = $chr;
	foreach my $locus (@block) {
		$block{$locus} = $comp_info{$locus}[2];
	}
	push @blocks, \%block;

	# Record the chromosomal position of the first and last locus in the block (with the smaller value first)
	my @coords = sort {$a <=> $b} ($comp_info{$block[0]}[2], $comp_info{$block[-1]}[2]);
	push @block_coords, \@coords;
}
close BLOCK;


# Check to see if each locus is within a block
open(OUT, ">", $genome . '.synteny.mapped.out') or die $!;
print OUT join("\t", 'LOCUS', 'MAP_LG', 'LEFT', 'RIGHT', 'MAP_LEFT_POS', 'MAP_RIGHT_POS'), "\n";

foreach my $locus (keys %single_hits) {
	my $pos = $single_hits{$locus}[2];
	for (my $i = 0; $i < scalar(@blocks); $i++) {
		if ($single_hits{$locus}[0] == $blocks[$i]{'chr'}) { # next unless the locus is on the same chromosome as the block
		} else {
			next;
		}
		if ($pos >= $block_coords[$i][0] && $pos < $block_coords[$i][1]) { # The locus is within the block
			## Find the smallest interval that the locus can map to
			# make arrays with the sorted coordinates and loci
			my @loci;
			my @coords;
			my %temp_block = %{$blocks[$i]};
			delete $temp_block{'chr'};
			foreach my $marker (sort {$temp_block{$a} <=> $temp_block{$b}} keys %temp_block) {
				push @loci, $marker;
				push @coords, $temp_block{$marker};
			}

			my $left;
			my $right;
			for (my $j = 0; $j < scalar(@loci); $j++) {
				$left = $loci[$j];
				if ($pos > $coords[$j + 1]) {
					next;
				} else {
					$right = $loci[$j + 1];
					last;
				}
			}
			print OUT join("\t", $locus, $map_info{$left}[0], $left, $right, $map_info{$left}[2], $map_info{$right}[2]), "\n";

		} else { # the locus is not within the block
			next;
		}
	}
}

close OUT;


# Subroutines

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
