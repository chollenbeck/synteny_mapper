#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Copy;
use Cwd;

pod2usage(-verbose => 1) if @ARGV == 0;

# Get the working directory

my $cwd = getcwd();

# Read in the config file
my $config = shift;

my %config;
open(CON, "<", $config) or die $!;
while(<CON>) {
	next if $_ =~ /^\s/;
	next if $_ =~ /^#/;
	chomp;
	my ($key, $value) = split;
	$config{$key} = $value;
}


# Determine which comparison species are needed
my $comps = add_comps(\%config);
my %comps = %{$comps};


# Run the local BLAST
my @reports;
print "Performing local BLAST search...\n";
foreach my $species (keys %comps) {
	system("perl map_blast.pl -f $config{'fasta_file'} -a $config{'blast_args'} -b $comps{$species} -o $species.report.out -z $config{'db_dir'}");
	push @reports, "$species.report.out.tab";
}

# Parse the blast reports and combine locus position information into a single file ('synteny.out')
print "Parsing BLAST reports...\n";
system("perl parse_blast_report.pl -m $config{'map_file'} -l loci.txt -o synteny.out -f " . join(' ', @reports));

# Identify syntenic blocks for each species
print "Identifying syntenic blocks for each species...\n";
foreach my $species (keys %comps) {
	print "\tProcessing $species...\n";
	system("perl id_syntenic_blocks.pl -m $config{'map_file'} -i synteny.out -g $species");
}

# Combine syntenic block data
open(ALL, ">", 'all.blocks.tab') or die $!;
print ALL join("\t", 'COMP_SPECIES', 'START_LOCUS', 'END_LOCUS', 'NUM_LOCI', 'COMP_CHR', 'COMP_START', 'COMP_END', 'COMP_SIZE', 'MAP_LG', 'MAP_START', 'MAP_STOP', 'MAP_SIZE'), "\n";
foreach my $species (keys %comps) {
	open(TAB, "<", "${species}.synteny.tab") or die $!;
	<TAB>;
	while(<TAB>) {
		print ALL "$species\t", $_;
	}
	close TAB;
}
close ALL;

# Synteny map loci, if specified
if ($config{'synteny_map_loci'}) {
	print "Synteny mapping loci...\n";

	my @synteny_outfiles;
	foreach my $species (keys %comps) {
		print "\tMapping loci for $species...\n";
		system("perl map_loci.pl -i $config{'synteny_map_loci'} -s synteny.out -l ${species}.synteny.sum -g $species -b $comps{$species} -z $config{'db_dir'} -m $config{'min_hit_length'} -a $config{'blast_args'}");
		push @synteny_outfiles, $species . '.synteny.mapped.out';
	}

	system("perl consol_synteny_map_loci.pl -m $config{'map_file'} -o all.synteny.mapped.out -f " . join(' ', @synteny_outfiles));
}

# Check to see if there is a results directory

my $main_dir = getcwd();

opendir(CWD, $main_dir);
my @res = grep { /^results_/ && -d "$main_dir/$_" } readdir(CWD);

my $highest = 0;
foreach my $dir (@res) {
	if ($dir =~ /results_(\d+)/) {
		if ($1 > $highest) {
			$highest = $1;
		}
	}
}
my $num = $highest + 1;
mkdir 'results_' . $num;

my $res_dir = $main_dir . '/results_' . $num;

chdir $res_dir;

foreach my $species (keys %comps) {
	mkdir $species;
	my @files = glob "$main_dir/${species}*";
	foreach my $file (@files) {
		File::Copy::move($file,"$res_dir/$species/") or die "Could not move $file: $!\n";
	}
}

chdir $main_dir;

# Organize the results files
File::Copy::copy('config.txt',"$res_dir/") or die "Could not move file: $!\n";
File::Copy::move('loci.txt',"$res_dir/") or die "Could not move file: $!\n";
File::Copy::move('all.blocks.tab',"$res_dir/") or die "Could not move file: $!\n";
File::Copy::move('all.synteny.mapped.out',"$res_dir/") or die "Could not move file: $!\n";
File::Copy::move('synteny.out',"$res_dir/") or die "Could not move file: $!\n";


print "Program Finished\n";



sub add_comps {

	my %config = %{$_[0]};
	my %comps;

	# Danio rerio
	if ($config{'drer'} eq 'TRUE') {
		$comps{'drer'} = 'Drer_v9';
	}

	# Gasterosteus aculatus
	if ($config{'gacu'} eq 'TRUE') {
		$comps{'gacu'} = 'Gacu_v1';
	}

	# Oreochromis niloticus
	if ($config{'onil'} eq 'TRUE') {
		$comps{'onil'} = 'Onil_v1_1';
	}

	# Tetraodon nigroviridis
	if ($config{'tnig'} eq 'TRUE') {
		$comps{'tnig'} = 'Tnig_v8';
	}

	# Takifugu rubripes
	if ($config{'trub'} eq 'TRUE') {
		$comps{'trub'} = 'Trub_v5';
	}

	# Lates calcarifer
	if ($config{'lcal'} eq 'TRUE') {
		$comps{'lcal'} = 'Lcal_v1';
	}

	# Dicentrarchus labrax
	if ($config{'dlab'} eq 'TRUE') {
		$comps{'dlab'} = 'Dlab_v1';
	}

	return \%comps;

}
