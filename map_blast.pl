#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;

pod2usage(-verbose => 1) if @ARGV == 0;

my $fastafile;
my $argfile;
my $outfile = 'results.out';
my $blast_db;
my $db_dir = '/home/blast/db/';
my $min_length = 50;


GetOptions(		'fastafile|f=s' => \$fastafile,
				'argfile|a=s' => \$argfile,
				'outfile|o=s' => \$outfile,
				'blast_db|b=s' => \$blast_db,
				'db_dir|z=s' => \$db_dir,
				'min_length|m=s' => \$min_length,
	);

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

my $seq_in = Bio::SeqIO->new(-file => $fastafile,
								-format => 'fasta');

open(LOC, ">", 'loci.txt') or die $!;								
my @seqs;								
while (my $seq = $seq_in->next_seq) {
	push @seqs, $seq;
	print LOC $seq->display_id, "\n";
}
close LOC;

my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
    -db_name => $blast_db, -db_dir => $db_dir
);

open(OUT, ">", $outfile) or die $!;
open(TAB, ">", $outfile . '.tab') or die $!;

foreach my $seq (@seqs) {
	my $result = $fac->blastn( -query => $seq,
                        -method_args => \@args );
	
		while (my $hit = $result->next_hit ) {
			while (my $hsp = $hit->next_hsp ) {
				if ($hsp->length('total') >= $min_length ) {
					print OUT join("\t", $result->query_name, $hit->name, $hit->description, $hit->significance, $hsp->start('hit'), $hsp->end('hit')), "\n";
					print TAB join("\t", $result->query_name, $hit->name, $hit->description, $hit->significance, $hsp->start('hit'), $hsp->end('hit')), "\n";
					print OUT $hsp->seq_str('query'), "\n";
					print OUT $hsp->seq_str('match'), "\n";
					print OUT $hsp->seq_str('sbjct'), "\n";
					print OUT "=========================================================================\n\n";
					#print join("\t", $result->query_name, $hit->name, $hit->description, $hit->significance), "\n";
				}
			}  
		}
}
$fac->cleanup;
