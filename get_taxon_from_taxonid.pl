#!/bin/env perl
#AUTHOR:GIOVANNI MARQUES DE CASTRO
#TESTED on illumina paired end reads (hiseq)

use Bio::Taxon;
use Bio::Tree::Tree;
use Getopt::Long;
use warnings;
use strict;
my $binningdir="/home/giovannimc/projects/Carol_metagenomics-cacao/work/noadapter/version2/KAIJU/";
my $datadir="/home/giovannimc/projects/Carol_metagenomics-cacao/data/cutadapter/Filtered/";
my $classified_reads="";
my ($read1,$read2)="";
my $usage="USAGE:\nperl $0 -c CLASSIFICATION_FILE(KAIJU_LIKE) -r1 READ_1 -r2 READ_2 [optional: -rank phylum]\n";
my $rank='phylum';
#my $tree_functions = Bio::Tree::Tree->new();
GetOptions (#"length=i" => \$length,    # numeric
              "classification_file=s" => \$classified_reads,      # string
              "r1=s" => \$read1,      # string
              "r2=s" => \$read2,      # string
	      "rank=s"=>,\$rank		#string (defaults to phylum)
#              "verbose"  => \$verbose   # flag
)||die "$usage";

open(CLASSIFIED_READS,"<$classified_reads")||die "can not open '$classified_reads'\n$usage";
open(READ1,"<$read1")||die "cannot open read1 '$read1'\n$usage";
open(READ2,"<$read2")||die "cannot open read2 '$read2'\n$usage"; #BUG: pass if named pipe is empty

#f in /home/giovannimc/projects/Carol_metagenomics-cacao/work/noadapter/KAIJU/[^ABCD].combined.tsv

# Get one from a database
my $dbh = Bio::DB::Taxonomy->new(-source   => 'flatfile',
                                 -directory=> '/home/giovannimc/genomas_and_databases/NCBItaxonomy/BioTaxon',
                                 -nodesfile=> '/home/giovannimc/genomas_and_databases/NCBItaxonomy/nodes.dmp',
                                 -namesfile=> '/home/giovannimc/genomas_and_databases/NCBItaxonomy/names.dmp');
my $taxon = $dbh->get_taxon(-taxonid => '9606');
my %hashing;#accelerate the search
my %taxa_read;
my %filehandles; # to write multiple files 
while(<CLASSIFIED_READS>){
	if( $_=~/^C\t(.*)\t(.*)$/){
		my $read=$1;
		my $taxid=$2;
		if($hashing{$taxid}){ #has seen the taxid, avoid bio:tree fucntions 
			$taxa_read{$read}=$hashing{$taxid};
		}else{
			$taxon = $dbh->get_taxon(-taxonid => $taxid);
			my $tree = Bio::Tree::Tree->new(-node => $taxon);
			my $lineage = $tree->find_node(-rank => $rank);
			if($lineage){ #Lineage has phylum, or any other desired rank
				my $TAXON=$lineage->scientific_name;
				$hashing{$taxid}=$TAXON; #hashing to use Bio::Taxon and BIO::Tree::Tree once per desired taxon
				$taxa_read{$read}=$TAXON;
				if(!$filehandles{"$TAXON.R1"}){
					open($filehandles{"$TAXON.R1"},">$TAXON.R1.fastq")||die "CANNOT open $TAXON.R1.fastq to write\n";
					open($filehandles{"$TAXON.R2"},">$TAXON.R2.fastq")||die "CANNOT open $TAXON.R2.fastq to write\n";
				}
			}
		}
	}
}
my$c=0;
while(<READ1>){
	my$r1=$_;
	$r1=~s/^@|\/1$//g;
	chomp($r1);
	if($taxa_read{$r1}){	
		print {$filehandles{"$taxa_read{$r1}.R1"}} (   $_  .<READ1>.<READ1>.<READ1>);#dereference the filehandle to print
		print {$filehandles{"$taxa_read{$r1}.R2"}} (<READ2>.<READ2>.<READ2>.<READ2>);
	}else{ #read pair not classified
		   	<READ1>.<READ1>.<READ1>;#jump 3 lines
 		<READ2>.<READ2>.<READ2>.<READ2>;#jump 4 lines
	}

}	

