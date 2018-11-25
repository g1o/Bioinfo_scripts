#!/bin/env perl
use Bio::Taxon;
use Bio::Tree::Tree;
use Getopt::Long;
use warnings;
use strict;
my $classified_reads="PATH/version2/KAIJU/A.combined.tsv.head";
#my $tree_functions = Bio::Tree::Tree->new();
GetOptions (#"length=i" => \$length,    # numeric
              "input=s" => \$classified_reads,      # string
#              "verbose"  => \$verbose   # flag
);

my $binningdir="PATH/noadapter/version2/KAIJU/";
my $datadir="PATH/data/cutadapter/Filtered/";
open(CLASSIFIED_READS,"<$classified_reads")||die "cannot open $classified_reads\n";

# Get one from a database
my $dbh = Bio::DB::Taxonomy->new(-source   => 'flatfile',
                                 -directory=> '/home/giovannimc/genomas_and_databases/NCBItaxonomy/BioTaxon',
                                 -nodesfile=> '/home/giovannimc/genomas_and_databases/NCBItaxonomy/nodes.dmp',
                                 -namesfile=> '/home/giovannimc/genomas_and_databases/NCBItaxonomy/names.dmp');

my $taxon = $dbh->get_taxon(-taxonid => '9606'); 
my %hashing; #keep taxid => 'other rank' for faster search.
while(<CLASSIFIED_READS>){
        if( $_=~/^C\t(.*)\t(.*)$/){
my $read=$1;
my $taxid=$2;
        if($hashing{$taxid}){
                print "$read\t tx:$taxid\tphylum is ",$hashing{$taxid},"\n";
        }else{

        $taxon = $dbh->get_taxon(-taxonid => $taxid);
        my $tree = Bio::Tree::Tree->new(-node => $taxon);
        my $lineage = $tree->find_node(-rank => 'phylum');
        #print "tx:$taxid\t".$taxon->rank.":rank\t". $tree->number_nodes . " nodes\n";
        if($lineage){ #Lineage has phylum, or any other desired rank
        print "$read\t tx:$taxid\tphylum is ",$lineage->scientific_name,"\n";
    }
    }
}
}

