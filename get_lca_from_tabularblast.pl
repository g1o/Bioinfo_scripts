#!/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Taxon;
use Bio::Tree::Tree;
my $tree_functions = Bio::Tree::Tree->new();
my $HOME=$ENV{"HOME"};
my ($help,$blast,$simple);
my $maxbitscore=0;
my $topp=1;
my $minbitscore=50;
my $query='';
my @taxids;
my $usage="$0 -in file.blast
Mandatory argument:
-in (FILE)            The resulting Blast+ file in tabular (oufmt 6) format with the 12 colunm for taxids
Optional arguments: 
-top_percent ($topp)      A number (0-100). Of the maximun bitscore, allow this much percent of variation to include in the LCA search
-minscore ($minbitscore)        A number. Minimum bitscore to include in the LCA search.\n";

GetOptions (
	#"length=i" => \$length,    # numeric
	"in=s"   => \$blast,      # string
	"simple"  => \$simple,    # flag
	"top_percent=i" => \$topp, # Of the maximun bitscore, allow this much of variation to include in the LCA search
	"minscore=i" => \$minbitscore, # minimum bitscore to include in the LCA search 
	"help" => \$help
);
if (!$blast || $help) {die "$usage"};
my $dbh = Bio::DB::Taxonomy->new(-source   => 'flatfile', #SET DIR TO NCBI TAXONDUMP
	-directory=> "$HOME/Documentos/projects/ncbi-taxon/ncbinodeindex/",
	-nodesfile=> "$HOME/Documentos/projects/ncbi-taxon/nodes.dmp",
	-namesfile=> "$HOME/Documentos/projects/ncbi-taxon/names.dmp"
);

open(OUT,'>',"OUTtxid.tsv")||die"COULD not open out to write\n";
open(BLASTOUT,"sort --parallel=4 -k1,1 -k12,12gr -k3,3gr  $blast  |grep '^[^#]'|" )||die "Could not open $blast \n $usage";
#open(BLASTOUT,"head -n 2000 $blast |sort --parallel=4 -k1,1 -k13,13gr -k3,3gr  |grep '^[^#]'|" )||die "Could not open $blast \n";
my $tab='NULL';
while ($tab){
	$tab = <BLASTOUT>;
	my @hit=(0)x14;
	if($tab){
		@hit=split(/\t/,$tab);
		if ($hit[12] eq "N/A\n" ){ next; } #ignore lines with no taxid
	}
	if ($query eq $hit[0]){
		#Add new taxid if pass
		if ($hit[11]+$maxbitscore*$topp/100 >= $maxbitscore && $hit[11] >= $minbitscore ){ #&& $hit[2] >= 90 ){ 
			chomp($hit[12]);
			push(@taxids,$hit[12]);
		}
	}
	else{
		if ($query ne '' ){
			my @taxons=();
			foreach my $id (@taxids){
				if ($id == 410659 || $id == 256318){ next;} #ignore unclassified and metagenomic sequences.  
				my $taxon=$dbh->get_taxon(-taxonid => $id);
				#			print "$id\n";
				if (!$taxon ){
					warn "$id not in nodes\n";
					next;
				}
				push (@taxons,$taxon);
				#				print $taxon."\t".$id."\n";
			}
			if (@taxons > 1){
				my $lca = $tree_functions->get_lca(@taxons);
				my @lineage = $tree_functions->get_lineage_nodes($lca);
				push(@lineage,$lca);
				$lca= "$query\t" . $tree_functions->get_lineage_string($lca);
				$lca =~ s/(.*);(.*?)$/$1;$2;Other $2;$2($query)/;
				foreach my $node (@lineage){
					if($node->rank ne "no rank" && $node->rank=~/phylum|^class|^family/){
						print OUT "$query\t".$node->rank.";".$node->scientific_name."\n";
					}
				}
				print "\n";
				if($lca=~/cellular organisms$/){$lca =~ s/(.*)\t(.*?)$/$1\t$2;other $2;$2($query)/;}
				print $lca."\n";
			}elsif(@taxons==1){
				my $lca="$query\t" . $tree_functions->get_lineage_string($taxons[0]);
                                my @lineage = $tree_functions->get_lineage_nodes($taxons[0]);
                                push(@lineage,$taxons[0]);
				$lca =~ s/(.*);(.*?)$/$1;$2;Other $2;$2($query)/;
				if($lca=~/cellular organisms$/){$lca =~ s/(.*)\t(.*?)$/$1\t$2;other $2;$2($query)/;}
                                foreach my $node (@lineage){
					if($node->rank ne "no rank" && $node->rank=~/phylum|^class|^family/){
                                                print OUT "$query\t".$node->rank.";".$node->scientific_name."\n";
                                        }
                                }
				print $lca."\n";
			}else{
				print "$query\t" . "Unindentified;$query" ."\n"; #Hits muito divergentes ou que não enão no nodes. (root, life) [taxid:1] 
			}
		}
		$query=$hit[0];
		$maxbitscore=$hit[11];
		@taxids=();
		if ($maxbitscore >= $minbitscore){
			chomp($hit[12]);
			push(@taxids,$hit[12]);
		}
	}
}
