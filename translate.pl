#!/bin/env perl
#
#Translate a fasta file in the specified frame or in all six frames if no frame was specified
#The first argument is the Fasta file, the second is the minimum length of the tranlation (to exclude lots of  small proteins), and the third is the frame [0,1,2] 
#If you wish to use another codon table, edit this code (line 53: -codontable_id )
# 
# Require BioPerl
#
#Usage: perl ThisScript.pl file.fasta 0 1
#
#	AUTHOR: Giovanni Marques de Castro
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>
#
#
use warnings;
use strict;
use Bio::SeqIO;
my $file=shift;
my $minlength=shift;
my $frame=shift;
my $sequences = Bio::SeqIO->new( 
    -file   => $file,
    -format => "fasta",
);
if(! $minlength){
	$minlength=2;
}

while ( my $dna = $sequences->next_seq ){
	if(defined($frame)) {
		my $protein = $dna->translate( 
				-codontable_id => 1, # standard genetic code
				-frame         => $frame, #reading-frame offset 0
				);
		print ">",$dna->display_id, "\n";
		print $protein->seq, "\n\n";
	}
	else{
		for( $frame=0;$frame<3;$frame++ ){
			my $protein = $dna->translate(
					-codontable_id => 1, # standard genetic code
					-frame         => $frame, #reading-frame offset 0
					);
			my $id = ">".$dna->display_id . " frame_$frame";
			my @orfs = split(/\*/,$protein->seq);
#			my @orfs=split(/\*/,$orfs);
			my $c=0;
			foreach my $l (@orfs){
				if (length($l)<=$minlength){next};
				print $id."_ORF.$c\n";
				print $l."\n";
				$c++;
			
			}
			$protein = $dna->revcom->translate(
					-codontable_id => 1, # standard genetic code
					-frame         => $frame, #reading-frame offset 0
					);
			$id = ">".$dna->display_id . " revframe_$frame";
			@orfs = split(/\*/,$protein->seq);
			$c=0;
			foreach my $l (@orfs){
				if (length($l)<=$minlength){next};
				print $id."_ORF.$c\n";
				print $l."\n";
				$c++;
			}
		}
		$frame=undef;
	}

}
