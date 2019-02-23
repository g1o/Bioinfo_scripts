#!/usr/bin/perl
=aheader
*************************************************************************
 * fastq_recover - A program to recover reads from corrupted fastq files
 *
 * 2019 Giovanni Marques de Castro (gmdc@outlook.com)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation
 * Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307 USA
 ************************************************************************
=cut

use warnings;
use strict;
use Getopt::Long;
# TAKES an fastq as input and output all lines that agree to the fastq format. the first 4 lines must be OK. 
my $c=0;
my $ok=1;
my $verbose;
my $corrupted_fastq;
my $additional_fastq;
my %read_ids;
my $acomment='';
my $usage= "Usage:
perl $0 -corrupted FILE.fastq [optional: -a additional.fastq] [optional: -verbose] [optional: -comment 'Something to add to the comment line of the additional reads that will be used to fill in the gaps']
#\tBy default, it outputs reads in STDOUT\n";

GetOptions ("corrupted_fastq=s" => \$corrupted_fastq,         # File to be checked and have its reads extracted
              "additional_fastq=s"   => \$additional_fastq,   # Optional Fastq which might have reads that were missing from the corrupted file. This file might be one after reads were processed or a copy of the original.
              "comment=s" => \$acomment,
              "verbose"  => \$verbose);   # flag

if(!$corrupted_fastq){die "$usage"};

warn "OPENING FILES...\n";
open(CORRUPTED_FASTQ,"<$corrupted_fastq")|| die "Could not open the required corrupted fastq file: '$corrupted_fastq' \n$usage";
if($additional_fastq){
	open(ADDITIONAL_FASTQ,"<$additional_fastq")|| die "Could not open the optional fastq file: '$additional_fastq' \n$usage";
}
##INITIALIZE - the first 4 lines must be a good read. 
my $id=<CORRUPTED_FASTQ>;
my $sequence=<CORRUPTED_FASTQ>;
my $comment=<CORRUPTED_FASTQ>;
my $qual=<CORRUPTED_FASTQ>;
#for illumina (eg. @E00576:98:HFC2FCCXY:1:1101:6948:1608)
my $MachineID=$id;
$MachineID=~s/^(@[-A-Z0-9]+:\d+:\w+).*([\s\/]+.*\n)/$1/ || die "Error at the first line! First read id not recognized! (Or is it the regex?)\n";
my $idEND=$2;
print $id.$sequence.$comment.$qual;
$id=~s/[\/\s]+.*\n//;
$read_ids{$id}=0; #changing 0 to undef may save RAM, as it is around 3 GB of RAM is used for a 7 GB fastq.
warn "LOOPING ...\n";
##Now begin to output only nice and good reads; 
while(<CORRUPTED_FASTQ>){
	if($c==0 || $ok == 0 ){
		if($_=~/^$MachineID:\d+:\d+:\d+:\d+[\s\/]/){
			$id=$_;
			$ok= 1;
			$c = 0;
		}
		else{
			$ok= 0;
		}
		$c++;
		next;
	}
	if($c==1 && $ok == 1){
		$sequence=$_;	
		chomp($sequence);
		if($sequence=~/[^ATCGNatcgn]/){ # if has any non nucl char, it is wrong
			$ok=0;
			if($verbose){
				warn "wrong sequence at $.\n";
			}
		}
		else{
			$sequence.="\n"; 
		}
		$c++;
		next;
	} 
	if($c==2 && $ok == 1 ){
		if(/^\+.{0,100}/){ # In very rare cases something is written here, so, this check may be changed 
			$comment=$_;
		}
		else{ 
			$ok=0;
			if($verbose){
                                warn "Comment line have over 100 chars, probably wrong, skipping read at line $.\n";
                        }
		}
		$c++;
		next;
	}
	if($c==3 && $ok == 1 ){
		if ( length($_) == length($sequence) ){
			$qual=$_;
			print $id.$sequence.$comment.$qual;
			if($additional_fastq){
				$id=~s/[\/\s].*\n//;
				$read_ids{$id}=0; #changing 0 to undef may save RAM
			}
			$c=-1;
		}
		else{ 
			$ok=0;
			if($verbose){
				warn "Length of quality does not mach sequence length at $.\n";
			}
		}
	}
	$c++;
	if( $. % 2000000 == 0 ){
		warn "$. lines processed (about ". $./4 ." reads) \n";
	} 
}
close(CORRUPTED_FASTQ);
warn "Finished looping $corrupted_fastq\n";
if($additional_fastq){ #Obviously, the additional fastq should be a normal one
	warn "Looping in $additional_fastq\n";
	$c=0;
	while(<ADDITIONAL_FASTQ>){
		my $line=$_;
		if($c!=0){$c--;next}; # skip if it is not a header;
		$line=~s/\/\d.*\n|\s+.*\n//;
		if( !exists($read_ids{$line})){
			print $line.$idEND;	 #@id
			$line=<ADDITIONAL_FASTQ>;#sequ
			$line.=<ADDITIONAL_FASTQ>;
			chomp($line);            #+
			$line.=$acomment."\n";
			$line.=<ADDITIONAL_FASTQ>;#qual
		print $line;	 
		}else{
			$c=3;
		}
		if( $. % 2000000 == 0 ){
			warn "$. lines processed (about ". $./4 ." reads) \n";
		}
	}
}
close(ADDITIONAL_FASTQ);
