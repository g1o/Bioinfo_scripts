#!/bin/env perl 
#########################################
#updated: 2018/04/04
#version: 0.9.2
#roadmap: check and join with last downloaded database to complete the update
#roadmap: CHECK obsolte versions and remove 
#roadmap: auto check date timestamp
#1-problems: WP sequences are not in NT (has no link to CDS from efetch) even downloading entire nt would not find such sequences
#1.1-check which has not being downloaded and download the contig/scaffold position it is in
#License: GLP3
#Author: Giovanni Marques de Castro 
#########################################
use warnings;
use strict;
use Time::HiRes qw(time);

#if ( -e "TINY.XML") {die "ABORTING: TINY.XML exists!!!\n"};

my $lastupdate="2018/04/06";
my $query = '%22+'; #%22 = " ; + = space ;

#assemble the esearch URL
my ($dia,$mes,$ano)= (localtime)[3,4,5];
my $hoje	= sprintf ("%.4d/%.2d/%.2d", $ano+1900, $mes+1, $dia+2); #mais dias para ter certeza do fuso horario
my $date	= "mindate=$lastupdate&maxdate=$hoje&datetype=pdat";
my $db		= "protein";
my $base	= 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url		= $base . "esearch.fcgi?db=$db&term=$query&$date&usehistory=y";

#post the esearch URL
my $output = get($url);
#print $url ."\n";
#print $output;
#parse WebEnv, QueryKey and Count (# records retrieved)

my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

#Open output file for writing
open(TINYXML,">TINY.XML")||die "Can't open tiny.xml to write\n";
open(CDSFASTA,">CDS_COI.fasta")||die "Can't open CDS_COI.fasta to write\n";
#open(COINUC,">COI_NUCCORE.fasta")||die "Can't open COI_NUCCORE.fasta to write\n";
#open(DEBUG,">DEBUG.XML")||die "Can't open to DEBUG\n";
#Retrieve data in batches of at most 500 (10000 is max from NCBI, but no so efficient and I am not sure why)
my $retmax = 150;
my $children = 0;
print "Process ID: $$ \n";
for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
	sleep 1; #BE CAREFULL, MORE THAN 3 SIMULTANEOUS ACCESS MAY BAN YOU FROM NCBI
	if ($children == 2) { #Control amount of forks (simultaneous downloads); 
		my $finished = wait();
		$children--;
		#		print "PID $$ finished $finished\n";
	}
	my $pid = fork;
	if (not defined ($pid) ) {
		warn "could not fork\n";
		next;
	}
	#Parent process
	if ($pid) { 
		$children++;
		warn "In the parent process PID ($$), Child pid: $pid Num of fork child processes: $children \t ret=$retstart\n";
		#			$retstart += $retmax
	} 
	#Child process
	else{
		print "In the child process PID ($$) \t Downloading $retstart to ".($retstart + $retmax)." of $count\n" ;
		my $TSeq=0 ;
		while ($TSeq != $retmax){
			sleep 1;
			my $dtime=time;
			my $efetch_url = $base ."efetch.fcgi?db=$db&WebEnv=$web";
			$efetch_url .= "&query_key=$key&retstart=$retstart";
			my $cds_url=$efetch_url;
			$efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=xml";
			$cds_url .= "&retmax=$retmax&rettype=fasta_cds_na";
			my @efetch_out = get($efetch_url);
			my @cds_out = get($cds_url);
			$dtime=time-$dtime;
			print "download has taken $dtime seconds \n";
			$dtime=0;
			$TSeq = grep (/<TSeq>/,@efetch_out );
			if ($TSeq != $retmax) {
				warn "WARNING:  downloaded $TSeq of $retmax\t...retrying if this isn't the last batch\n";
				sleep 1;
			}else{
				print TINYXML @efetch_out;
				print CDSFASTA @cds_out;
			}
			if( $retmax>($count-$retstart) ){
				print TINYXML @efetch_out;
				print CDSFASTA @cds_out;
				last;
			}
		}
		exit;
	}

}

for (1 .. $children) {
   my $pid = wait();
   print "Parent saw $pid exiting\n";
}

print "Parent ($$) ending\n";

close (TINYXML);
close (CDSFASTA);

warn "Downloaded $count of $count\n";


##################################################################################


sub get {
	my $url=shift @_;
	return `wget -q -O - "$url"`;
	#return `wget --header='Accept-Encoding: gzip' -q -O - "$url"`;
}

