#/bin/env perl
use strict;
use warnings;

open(my $FILEOFFILES,shift)||die "Can't open file of files!\n";
open(my $GITOTAXID,shift)||die "Can't open gi_taxid_nucl.dmp !\n";
my %gitotax;
my $gi="";

while (<$FILEOFFILES>){  
        my$file=$_;
        chomp$file;
        open(my$FH,$file)||die"Could not open file $file ! Is it a file of files?\n";
        $gi=<$FH>;
        if($gi!~/gi/){
		warn "$file does not contain a gi in the first identifier! skipping it...\n"; 
		next;	
	}
        $gi=~s/>gi\|(\d+)\|.*/$1/;
	chomp$gi;
        close($FH);
	if($gitotax{$gi}){
	$gitotax{$gi} .= "$file\tTAXID";
	}else{
	$gitotax{$gi}  = "$file\tTAXID";
	}
}

while (<$GITOTAXID>){
        $gi = $_;
	$gi=~s/\t.*//;
	chomp$gi;
	if($gitotax{$gi}){
		$gitotax{$gi} =~ s/TAXID/$_/g;
		print $gitotax{$gi};
	}
}
