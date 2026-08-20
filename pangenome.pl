#!/usr/bin/perl

open (IN, $ARGV[0]); #open Orthogroups.csv
open (IN2, $ARGV[1]); #open Orthogroups_UnassignedGenes.csv

$first_line = <IN>;
@header = split (/\t/, $first_line);
$numb_of_genomes = scalar (@header) - 1;
$first_column = 1;
$last_column = $numb_of_genomes;
$first_line = <IN2>;

@ortho = <IN>;
@single = <IN2>;
$c=0;
open (OUT, ">core.txt");
open (OUT2, ">shared.txt");
open (OUT3, ">singletons.txt");

for $i (@ortho){
	@line = split (/\t/, $i);
	@result="";
	for ($i2 = $first_column; $i2 <= $numb_of_genomes; $i2++){
		$line[$i2] =~ s/\s//g;
		push @result, $line[$i2];
		if ($line[$i2] =~ m/./){
			$c++;
		}
	}
	if ($c == $numb_of_genomes){
		shift @result;
		$core = join "\t", @result;
		print OUT $core."\n";
	}else{
		shift @result;
		$shared = join "\t", @result;
		print OUT2 $shared."\n";
	}
	$c=0;
	
}
close IN;
close OUT;
close OUT2;

for $i (@single){
	@line = split (/\t/, $i);
	@result = "";
	for ($i2 = $first_column; $i2 <= $numb_of_genomes; $i2++){
		$line[$i2] =~ s/\s//g;
		push @result, $line[$i2];
		
	}
	shift @result;
	$singleton = join "\t", @result;
	print OUT3 $singleton."\n";
}

