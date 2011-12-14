#!/usr/bin/perl -w

my $usage = <<END;

	$0 - call omssa and gzip the result
Usage:
	
	$0 <working dir> <ommsa command> <options> <params file> <result file> <user mods file>
	
Requires:

	ommsa

END

if (scalar @ARGV < 6) {
	print $usage;
	exit 1;
}

use strict;



# inputs are <working dir> <ommsa command> <options> <params file> <result file> <user mods file>

my $workingDir = shift @ARGV;
my $omssaCl = shift @ARGV;
my $options = shift @ARGV;
my $paramsFile = shift @ARGV;
my $expected_resultFile = shift @ARGV;
my $user_modsFile = shift @ARGV;

print "input:\n"
."\tcwd        = $workingDir\n"
."\tomssaCl    = $omssaCl\n"
."\toptions    = $options\n"
."\tparamsFile = $paramsFile\n"
."\tresult     = $expected_resultFile\n"
."\tuser_mods  = $user_modsFile\n\n";

# call omssa
my $com = "cd $workingDir && $omssaCl $options $paramsFile -mux $user_modsFile";

my $result = "";

print "command:\n";
print "$com\n\n";

$result = system $com;

# see if the result is non zero
if($result!=0){
	print("ERROR: omssa call failed with code=$result\n");
	exit $result;
}
# gzip the result
$com = "gzip  $expected_resultFile";
print "gzipping the result. Command:\n$com\n\n";
$result = system $com;
if($result != "0"){ 
	print("ERROR: gzipping failed with code=$result\n");
	exit $result;
}





