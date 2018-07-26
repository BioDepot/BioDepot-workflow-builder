#!/usr/bin/perl

#organizes the widgets into directories
my ($dir) =@ARGV;

print STDERR "working on $dir\n";
my (@files)= split(' ',`ls $dir/OW*.py`);
foreach my $file (@files) {
	print STDERR "working on python $file\n";
	my ($d1,$iconFile)=split('"',`fgrep icon $file`);
	my (@parts)=split(/\//,$iconFile);
	my $bareIconFile=$parts[-1];
	printf "icon file is %s\n",$iconFile;
	my ($d1,$jsonFile)=split('"',`fgrep \json $file`);
	my @parts=split(/\//,$jsonFile);
	my $bareJsonFile=$parts[-1];
	printf "json file is %s\n",$jsonFile;
	my $wDir=substr($file,0,-3);
	$wDir=~s/\/OW/\//;
	my @parts=split(/\//,$wDir);
	my $barePyFile="$parts[-1].py";
	my $bareDir=$parts[-1];
	system ("mkdir -p $wDir");
	system ("mv $dir/icons/$bareIconFile $wDir/$bareIconFile");
	system ("mv $dir/json/$bareJsonFile $wDir/$bareJsonFile");
 system ("mv $file $wDir/$barePyFile");
 system ("ln -s ../$bareDir/$bareIconFile  $dir/icons/$bareIconFile" );
 system ("ln -s ../$bareDir/$bareJsonFile $dir/json/$bareJsonFile");
 system ("ln -s  $bareDir/$barePyFile $file");	
}
