#!/bin/perl
my @dirs= split(' ',`find */widgets/*/* -maxdepth 0 -type d`);
foreach my $dir (@dirs){
	print STDERR "working on $dir\n";
	my @parts=split('/',$dir);
	my $name=$parts[-1];
	if ($name ne "icon") {
	 print STDERR "searching for widget $name\n";
  $copyDir=`find ../widgets/*/$name/Dockerfiles -maxdepth 0 -type d`;
  chomp($copyDir);
  #print "copyDir is $copyDir\n";
  if ($copyDir){
   print "cp -r $copyDir\ $dir\/.\n";
   system("cp -r $copyDir\ $dir\/.");
 	}
	}
}
