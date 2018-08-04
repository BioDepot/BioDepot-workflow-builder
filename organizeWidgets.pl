#!/usr/bin/perl

my ($topDir)=@ARGV;
my $dest="biodepot/RNA_seq";
my @dirs = split(' ',`ls $topDir`);
foreach my $dir (@dirs){
	my $fullDir = "$topDir/$dir";
	print STDERR "working on $fullDir\n";
	my $py= "$fullDir/$dir.py";
	system ("ln -s ../../$py $dest/OW$dir.py")
}
