#!/usr/bin/perl

my $repo="/home/lhhung/ssh_repos/workflows";
my @dirs=split(/\n/,`find $repo -type f -name *.ows -exec dirname {} \\;`);
foreach my $dir (@dirs){
    print "$dir\n";
}