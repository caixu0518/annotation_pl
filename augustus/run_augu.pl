#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##-
my $command = "commands.list";
open IN0, $in0;
while(<IN0>){
  chomp;
  my $line = readpipe("head -1 $_");
  print  $line, "\n";
}
close IN0;
