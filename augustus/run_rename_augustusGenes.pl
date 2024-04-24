#!/usr/bin/perl -w
use strict;

##- format gff3
open IN0, "augustus.gff3";
open OUT0, ">augustus.form.gff3";
while(<IN0>){
   next, if(/^\s+$/);

  if(/^[^#]/){
     print OUT0 $_;
  }

}
close IN0;
close OUT0;
