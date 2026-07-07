#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($list, $add, $discard);

GetOptions(
           'l:s'           => \$list, ##list of reads to delete
           'a'              => \$add,   #if this option is given the the reads in the list are going to be kept
           'b'               => \$discard, #If this option is given, the reads in the list are going to be discarded
           );

open LIST, "<", "$list";
my %ids;
while (<LIST>){
  chomp;
  if ($_ =~ m/@/){
    $_ =~ s/@//;
  }
  $ids{$_}++;
}
close LIST;

my $reads;
my $l = 1;
my $keep = 0; 
while (<STDIN>){
  chomp;
  if ($l == 1) {
    my @line = split /\s+/, $_;
    if ($line[0] =~ m/^@/){
      $line[0] =~ s/@//; 
    }
    if ($discard && !exists $ids{$line[0]}){
      print "$_\n";
      $keep = 1;
    }
    elsif ($add && exists $ids{$line[0]}) {
      print "$_\n";
      $keep = 1;
    }
    else {
      $keep = 0;
    }
  }
  else {
    print "$_\n" if ($keep == 1); 
    $l = 0 if ($l == 4);
  }
  $l++;
}
