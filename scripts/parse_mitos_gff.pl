#!/usr/bin/perl

############################################################
#
# script to reformat the mitogenome annotation previously obtained with mitos. 
#
# Author: J. Gómez Garrido
############################################################

use warnings;
use strict;
use Getopt::Long;

my ($gff, $pref, $add, $discard);

GetOptions(
           'g:s'           => \$gff, 
           'p:s'           => \$pref,
        #    'a'              => \$add,   #if this option is given the lines with the id are going to be kept
        #    'b'               => \$discard, #If this option is given, the lines with the id are goin to be discarded
           );

open GFF, "<", "$gff";
my $c = 1;
my $id;
my $e;
while (<GFF>){
  chomp;
  my @line = split /\t/, $_; 
  if ($line[2] eq "gene" || $line[2] eq "ncRNA_gene"){
    $id = sprintf ($pref . "_g%02d",$c);
    my @gene = @line;
    $gene[2] = "gene";
    $gene[8]=~ s/ID=([^;]+);gene_id=/ID=$id;Name=/;
    my $gene = join "\t", @gene;
    print "$gene\n";
    $c++;
    $e=0;
  }
  elsif ($line[2] eq "exon"){
    $e++;
    $line[8] =~ s/Parent=([^;]+)/ID=$id.exon$e;Parent=$id.1/; 
    $line[7] = ".";
    my $exon = join "\t", @line;
    print "$exon\n";
  }
  else {
    $line[8] =~ s/ID=([^;]+);Parent=([^;]+);gene_id/ID=$id.1;Parent=$id.1;Name/;
    my $trans = join "\t", @line; 
    print "$trans\n";  
  }
}
close GFF;