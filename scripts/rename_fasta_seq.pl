#! /usr/bin/perl

# Program to rename the sequence of a single FASTA (non-multifasta, just with just a sequence entry  
 
# Fernando Cruz, 25-05-2023

use strict;
use Getopt::Long;
use File::Basename;

# define variables
my $file;
my $seqname;

my @data=();


my $line;
my $scaf_id;

GetOptions(

           'fasta|f:s'   => \$file,
           'seqname|n:s'   => \$seqname  #
);


my ($basename,$path,$ext) = fileparse($file,qw(\.fa));

#MAIN FUNCTIONS
#open the list file

rename_sequence();

#SUB-ROUTINES

sub rename_sequence {
 
open (IN, "$file") or die print "I cannot open file $file \n";
{

  while(<IN>)
  {
      chomp;
      $line=$_;
      
      if ($line=~/^>/) 
      {   
          print "\>$seqname\n";
         
      }# FASTA header

      else{
          
          print "$line\n";   #print sequence line by line      
      
      }#FASTA sequence

  }
    

  close (IN) or die print "cannot close the file $file\n";

}#FASTA closed


}#clean_fasta sub-routine
