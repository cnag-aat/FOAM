#! /usr/bin/perl

# Program to get the specified contig/scaffold from an assembly FASTA
 
# Fernando Cruz, 01-02-2023

use strict;
use Getopt::Long;
use File::Basename;

# define variables
my $file;
my $target;

my @data=();


my $line;
my $scaf_id;

GetOptions(

           'fasta|f:s'   => \$file,
           'target|t:s'   => \$target  #
);


my ($basename,$path,$ext) = fileparse($file,qw(\.fa));

#MAIN FUNCTIONS
#open the list file

select_sequence();

#SUB-ROUTINES

sub select_sequence {
  # cleaning Input FASTA
   # REDIRECT OUTPUT TO STDOUT

open (IN, "$file") or die print "I cannot open file $file \n";
{

  while(<IN>)
  {
      chomp;
      $line=$_;
      
      if ($line=~/^>/) 
      {   $scaf_id=$line;
          $scaf_id =~ s/\>//g;
          if ($scaf_id eq $target) {
                       print "\>$scaf_id\n";
         }
      }# FASTA header

      else{
            if ($scaf_id eq $target){ 
               print "$line\n";    #print sequence
            }
      
      }#FASTA sequence

  }
    

  close (IN) or die print "cannot close the file $file\n";

}#FASTA closed


}#clean_fasta sub-routine
