#! /usr/bin/perl

# Fernando Cruz, 09-01-2020
# Program to select ONT reads matching certain contaminant in the preprocessed fastq reads (*.contam.fastq.gz)
# By default is just keeping those matching the standard Nanopore spike-in: the DNA Control Sequence 3.56Kb of the lambda phage genome

use strict;
use Getopt::Long;
use File::Basename;

# define variables
my $contaminant="DNA_CS";# default contaminant lambda phage
my $count=0;
my $read_id;
my $seq;
my $sign;
my $qual;
my $seq_len;

GetOptions(

           'contaminant:s'   => \$contaminant # it's necessary to define the kind of variable :s means string. but perl automatically turns this into a number if its required
           
);

if (exists $ARGV[0] ){

    $contaminant=$ARGV[0];

} 

# reading STDIN
while(<STDIN>)
{
 chomp;
 $count++;
 
      if ($count == 1){
          $read_id=$_;
      }
      if ($count == 2) {
          $seq=$_;
      }
     if ($count == 3) {
         $sign=$_;
     }
     if ($count == 4) {
         $qual=$_;
         if ($read_id=~/$contaminant/)
         {
           print"$read_id\n$seq\n$sign\n$qual\n";
         }
         else{
          #do nuthin'
         }
      $count=0;
    }
    

}


exit;
