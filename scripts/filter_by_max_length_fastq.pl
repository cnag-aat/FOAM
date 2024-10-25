#! /usr/bin/perl

# Program to select reads of certain length or longer from a fastq file
# Fernando Cruz, 26-08-2021

use strict;
use Getopt::Long;
use File::Basename;

# define variables
my $max_len=0;# to enforce giving an integer, see coditional below
my $count=1;
my $read_id;
my $seq;
my $sign;
my $qual;
my $seq_len;

GetOptions(

           'max_len:s'   => \$max_len # it's necessary to define the kind of variable :s means string. but perl automatically turns this into a number if its required
           
);

if ($max_len == 0) { 
die "aborted please specify, maximum read length\!\nsyntax: zcat -f fastq_infile \| filter_by_max_length_fastq.pl \-max_len maximum_read_length\n";

}


# reading STDIN
while(<STDIN>)
{
 chomp;
 $count++;
 

 if ($count == 2){
     $read_id=$_;
 }
 if ($count == 3) {
     $seq=$_;
     $seq_len=(length($seq));
    # print "seq_len $seq_len max_len $max_len\n";
 }
 if ($count == 4) {
     $sign=$_;
 }
 if ($count == 5) {
     $qual=$_;
     
     if($seq_len <= $max_len) {
       print"$read_id\n$seq\n$sign\n$qual\n";
      }
     else{
     #do nuthin'
     }
      $count=1;
 }
    

}


exit;
