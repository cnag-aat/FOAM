#! /usr/bin/perl

# Program to order and orient mitogenomes based on MITOS annotation. (currently mitos version 2.1.3)
 
# Fernando Cruz, 25-04-2023

use strict;
use Getopt::Long;
use File::Basename;

# define variables
my $fasta;
my $bed;
my $count=0;
my $line;
my $seq_id;
my $dna_seq;
my $totlen;

my $strand="";
my @data=();
my $F_start;
my $F_end;
my $F_strand;

my $begin;
my $end;
my $newseq;
my $rev;

GetOptions(

           'fasta|f:s'   => \$fasta, # fasta file with circular contig or reference mitogenome 
           'bed|b:s'   => \$bed # annotation bedfile
);


my ($basename,$path,$ext) = fileparse($fasta,qw(\.fa));

my $outfile=$basename.".oriented.fa";
print "outfile will be $outfile\n";


# Main Functions
load_fasta();
parse_annotation_bed();
reorient_fasta();

# SUB-ROUTINES

# reading Input FASTA

sub load_fasta {

open (FASTA, "$fasta") or die print "I cannot open file $fasta \n";
{

  while(<FASTA>)
  {
      chomp;
      $line=$_;
      
      if (($line=~/^>/) && ($count==0))
      {   $count++;
          $seq_id=$line;
          $seq_id=~s/\>//g;             
      }
      else{
          unless ($count > 1) { 
	      $dna_seq.=$_;
	  }
      }

  }
  print "$seq_id\t".length($dna_seq)."\n";
  $totlen=length($dna_seq);

  close(FASTA) or die print "cannot close the file $fasta\n";

}#FASTA closed

}# END subroutine 

sub parse_annotation_bed {

    open (BED, "$bed") or die print "I cannot open file $bed \n";
{

  while(<BED>)
  {
      chomp;
      $line=$_;

      if ($line=~/trnF/)
      {   
          @data=split/\t/,$line;
          $F_start=$data[1];
          $F_end=$data[2];
          $F_strand=$data[5];
          print "$line\n";
      }
      else{
      #do nothin'
          }
      

  }#while
  print "$F_start $F_end $F_strand\n";
 

  close(BED) or die print "cannot close the file $bed\n";

}#BED closed

}# END read_annotation subroutine

sub reorient_fasta {
# ensure gene order starting with trnF in +
    if ($F_strand eq "+")
    {
      print "trnF strand is $F_strand\n";  
      #Syntax: substr(string, index, length, replacemen
      $begin=substr($dna_seq,$F_start,($totlen-$F_start));
      $end=substr($dna_seq,0,($F_start));
      $newseq=$begin.$end;
      
      print "newseq length:".length($newseq)."\n";
      open (OUT, ">$outfile") or die print "cannot open $outfile\n";
      print OUT ">".$seq_id."\|fwd\n";
      print OUT "$newseq\n";
      close (OUT);
    }
    else {# negative strand
     print "trnF strand is $F_strand\n";  
     #We will place the trnF gene at the end and then reverse-complement it
      $begin=substr($dna_seq,$F_end,($totlen-$F_end+1));# end coords are always 1-base 
      $end=substr($dna_seq,0,($F_end));# end will become the begin once reversed. It actually goes from 0 to F_end-1
      $newseq=$begin.$end;
      print "newseq length: ".length($newseq)."\n";
     
     #reverse to place trnF in forward sense   
     $newseq=~tr/acgtACGT/tgcaTGCA/;
     $rev=reverse($newseq);

     open (OUT, ">$outfile") or die print "cannot open $outfile\n";
     print OUT ">".$seq_id."\|fwd\n";
     print OUT "$rev\n";
     close (OUT);
   
    }


    
}###
