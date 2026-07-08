#!/usr/bin/env perl
use strict;
use lib "/home/groups/assembly/talioto/myperlmods";
use lib "/home/groups/assembly/talioto/myperlmods/Bio";
use SeqOp;
use Getopt::Long;
use File::Basename qw( fileparse );
my $fname = 0;
my $num_n = 1;
my $species = '';
my $taxid = '';
my $assembly_name = 0;
my $outfilebase = 0;
my $prefix = '';
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my @abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my $realyear = 1900+$year;
my $date = "$mday-$abbr[$mon]-$realyear";
my $center = "CNAG";
my $description = "";
my $min_scaff_size = 20;
my $min_contig_size = 20;
my $sp = 0;
my $cp = 0;
my $lookuptable = 1;
my $keep = 1;
my $sort = 1;
my $gap_cutoff = 0;
my $rename_tsv  = 0;
GetOptions(
	   'f|s|i:s'        => \$fname,
	   'k|keep!'        => \$keep,
	   'ns:s'           => \$num_n,
	   'name:s'         => \$assembly_name,
	   'taxid:s'        => \$taxid,
	   'species:s'      => \$species,
	   'center:s'       => \$center,
	   'date:s'         => \$date,
	   'mins:s'         => \$min_scaff_size,
	   'minc:s'         => \$min_contig_size,
	   'sp:s'           => \$sp,
	   'cp:s'           => \$cp,
	   'lookup!'        => \$lookuptable,
	   'sort!'          => \$sort,
     'g|gap:i'        => \$gap_cutoff,
     'rename:s'       => \$rename_tsv
	);
my $options = <<HERE;
           -f,-s,-i	<file.fasta>  name of fasta file, required
           -n	<int>	              number of Ns on which to break (1)
           -g	<int>	              threshold unestimated gap size  (0)
           -name <assembly_name>  assembly name used for output filename and seq prefix
           -taxid <int>           taxonomy id from NCBI
           -species <species>     species
           -center <center>	      center (CNAG)
           -date <date>	          date of assembly (current local time)
           -mins <int>	          minimum scaffold size to keep (20)
           -minc <int>	          minimum contig size to keep (20)
           -sort                  sort scaffolds by size longest to shortest (default)
           -nosort                don't sort/original order
           -rename <newnames.tsv> tab-delimited text file in format:
                                    old_scaff_name <tab> new_scaff_name
HERE
die "usage: $0 -f <scaffolds.fa> [options]\n$options\n" if !$fname;
print STDERR "'-n 0' does not make sense. Running with '-n 1'\n" if !$num_n;
my ($base,$path,$ext) = fileparse($fname,qw(\.fa \.fasta));
if ($assembly_name){
  $prefix = $assembly_name.'_';
	$outfilebase = $assembly_name;
}else{
  $prefix = '';
	$outfilebase = $base;
}
my $cout = "$outfilebase.contigs.fa";
my $sout = "$outfilebase.scaffolds.fa";
my $offsets = "$outfilebase.contigs.offsets";
my $soffsets = "$outfilebase.scaffolds.offsets";
my $gaps = "$outfilebase.scaffolds.gaps.bed";
my $agp = "$outfilebase.agp";
my $lup = "$outfilebase.lookup.txt";
my @scaffolds;

my %newnames = ();
if ($rename_tsv){
	open (TSV,"<$rename_tsv") or die "Couldn't open $rename_tsv!";
	while(my $l = <TSV>){
		chomp $l;
		my @f = split "\t",$l;
		$newnames{$f[0]}=$f[1];
	}
	close TSV;
}

open IN, "FastaToTbl $fname |";
my $maxc=1;
my $nums=0;
my $assembly_span = 0;
print STDERR "Reading $fname...\n";
while (<IN>) {
  $nums++;
  my $c=1;
  my @F=split;
	my $slength = length($F[1]);
	$assembly_span += $slength;
  push @scaffolds,{'len'=>$slength,'id'=>$F[0],'seq'=>$F[1]};
    while ($F[1]=~/N+/g) {
      $c++;
    }
  if ($c>$maxc) {
    $maxc=$c;
  }
}
close IN;
if (!$sp) {
  $sp = length($nums);
}
if (!$cp) {
  $cp = length($maxc);
}
my $progress_chunksize = int($assembly_span/20);
#open IN, "FastaToTbl $fname |";
unlink $sout if -e $sout;
unlink $cout if -e $cout;
open (AGP,">$agp");
print AGP '##agp-version 2.1',"\n";
print AGP '# ',"ORGANISM: $species\n";
print AGP '# ',"TAX_ID: $taxid\n";
print AGP '# ',"ASSEMBLY NAME: $assembly_name\n";
print AGP '# ',"ASSEMBLY DATE: $date\n";
print AGP '# ',"GENOME CENTER: $center\n";
print AGP '# ',"DESCRIPTION: $description\n";


open SOUT, "| TblToFasta >> $sout";
open COUT, "| TblToFasta >> $cout";
open OFF, ">$offsets";
open SOFF, ">$soffsets";
open GAP, ">$gaps";
if ($lookuptable){open LU,">$lup";}
my $scount = 1;
print STDERR "Processing scaffolds\n";
#while (<IN>){
my @sorted_scaffolds=();
if ($sort){@sorted_scaffolds = sort lensort @scaffolds;@scaffolds=@sorted_scaffolds;}
my $cum_length = 0;
my $thresh = $progress_chunksize;
foreach my $scaff (@scaffolds){
  chomp;
  my $sid = $scaff->{id};
  my $seq = $scaff->{seq};
	next if $seq!~/[ACGTacgt]/;
  my $seq_length = $scaff->{len};
	### Progress bar
	$cum_length += $seq_length;
	if ($cum_length>=$thresh){
		print STDERR "." x int(($cum_length-$thresh)/$progress_chunksize);
		$thresh = int($cum_length/$progress_chunksize)+1;
	}
  my $oldid = $sid;
  #my $seq_up = uc($seq);
  #$seq = $seq_up;
  #print STDERR "$sid\n";
  my $num = $min_contig_size - 1;
  if ($num>0){
  	$seq=~s/^([^N]{1,$num})N/"N" x (length($1)+1)/ge;
  	$seq=~s/N([^N]{1,$num})N/"N" x (length($1)+2)/ge;
  	$seq=~s/N([^N]{1,$num})$/"N" x (length($1)+1)/ge;
  }
  if ($seq=~s/^(N+)//){ ### this creates new  offset
    my $scaffoffset = length($1);
    print SOFF "$oldid\t$scaffoffset\n";
  }else{
    print SOFF "$oldid\t0\n";
  }
  $seq=~s/N+$//;
  next if length($seq)<$min_scaff_size;

  #my $outseq = $seq;
  #$outseq=~s/[^ACGTacgt]/N/g;
  if($rename_tsv){
		$keep = 1;
		my $suffix = '';
		if($sid =~s/(_unloc)$//){
				$suffix = '_unloc';
		}
		if (exists $newnames{$sid}){
			$sid = $newnames{$sid};
		}
		$sid = $sid.$suffix;
	}
	if ($sid=~/SUPER_/){
      $sid=~s/SUPER_/Chr/;
  }
	if($keep){
		$sid = $prefix.$sid;
	}else{
		$sid = sprintf("%s"."s%0$sp"."i",$prefix,$scount++);
	}



  if($lookuptable){print LU "$sid\t$oldid\n";}
  my $count = 1;
  # my @contigs = split /N+/,$seq;
  # foreach my $contig (@contigs){

  #   my $cid = sprintf("%s.contig%04d",$sid,$count++);
  #   print COUT "$cid\t$contig\n";
  # }
  my $END = 0;
  my $offset = 0;
  my $element = 0;
  my $part_number = 1;
  my $evidence = 'paired-ends';
  my $component_type = "N";
	#print SOUT "$sid\t$seq\n";
	print SOUT "$sid\t";
  foreach my $contig_or_gap (split(/(N+)/i,$seq)){# foreach my $contig_or_gap (split(/(N{$num_n,})/i,$seq))
    #my $contig = $1;
    if ($element % 2){
      my $gapsize = length($contig_or_gap);
      if ($gapsize >100){
				$gapsize = 100;
				$evidence = "proximity_ligation";$component_type="U";
      }elsif($gapsize > $gap_cutoff){
				$gapsize = 100;
				$evidence = "paired-ends";$component_type="U";
      }else{
				$evidence = "paired-ends";$component_type="N";
      }
			print SOUT "N" x $gapsize;
      print AGP join("\t",($sid,$offset+1,$offset+$gapsize, $part_number++,$component_type,$gapsize,'scaffold','yes',$evidence)),"\n";
      $offset+=$gapsize;
    }else{
      #$contig_or_gap=~s/[^ACGTacgt]/N/g; #replace IUPAC ambiguity codes;
      #my $offset = $end - length($contig);
      print GAP  "$sid\t$END\t$offset\n" if $count > 1;;
      my $cid = sprintf("%s_c%0$cp"."d",$sid,$count++);
      my $contig_len = length($contig_or_gap);
      print COUT "$cid\t$contig_or_gap\n";
			print SOUT "$contig_or_gap";
      print OFF  "$cid\t$offset\n";
      print AGP join("\t",($sid,$offset+1,$offset + $contig_len, $part_number++,'W',$cid,'1',$contig_len,'+')),"\n";
      $offset+=$contig_len;
      $END = $offset;
    }
    $element++;
  }
	print SOUT "\n";
}
#close IN;
close COUT;
close SOUT;
close GAP;
close AGP;
close OFF;
if ($lookuptable){close LU;}
print STDERR "\nFinished OK.\n";
sub lensort
  {
    $b->{len} <=> $a->{len}
  }
