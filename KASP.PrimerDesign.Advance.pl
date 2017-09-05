use strict;
use warnings;
use POSIX qw(tmpnam);
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use List::MoreUtils qw(uniq);
use Time::localtime;
use List::Util qw( min max );

## ======================================
## Usage: see -h
## ======================================

sub usage{
  warn <<END;
  Usage:
  Run by typing: perl KASP.PrimerDesign.Advance.pl -cfg [Parameter config file] -input [Flanking Seq fasta file] -out [Primer output] 
    Required params:
	-c|cfg							[s]	Parameter config file
	-i|input						[s]	Flanking Seq fasta file
	-o|out							[s]	Primer output
    Example: perl KASP.PrimerDesign.Advance.pl -cfg [Parameter config file] -input [Flanking Seq fasta file] -out [Primer output] 
END
  exit;
}
## ======================================
## Get options
## ======================================

my %opt;
%opt = (
	'help'				=> undef,
	'debug'				=> undef,
	'cfg'				=> undef,
	'input'				=> undef,
	'out'				=> undef
);

die usage() if @ARGV == 0;
GetOptions (
  'h|help'				=> \$opt{help},
  'debug'				=> \$opt{debug},
  'c|cfg=s'				=> \$opt{cfg},
  'i|input=s'			=> \$opt{input},
  'o|out=s'				=> \$opt{out}
) or die usage();

#check input paramaters
die usage() if $opt{help};
die usage() unless ( $opt{cfg} );
die usage() unless ( $opt{input} );
die usage() unless ( $opt{out} );

########
#Main Function
########
my %complements;
$complements{"A"} = "T";
$complements{"C"} = "G";
$complements{"G"} = "C";
$complements{"T"} = "A";
$complements{"N"} = "N";
$complements{"a"} = "t";
$complements{"c"} = "g";
$complements{"g"} = "c";
$complements{"t"} = "a";
$complements{"n"} = "n"; 

if (-e "Primer3.Log.txt"){
		system "rm Primer3.Log.txt";
}
if (-e "Blast.log.txt"){
		system "rm Blast.log.txt";
}
if (-e $opt{out}){
		system "rm $opt{out}";
}
if (-e "tmp.dir"){
	system "rm -r tmp.dir";
}
system "mkdir tmp.dir";
my %params = &CFG_Parser($opt{cfg});
print "$opt{cfg} has been read! \n";
#print $params{'max_dist_for_sequencing'}."\n";
my %seq = &Seq_Parser($opt{input});
print "$opt{input} has been parsed! \n";
open(VCF,$opt{input});
my @vcf; my %VCFs;
my $n=0;
while(<VCF>){
	push @vcf,$_;
	next if $_ =~ /^#/;
	chomp($_);$_=~s/\r//;
	$n++;
	my @g = split /\t/,$_;
	my $subchr = $g[1];
	my $pos = $g[2];
	$VCFs{$n}{"Chr"} = $subchr;
	$VCFs{$n}{"Pos"} = $pos;
}
close VCF;

foreach my $snpid (sort keys %seq){
	print "                      \n";
	print "                      \n";
	print "########################################################\n";
	print "           Designing Primer for $snpid                  \n";
	print "########################################################\n";
	my @alleles = split /\//,$seq{$snpid}{"Alleles"};
	my $allele1 = $alleles[0];
	my $allele2 = $alleles[1];
	my $Left_seq = $seq{$snpid}{"Left_Seq"};
	my $Right_seq = $seq{$snpid}{"Right_Seq"};
	my $input1 = $Left_seq.$allele1.$Right_seq; ### Sequence with one SNP allele 
	my $input2 = $Left_seq.$allele2.$Right_seq; ### Sequence with the other SNP allele
	(my $hout1_forward,my $hout1_reverse) = &Primer3_Run($snpid,$input1);
	(my $hout2_forward,my $hout2_reverse) = &Primer3_Run($snpid,$input2);
	my @forwards = uniq(@{$hout1_forward},@{$hout2_forward});
	my @reverses = uniq(@{$hout1_reverse},@{$hout2_reverse});
	#print "Step1: \n$forwards[0]\n$reverses[0]\n";
	die "No Forward Primers can be designed by Primer3! \n" if (!@forwards);
	die "No Reverse Primers can be designed by Primer3! \n" if (!@reverses);
	print "Forward and Reverse Primers Candidates has been designed by Primer3! \n";
	###Find Competitive Primers###
	my $len_flanking = $params{"len_flanking"};
	(my $hcompetitive_forwards,my $hcommon_forwards) = &find_competitive(\@forwards,$len_flanking,"forward"); ### Find the primer with SNP allele in 3'end in positive direction.
	(my $hcompetitive_reverses,my $hcommon_reverses) = &find_competitive(\@reverses,$len_flanking,"reverse"); ### Find the primer with SNP allele in 3'end in reverse direction.
	#my @p1 = @{$hcompetitive_forwards};my @p2 = @{$hcommon_forwards};
	#print "Step2: \n$p1[0]\n";
	#print "$p2[0]\n";
	print "No Competitive Forward Primers can be designed by Primer3! \n" if (!@{$hcompetitive_forwards});
	print "No Competitive Reverse Primers can be designed by Primer3! \n" if (!@{$hcompetitive_reverses});
	print "No Common Forward Primers can be designed by Primer3! \n" if (!@{$hcommon_forwards});
	print "No Common Reverse Primers can be designed by Primer3! \n" if (!@{$hcommon_reverses});
	if ( (!@{$hcompetitive_forwards}) and (!@{$hcompetitive_reverses}) ){
		print "No Competitive Primers can be designed by Primer3! \n";
		next;
	}else {
		print "Competitive Primers have been designed by Primer3! \n";
	};
	 if ((!@{$hcommon_forwards}) and (!@{$hcommon_reverses})){
		print "No Common Primers can be designed by Primer3! \n";
		next;
	 }else {
		print "Common Primers have been designed by Primer3! \n";
	 };
	####Find Competitive Primer Pairs###
	my $hcompetitive_forwards_pairs;my $hcompetitive_reverses_pairs;
	$hcompetitive_forwards_pairs = &find_competitive_pairs($hcompetitive_forwards,"forward"); ### Find Forward Competitive Primer Pair with different SNP slleles.
	$hcompetitive_reverses_pairs = &find_competitive_pairs($hcompetitive_reverses,"reverse"); ### Find Reverse Competitive Primer Pair with different SNP slleles.
	#my @p1 = @{$hcompetitive_forwards_pairs};my @p2 = @{$hcompetitive_reverses_pairs};
	#print "Step3: \n$p1[0]\n$p2[0]\n";
	print "No Competitive Forward Primers can be paired by Primer3! \n" if (!@{$hcompetitive_forwards_pairs});
	print "No Competitive Reverse Primers can be paired by Primer3! \n" if (!@{$hcompetitive_reverses_pairs});
	if ((!@{$hcompetitive_forwards_pairs}) and (!@{$hcompetitive_reverses_pairs})){
		print "No Competitive Primers can be paired by Primer3! \n";
		next;
	} else {
		print "Competitive Primers have been paired by Primer3! \n";
	}
	###Primers Blast###
	my $hcompetitive_forwards_blast_outs;my $hcompetitive_reverses_blast_outs;my $hcommon_forwards_blast_outs;my $hcommon_reverses_blast_outs;
	$hcompetitive_forwards_blast_outs = &qualified_blast_run($hcompetitive_forwards_pairs,"Competitive","forward",$allele1) if (@{$hcompetitive_forwards_pairs}); ### Run Blastn and filter output according to mismatch number and SNP overlap.
	$hcompetitive_reverses_blast_outs = &qualified_blast_run($hcompetitive_reverses_pairs,"Competitive","reverse",$allele1) if (@{$hcompetitive_reverses_pairs}); ### Run Blastn and filter output according to mismatch number and SNP overlap.
	$hcommon_forwards_blast_outs = &qualified_blast_run($hcommon_forwards,"Common") if (@{$hcommon_forwards});
	$hcommon_reverses_blast_outs = &qualified_blast_run($hcommon_reverses,"Common") if (@{$hcommon_reverses});
	#my $hcompetitive_reverses_blast_outs;my $hcommon_forwards_blast_outs;my $hcommon_reverses_blast_outs;
	print "No Competitive Forward Primers Blast outs are qualified! \n" if (!defined($hcompetitive_forwards_blast_outs));
	print "No Competitive Reverse Primers Blast outs are qualified! \n" if (!defined($hcompetitive_reverses_blast_outs));
	print "No Common Forward Primers Blast outs are qualified! \n" if (!defined($hcommon_forwards_blast_outs));
	print "No Common Reverse Primers Blast outs are qualified! \n" if (!defined($hcommon_reverses_blast_outs));
	if (!defined($hcompetitive_forwards_blast_outs) and !defined($hcompetitive_reverses_blast_outs) ){
		print "No Competitive Primers Blast outs are qualified! \n";
		next;
	} else {
		print "Competitive Primers Blast outs are qualified! \n";
	}
	if ( !defined($hcommon_forwards_blast_outs) and !defined($hcommon_reverses_blast_outs) ){
		print "No Common Primers Blast outs are qualified! \n" ;
		next;
	} else {
		print "Common Primers Blast outs are qualified! \n" ;
	}
	###Determine Unique Primer Pairs###
	my %Reverse_Unique_Primer_Pairs;
	if ($hcommon_forwards_blast_outs and $hcompetitive_reverses_blast_outs){
		if ( (%{$hcommon_forwards_blast_outs} and %{$hcompetitive_reverses_blast_outs})){
			%Reverse_Unique_Primer_Pairs = &Unique_Primers($hcommon_forwards_blast_outs, $hcompetitive_reverses_blast_outs,$seq{$snpid}{"Chr"}); ###Filter and output unique Primer Pair 
			print "Reverse Competitive Primer Pairshas been uniqued! \n";
		} 
	}else {
		print "No input for common_forwards_blast and competitive_reverses_blast! \n";
	}
	print "No Unique Competitive Reverse Primers Pairs are found! \n" if (!%Reverse_Unique_Primer_Pairs);
	my %Forward_Unique_Primer_Pairs;
	if ($hcommon_reverses_blast_outs and $hcompetitive_forwards_blast_outs){
		if ((%{$hcompetitive_forwards_blast_outs} and %{$hcommon_reverses_blast_outs}) ){
			%Forward_Unique_Primer_Pairs = &Unique_Primers($hcompetitive_forwards_blast_outs, $hcommon_reverses_blast_outs,$seq{$snpid}{"Chr"}) ;
			print "Forward Competitive Primer_Pairshas been uniqued! \n";
		} 
	}else {
		print "No input for competitive_forwards_blast and common_reverses_blast! \n";
	} 
	print "No Unique Competitive Forward Primers Pairs are found! \n" if (!%Forward_Unique_Primer_Pairs);
	
	if ((!%Forward_Unique_Primer_Pairs) and (!%Reverse_Unique_Primer_Pairs) ){
		print "No Unique Primers Pairs are found! \n";
		next;
	}
	open(OUT,">>$opt{out}");
	print "Outputing .... >> $opt{out} \n";
	print OUT "##The KASP Primers of $snpid \n";
	print OUT "Forward_Primer\tStart\tEnd\tGC\tTm\tReverse_Primer\tStart\tEnd\tGC\tTm\tProduct_Size\tCompetitive_Side\n";
	if (%Forward_Unique_Primer_Pairs){
		my $i=0;
		foreach my $id1 (keys %Forward_Unique_Primer_Pairs){
			my @infos1 = &Extract_Print_Info($id1,$seq{$snpid}{"Pos"},$len_flanking);
			foreach my $id2 (@{$Forward_Unique_Primer_Pairs{$id1}}){
				my @infos2 = &Extract_Print_Info($id2,$seq{$snpid}{"Pos"},$len_flanking);
				foreach my $info2 (@infos2){
					foreach my $info1 (@infos1){
						my @g1 = split /\t/,$info1;
						my @g2 = split /\t/,$info2;
						my $product_size = $g2[2] - $g1[1] + 1;
						print OUT $info1."\t".$info2."\t".$product_size."\t"."Forward"."\n";
						$i++;
					}
				}
			}
		}
		my $num1 = $i/2;
		print "$num1 Primer Pairs with Forward Competitive Primer has been designed! \n";
	}
	
	if (%Reverse_Unique_Primer_Pairs){
		my $j=0;
		foreach my $id1 (keys %Reverse_Unique_Primer_Pairs){
			my @infos1 = &Extract_Print_Info($id1,$seq{$snpid}{"Pos"},$len_flanking);
			foreach my $id2 (@{$Reverse_Unique_Primer_Pairs{$id1}}){
				my @infos2 = &Extract_Print_Info($id2,$seq{$snpid}{"Pos"},$len_flanking);
				foreach my $info1 (@infos1){
					foreach my $info2 (@infos2){
						my @g1 = split /\t/,$info1;
						my @g2 = split /\t/,$info2;
						my $product_size = $g2[2] - $g1[1] + 1;
						print OUT $info1."\t".$info2."\t".$product_size."\t"."Reverse"."\n";
						$j++;
					}
				}
			}
		}
		my $num2 = $j/2;
		print "$num2 Primer Pairs with Reverse Competitive Primer has been designed! \n";
	}
	close OUT;
}









########
#Subfunction
########


sub CFG_Parser{
	#Parse the configure file and obtain parameters for Primer3 and Blastn.
	my $cfgfile = shift;
	open(IN,$cfgfile);
	my %params;
	while(<IN>){
		chomp;s/\r//;
		next unless /\=/;
		my @g = split /\s*\=\s*/,$_;
		$params{$g[0]} = $g[1];
	}
	return %params;
}

sub Seq_Parser{
	#Parse the input sequence fasta file and obtain the Left, Right flanking sequene and SNP alleles.
	my $fa = shift;
	my %seq;
	open(FA,$fa);
	while(<FA>){
		chomp;s/\r//;
		my @g = split /\t/,$_;
		my $snp_id = $g[0];
		my $chr = $g[1];
		my $pos = $g[2];
		my $seq = $g[3];
		my $left_seq; my $right_seq; my $allels;
		if ($seq =~ /(^\w+)\[/){
			$left_seq = $1;
		};
		if ($seq =~ /\](\w+$)/){
			$right_seq = $1;
		};
		if ($seq =~ /\[(.+)\]/){
			$allels = $1;
		};
		if ($left_seq && $right_seq && $snp_id && $pos && $chr){
			$seq{$snp_id}{"Left_Seq"} = $left_seq;
			$seq{$snp_id}{"Right_Seq"} = $right_seq;
			$seq{$snp_id}{"Alleles"} = $allels;
			$seq{$snp_id}{"Pos"} = $pos;
			$seq{$snp_id}{"Chr"} = $chr;
		}
	}
	close FA;
	return %seq;
}


sub Primer3_Run {
	#Output Primer3 input file and Run Primer3, and then read output as array.
	my $snpid = shift;
	my $seq = shift;
	die "The length of Right Flanking sequence is smaller than $params{'PRIMER_MIN_SIZE'}! " if (length($seq) < $params{"PRIMER_MIN_SIZE"});
	if (-e "Primer3.tmp.input"){
		system "rm Primer3.tmp.input";
	}
	my $tmp_input = "Primer3.tmp.input";
	my $len = length($seq);
	open (OUT,">$tmp_input");
	print OUT "SEQUENCE_ID"."=".$snpid."\n";
	print OUT "SEQUENCE_TEMPLATE"."=".$seq."\n";
	#print OUT "SEQUENCE_TARGET"."="."0\,0"."\n";
	print OUT "PRIMER_TASK"."="."generic"."\n";
	print OUT "PRIMER_PICK_LEFT_PRIMER"."="."1"."\n";
	print OUT "PRIMER_PICK_INTERNAL_OLIGO"."="."0"."\n";
	print OUT "PRIMER_PICK_RIGHT_PRIMER"."="."1"."\n";
	print OUT "PRIMER_OPT_SIZE"."=".$params{"PRIMER_OPT_SIZE"}."\n";
	print OUT "PRIMER_MIN_SIZE"."=".$params{"PRIMER_MIN_SIZE"}."\n";
	print OUT "PRIMER_MAX_SIZE"."=".$params{"PRIMER_MAX_SIZE"}."\n";
	print OUT "PRIMER_MAX_TM"."=".$params{"PRIMER_MAX_TM"}."\n";
	print OUT "PRIMER_MIN_TM"."=".$params{"PRIMER_MIN_TM"}."\n";
	print OUT "PRIMER_OPT_TM"."=".$params{"PRIMER_OPT_TM"}."\n";
	print OUT "PRIMER_MAX_GC"."=".$params{"PRIMER_MAX_GC"}."\n";
	print OUT "PRIMER_MIN_GC"."=".$params{"PRIMER_MIN_GC"}."\n";
	print OUT "PRIMER_MAX_NS_ACCEPTED"."="."1"."\n";
	print OUT "PRIMER_PRODUCT_SIZE_RANGE"."=".$params{"min_pcr_product_size"}."-".$params{"max_pcr_product_size"}."\n";
	print OUT "P3_FILE_FLAG"."="."1"."\n";
	print OUT "PRIMER_EXPLAIN_FLAG"."="."1"."\n";
	print OUT "PRIMER_THERMODYNAMIC_PARAMETERS_PATH"."=".$params{"PRIMER_THERMODYNAMIC_PARAMETERS_PATH"}."\n";
	print OUT "=\n";
	close OUT;
	print "The Primer3 input file 'Primer3.tmp.input' for Primer has been generated! \n";
	if (-e "$snpid.for"){
		system "rm $snpid.for";
	}
	if (-e "$snpid.rev"){
		system "rm $snpid.rev";
	}
	system "$params{'primer3'} Primer3.tmp.input > Primer3.Log.txt 2>&1 ";
	print "Primer3 is done! \n";
	my $forward = "$snpid.for";
	my @forward_out;
	my $i=0;
	open(IN,$forward);
	while(<IN>){
		chomp;s/\r//;
		$i++;
		if ($i==1){
			next;
		}
		next if /#/;
		$_ =~ s/\s+/\t/g;
		$_ =~ s/^\t//g;
		my @g = split /\t/,$_;
		shift @g;
		my $seq = $g[0];
		my $line =join "\t", @g;
		push @forward_out, $line;
	}
	close IN;
	
	my $reverse = "$snpid.rev";
	my @reverse_out;
	my $ii=0;
	open(IN,$reverse);
	while(<IN>){
		chomp;s/\r//;
		$ii++;
		if ($ii==1){
			next;
		}
		next if /#/;
		$_ =~ s/\s+/\t/g;
		$_ =~ s/^\t//g;
		my @g = split /\t/,$_;
		shift @g;
		my $start = $g[1]-$g[2]+1;###This is important! For reverse primer the start point is the 5'end in the reverse direction and 3'end in the positive direction. 
		$g[1] = $start;
		my $line =join "\t", @g;
		push @reverse_out, $line;
	}
	close IN;
	if (-e "$snpid.for"){
		system "mv $snpid.for tmp.dir";
	}
	if (-e "$snpid.rev"){
		system "mv $snpid.rev tmp.dir";
	}
	return (\@forward_out,\@reverse_out);
}


sub find_competitive {
	#Find competitive primers according to the SNP position.
	(my $hprimers, my $len_flanking, my $direction) = @_;
	my @out1; my @out2;
	if ($len_flanking < 0) {die "The flanking length is smaller than 0 ! \n"};
	my @primers = @{$hprimers};
	foreach my $v (@primers){
		my @g = split /\t/,$v;
		my $start = $g[1]+1;
		my $end = $g[1]+$g[2];
		if (($direction eq "forward") and ($end == ($len_flanking+1))){
			push @out1, $v;
		} elsif ( ($direction eq "reverse") and ($start == ($len_flanking+1)) ){
			push @out1, $v;
		} else {
			push @out2, $v;
		} 
	}
	return (\@out1,\@out2);
}

sub find_competitive_pairs{
	#Find competitive primers pair according to the SNP alleles.
	(my $hcandidates, my $direction) = @_;
	my @candidates = @{$hcandidates};
	my %news;
	my @outs;
	if ($direction eq "forward"){
		foreach my $v (@candidates){
			my @g = split /\t/,$v;
			$g[0] =~ s/\w$//i;
			push @{$news{$g[0]}}, $v;
		}
	} elsif ($direction eq "reverse"){
		foreach my $v (@candidates){
			my @g = split /\t/,$v;
			$g[0] =~ s/\w$//i;
			push @{$news{$g[0]}}, $v;
		}
	}
	if (%news){
		foreach my $h (keys %news){
			if (scalar(@{$news{$h}}) == 2){
				my $out = join "-",@{$news{$h}};
				push @outs, $out;
			}
		}
	}
	return (\@outs);
}


sub qualified_blast_run {
	#Make Blast query fasta file and Filter Blast Otput according to mismatch number and SNP overlap.
	my $query = shift;
	my $type = shift;
	my $direction = shift;
	my $allele1 = shift;
	my $queryfile = "Blast.input.fa";
	open(IN,">".$queryfile);
	my $out = "tmp.blast.out";
	if (-e "tmp.blast.out"){
		system "rm tmp.blast.out";
	}
	if ($type eq "Competitive"){
		my @competitive_pairs = @{$query};
		foreach my $line (@competitive_pairs){
			my @g = split /\-/,$line;
			my @gg = split /\t/,$g[0];
			if ($direction eq "forward"){
				$gg[0] =~ s/\w$/$allele1/;
			} else {
				my $allele = &Complement($allele1);
				$gg[0] =~ s/\w$/$allele/;
			}
			my $seq = $gg[0];
			my @vs = split /\t/,$line;
			my $id = join "_",@vs;
			print IN ">$id\n$seq\n";
		}
	} elsif ($type eq "Common") {
		my @common_pairs = @{$query};
		foreach my $line (@common_pairs){
			my @g = split /\t/,$line;
			my $id = join "_",@g;
			my $seq = $g[0];
			print IN ">$id\n$seq\n";
		}
	} else {
		die "The Type is wrong! \n";
	}
	close IN;
	system "$params{'blastn'} -db $params{'reference_genome'} -query $queryfile -out $out -evalue $params{'evalue'} -num_threads $params{'num_cpus'} -outfmt \"6 std gaps nident\" -dust no -gapopen 4 -gapextend 2 -penalty -2 -reward 2 -word_size $params{'word_size'} -max_target_seqs 500 > Blast.log.txt 2>&1";
	#print "$params{'blastn'} -db $params{'reference_genome'} -query $queryfile -out $out -evalue $params{'evalue'} -num_threads $params{'num_cpus'} -outfmt \"6 std gaps nident\" -dust no -gapopen 4 -gapextend 2 -penalty -2 -reward 2 -word_size $params{'word_size'} -max_target_seqs 500 > Blast.log.txt 2>&1\n";
	my %outs;
	open (OUT,"$out");
	while(<OUT>){
		chomp;s/\r//;
		my @g = split /\t/,$_;
		my $start = $g[8];
		my $end = $g[9];
		if ($end < $start) {
			$g[8] = $end;
			$g[9] = $start;
		}
		my @ots = &Filter_Blast(@g);
		if (@ots){
			my $id = shift @g;
			my $line = join "\t",@g;
			push @{$outs{$id}}, $line;
		} else {
			next;
		}
	}
	close OUT;
	return (\%outs);
}


sub Filter_Blast {
	#Filter Blast Otput according to mismatch number and SNP overlap
	my ($id, $chr, $a, $b, $c, $d, $e, $f, $pos1, $pos2, $score1, $score2, $g, $h) = @_;
	my @raw = ($id, $chr, $a, $b, $c, $d, $e, $f, $pos1, $pos2, $score1, $score2, $g, $h);
	my @v = split /\_/,$id;
	my $seq = $v[0];
	my $len = length($seq);
	if (($h>$len-$params{"min_mismatches_close3p"}) or ($f<=$len-$params{"min_dist3p"} and $h>$len-$params{"min_mismatches"})){
		my $j = 0;
		foreach my $d (keys %VCFs){
			my $subchr = $VCFs{$d}{"Chr"};
			my $pos = $VCFs{$d}{"Pos"};
			if ((($subchr eq $chr) and ($pos>$pos1) and ($pos<$pos2)) || ( ($chr eq $subchr) and ($pos>$pos2) and ($pos<$pos1) )){
				$j=1;
			}
		}
		if ($j == 0){
			return (@raw);
		} else {
			return (());
		}
	} else {
		return (());
	}
}



sub Unique_Primers {
	#Get Unique Primers and test the Tm difference and Cross-Dimer between Forward and Reverse Primers!.
	(my $hforward, my $hreverse, my $target_chr) = @_;
	my %outs;
	my %forwards = %{$hforward}; my %reverses = %{$hreverse};
	foreach my $id1 (keys %forwards){
		my $hforward_blast_outs = \@{$forwards{$id1}};
		foreach my $id2 (keys %reverses){
			my $hreverse_blast_outs = \@{$reverses{$id2}};
			if (&Unique_Blast($hforward_blast_outs,$hreverse_blast_outs,$target_chr)){
				#print "$id1\n$id2\n";
				my @TMs1 = &Extract_TM($id1);
				my @TMs2 = &Extract_TM($id2);
				my @Seqs1 = &Extract_Seq($id1);
				my @Seqs2 = &Extract_Seq($id2);
				my $max_TM_diff = 0;
				my $Dimer = 0;
				foreach my $TM1 (@TMs1){
					foreach my $TM2 (@TMs2){
						if ( abs($TM2-$TM1) > $max_TM_diff ){
							$max_TM_diff = abs($TM2-$TM1);
						}
					}
				}
				foreach my $seq1 (@Seqs1){
					foreach my $seq2 (@Seqs2){
						if ( &dimer($seq1,$seq2)==1 ){
							$Dimer = 1;
							print "Cross-Dimer between Forward and Reverse Primers exists! $seq1 $seq2\n";
						}
					}
				}
				if (($max_TM_diff < 5) and ($Dimer ==0)){
					push @{$outs{$id1}},$id2;
				} else {
					print "TM difference between Left Primer and Right Primer is larger than 5! \n";
					#print "TM difference between Left Primer and Right Primer is larger than 5! \n $id1\n$id2\n";
				}
			}
		}
	}
	return (%outs);
}

sub Unique_Blast {
	#Get the primers pairs with unique relative position in the whole genome.
	(my $hforward_blast_outs, my $hreverse_blast_outs,my $target_chr) = @_;
	my @forward_blast_outs = @{$hforward_blast_outs};
	my @reverse_blast_outs = @{$hreverse_blast_outs};
	my %test_unique;
	my %best_unique;
	foreach my $out1 (@forward_blast_outs){
		my ($chr1, $a1, $b1, $c1, $d1, $e1, $f1, $pos11, $pos12, $score11, $score12, $g1, $h1) = split /\t/,$out1;
		foreach my $out2 (@reverse_blast_outs ){
			my ($chr2, $a2, $b2, $c2, $d2, $e2, $f2, $pos21, $pos22, $score21, $score22, $g2, $h2) = split /\t/,$out2;
			my $max_dis = ($pos22-$pos11+1);
			if (($chr1 eq $chr2) and ($max_dis <= $params{"max_dist_for_sequencing"}) and ($max_dis > 0)){
				$test_unique{$chr2}++;
				if ( ( $max_dis <= $params{'max_pcr_product_size'}) and ($max_dis >= $params{'min_pcr_product_size'}) ){
					$best_unique{$chr2}++;
				}
			}
		}
	}
	if ( (scalar(keys %test_unique) == 1) and (scalar(keys %best_unique) == 1) and ($test_unique{$target_chr} == 1) and ($best_unique{$target_chr} == 1)){
		return 1;
	} else {
		return 0;
	}
}

sub Extract_TM {
	#Extract Tm values from the Primer3 results.
	my $id = shift;
	my @TMs;
	if ($id =~ /\-/){
		my @g1 = split /\-/,$id;
		foreach my $v (@g1){
			my @g2 = split /\_/,$v;
			push @TMs, $g2[5];
		}
	} else {
		my @g2 = split /\_/,$id;
		push @TMs, $g2[5];
	}
	return (@TMs);
}

sub Extract_Seq {
	#Extract Sequence from the Primer3 results.
	my $id = shift;
	my @seqs;
	if ($id =~ /\-/){
		my @g1 = split /\-/,$id;
		foreach my $v (@g1){
			my @g2 = split /\_/,$v;
			push @seqs, $g2[0];
		}
	} else {
		my @g2 = split /\_/,$id;
		push @seqs, $g2[0];
	}
	return (@seqs);
}


sub Extract_Print_Info {
	#Extract Sequence, position, GC content and Tm values from the Primer3 results. It can distinguish the forward and reverse primer.
	(my $id, my $snp_pos, my $len_flanking) = @_;
	my @outs;
	if ($id =~ /\-/){
		my @g1 = split /\-/,$id;
		foreach my $v (@g1){
			my @g2 = split /\_/,$v;
			my $seq = $g2[0];
			my $start = $snp_pos - $len_flanking + $g2[1];
			my $end = $start +$g2[2] -1;
			my $GC = $g2[4];
			my $TM = $g2[5];
			my $line = join "\t", ($seq,$start,$end,$GC,$TM);
			push @outs, $line;
		}
	} else {
		my @g2 = split /\_/,$id;
		my $seq = $g2[0];
		my $start = $snp_pos - $len_flanking + $g2[1];
		my $end = $start +$g2[2] -1;
		my $GC = $g2[4];
		my $TM = $g2[5];
		my $line = join "\t", ($seq,$start,$end,$GC,$TM);
		push @outs, $line;
	}
	return (@outs);
}

sub Complement {
	#Return the complemental sequence. 
	my $seq = shift;
	my @g = split //,$seq;
	my @out;
	foreach my $s (@g){
		push @out, $complements{$s};
	}
	my $line = join "",@out;
	return $line;
}

sub dimer {
	#Dimers: A primer self-dimer is formed by intermolecular interactions between two of the same primers, 
	#i.e., forward vs. forward or reverse vs. reverse. Cross-dimers are formed by intermolecular interactions
	#between a forward and a reverse primer. Primer dimers with more than 4 consecutive base pairings are not passed,
	#and will not be considered as potential primer pairs.
	(my $seq1, my $seq2) = @_;
	my $total_len1 = length($seq1);
	my $total_len2 = length($seq2);
	my $index = 0;
	my $reverse_complemennt_3end_seq1 = &Reverse_Complement(substr($seq1,($total_len1-1-4),5));
	my $end_seq2 = substr($seq2,($total_len2-1-4),5);
	if ($reverse_complemennt_3end_seq1 eq $end_seq2){
		$index=1;
	}
	if ($index==1){
		return 1;
	} else {
		return 0;
	}
}

sub Reverse_Complement {
	my $seq = shift;
	my @g = split //,$seq;
	my @out;
	foreach my $s (reverse @g){
		push @out, $complements{$s};
	}
	my $line = join "",@out;
	return $line;
}
