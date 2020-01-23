#!/usr/bin/perl
die  "Sorry\!You have (some) ARGV wrong\!\n" unless (@ARGV ==2);
my $idt=0.8;
my %basepair = ('A'=>'T','C'=>'G','G'=>'C','T'=>'A'); 
my $file=$ARGV[0];
$file=~s/\.fastq$//g;$file=~s/\.fq$//g;
my $file2=$ARGV[1];


sub levenshtein($$){
  my @A=split //, lc shift;
  my @B=split //, lc shift;
  my @W=(0..@B);
  my ($i, $j, $cur, $next);
  for $i (0..$#A){
 $cur=$i+1;
 for $j (0..$#B){
  $next=min(
   $W[$j+1]+1,
   $cur+1,
   ($A[$i] ne $B[$j])+$W[$j]
  );
  $W[$j]=$cur;
  $cur=$next;
 }
 $W[@B]=$next;
  }
  return $next;
}
sub min($$$){
  if ($_[0] < $_[2]){ pop @_; } else { shift @_; }
  return $_[0] < $_[1]? $_[0]:$_[1];
}

sub leng   #当2个数据进行模糊匹配的数据，用最长字符的数据底数相除
{
my ($s1,$s2) = @_;
    my $len1 = length $s1;
     my $len2 = length $s2;
  my $ss1=$len1>$len2?$len1:$len2;
   my $ss2=$len1<$len2?$len1:$len2;
   if ($ss2/$ss1>$idt){
	return ($s1,$s2)
   }

}
sub leng11
{   my ($s1,$s2) = @_;
    my $len1 = length $s1;
     my $len2 = length $s2;
  return $len1>$len2?$len1:$len2

}


my $file=$ARGV[0];
open IN1,"$file2" or die "$!";  ## input target file
open TAR1,"${file}_fullytarget20.txt" or die "$!";    ## output from step1_Get_sgRNA_fromPool.pl
open TAR3,"${file}_fullytarget20_negS.txt" or die "$!";  ## output from step1_Get_sgRNA_fromPool.pl

open OUT,">${file}_fullymatch.txt";
open OUT2,">${file}_mismatch.txt";
open OUT3,">${file}_nomatch.txt";

my (%data,%sgR,@ip,%count,%sgC);
while(<IN1>){
	chomp;s/\r//g;
	my @i=split("\t");
	#next unless $i[1]=~/UP/;
	$i[1]=~s/^[A|G|C|T]//g;
	#$i[1]=~s/SG-//g;$i[1]=~s/-UP//g;
	$data{$i[0]}=$i[1];
}

my $CT=0;
while(<TAR1>){
	chomp;s/\r//g;
	$sgR{$_}++;
}
$CT=keys %sgR;
print "$CT counted\n";


while(<TAR3>){
	chomp;s/\r//g;
	my @i=split("\t");
	my @rever_seq = ();
	my $sg2=$i[1];
	while($sg2){
		my $last = chop $sg2;
		push @rever_seq,$basepair{$last};
	}
	my $reverSg=join("",@rever_seq);
	$sgR{join("\t",$i[0],$reverSg)}++;
}
$CT=keys %sgR;
print "$CT counted\n";


my $CT2=0;
foreach my $m (sort keys %data){
	foreach my $n (sort keys %sgR){
		my @j=split("\t",$n);
		next unless ($j[1] eq $data{$m});
		$count{$m}++;$sgC{$j[1]}++;
		#print "perfect.. $m\n";
	}
}
$CT2=keys %count;
print "Equally is $CT2\n";

foreach my $p (sort keys %data){
	$count{$p}=0 if not exists $count{$p};
	print OUT "$p\t$data{$p}\t$count{$p}\n";
}


foreach my $h (sort keys %sgR){
	my @j=split("\t",$h);
	next if exists $sgC{$j[1]};
	@ip=(@ip,$j[1]);
}


for(my $i=0;$i<=$#ip;$i++){
	foreach my $j (sort keys %data){
	my $ID=$j;
	my $s2=$data{$j};
	if (leng($ip[$i],$s2)){
		if ((levenshtein($ip[$i],$s2)/leng11($ip[$i],$s2))<(1-$idt)){
			my $OBS_IDT=1-levenshtein($ip[$i],$s2)/leng11($ip[$i],$s2);
				print OUT2 "$j\t$i\t$ip[$i]\t$OBS_IDT\t$s2\n";
				$sgC{$ip[$i]}++;
		}
	}
	}
	#print "left.. $i\n";
}


foreach my $h (sort keys %sgR){
	my @j=split("\t",$h);
	next if exists $sgC{$j[1]};
	print OUT3 "$h\n";
}

close IN1;
close OUT;
