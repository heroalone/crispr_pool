#!/usr/bin/perl

die  "Sorry\!You have (some) ARGV wrong\!\n" unless (@ARGV ==1);
my $infilename=$ARGV[0];my $infilename2=$infilename;

my $LEFT="GCACCGAATTG";
my $RIGHT="GTTTTAGAGC";
%basepair = (A=>T,C=>G,G=>C,T=>A); 
my @reverleft = ();my @reverright = ();
my $LEFT0=$LEFT;my $RIGHT0=$RIGHT;
while($RIGHT0){
	my $last = chop $RIGHT0;
	push @reverleft,$basepair{$last};
}
my $revleft=join("",@reverleft);
while($LEFT0){
	my $last = chop $LEFT0;
	push @reverright,$basepair{$last};
}
my $revright=join("",@reverright);

################# for watson (+) strand #################

open(FQ,"$infilename") or die ("can not open $infilename\n");
$infilename=~s/\.fastq$//g;$infilename=~s/\.fq$//g;
open(OUT1,">${infilename}_fullytarget20_EXT20bp.txt") || die ;
open(OUT2,">${infilename}_fullytarget20.txt") || die ;
open(OUT3,">${infilename}_fullytarget19_21EXT20bp.txt") || die ;
open(OUT4,">${infilename}_fullytarget19_21.txt") || die ;
open(OUT5,">${infilename}_partial.txt") || die ;


my $T=9;my $ID;
while(<FQ>){
	if (/^@/){
		$T=9;
		chomp;
		$ID=$_;
		next;
	}
	if ($T==9){
		chomp;
		my $seq=$_;
		my $seq1=$_;
		my $seq2=$_;
		my $seq3=$_;
		my $seq4=$_;
		my $seq5=$_;
		my $seq6=$_;
		my $length=length($_);
	if ($seq=~/$LEFT/ig && $seq=~/$RIGHT/ig){
		if ($seq1=~/$LEFT[A|G|C|T|N]{19}$RIGHT/ig){
			while($seq2=~/$LEFT([A|G|C|T|N]{19})$RIGHT/ig){
				print OUT2 "$ID\t$1\n";
				my $G20="$LEFT"."$1"."$RIGHT";
				print OUT1 "$ID\t$G20\n";
			}
		}
		if ($seq3=~/$LEFT[A|G|C|T|N]{18,20}$RIGHT/ig){
			while($seq4=~/$LEFT([A|G|C|T|N]{18,20})$RIGHT/ig){
				next if length($1)==19;
				print OUT4 "$ID\t$1\n";
				my $G20="$LEFT"."$1"."$RIGHT";
				print OUT3 "$ID\t$G20\n";
			}
			next;
		}
	}
	elsif($seq=~/$LEFT/ig or $seq=~/$RIGHT/ig){
		if ($seq5=~/$LEFT/ig){
			while(/$LEFT/ig){
				my $end = pos();
				next if ($length-$end)<5;
				print OUT5 "$ID\t$seq\n";
			}
		}
		if ($seq6=~/$RIGHT/ig){
			while(/$RIGHT/ig){
				my $end = pos();
				next if ($end)<15;
				print OUT5 "$ID\t$seq\n";
			}
			next;
		}
	}
	}
	$T=0;
}

#print "\n";
close FQ;close OUT2;close OUT1;close OUT3;close OUT4;close OUT5;

################# for crick (-) strand #################

open(FQQ,"$infilename2") or die ("can not open $infilename2\n");
open(OUTT1,">${infilename}_fullytarget20_EXT20bp_negS.txt") || die ;
open(OUTT2,">${infilename}_fullytarget20_negS.txt") || die ;
open(OUTT3,">${infilename}_fullytarget19_21EXT20bp_negS.txt") || die ;
open(OUTT4,">${infilename}_fullytarget19_21_negS.txt") || die ;
open(OUTT5,">${infilename}_partial_negS.txt") || die ;


$T=9;
while(<FQQ>){
	if (/^@/){
		$T=9;
		chomp;
		$ID=$_;
		next;
	}
	if ($T==9){
		chomp;
		my $seq=$_;
		my $seq1=$_;
		my $seq2=$_;
		my $seq3=$_;
		my $seq4=$_;
		my $seq5=$_;
		my $seq6=$_;
		my $length=length($_);
	next if ($seq=~/$LEFT/ig or $seq=~/$RIGHT/ig);
		if ($seq=~/$revleft/ig && $seq=~/$revright/ig){
		if ($seq1=~/$revleft[A|G|C|T|N]{19}$revright/ig){
			while($seq2=~/$revleft([A|G|C|T|N]{19})$revright/ig){
				print OUTT2 "$ID\t$1\n";
				my $G20="$revleft"."$1"."$revright";
				print OUTT1 "$ID\t$G20\n";
			}
		}
		if ($seq3=~/$revleft[A|G|C|T|N]{18,20}$revright/ig){
			while($seq4=~/$revleft([A|G|C|T|N]{18,20})$revright/ig){
				next if length($1)==19;
				print OUTT4 "$ID\t$1\n";
				my $G20="$revleft"."$1"."$revright";
				print OUTT3 "$ID\t$G20\n";
			}
		}
	}
	
	elsif($seq=~/$revleft/ig or $seq=~/$revright/ig){
		if ($seq5=~/$revleft/ig){
			while(/$revleft/ig){
				my $end = pos();
				next if ($length-$end)<5;
				print OUTT5 "$ID\t$seq\n";
			}
			next;
		}
		if ($seq6=~/$revright/ig){
			while(/$revright/ig){
				my $end = pos();
				next if ($end)<15;
				print OUTT5 "$ID\t$seq\n";
			}
			next;
		}
	}
	}
	#$T=0;
}

#print "\n";
close FQQ;close OUTT2;close OUTT1;close OUTT3;close OUTT4;close OUTT5;
