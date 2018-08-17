#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opt;
GetOptions(\%opt,"trantab");
@ARGV || die"usage: perl $0 <in.table> [topline(10)] > out.sel.tab\n";
my ($intab,$line) = @ARGV;
$line ||= 10;
my %prof=();
open IN,$intab || die$!;
my $ti=<IN>;chomp $ti;
my @sam=split /\t/,$ti;
my %sum=();
shift @sam;
while(<IN>){
    chomp;
    my @l = /\t/ ? split/\t+/ : split;
    my $id=shift @l;
    for my $i(0..$#sam){
		$sum{$id}+=$l[$i];
		$prof{$id}{$sam[$i]}=$l[$i];
	}
}
close IN;
my %data=();
my @fun=();
if ($opt{trantab}){
	my $num=0;
	foreach my $k(sort{$sum{$b}<=>$sum{$a}} keys %sum){
			if ($k=~/Unclassified/){
				foreach my $k1(@sam){
					$data{Unclassified_KO}{$k1}+=$prof{$k}{$k1};
				}
			}else{
				$num+=1;
				if ($num>$line){
					foreach my $k1(@sam){
						$data{"Other"}{$k1}+=$prof{$k}{$k1};
					}
				}else{
					foreach my $k1(@sam){
						$data{$k}{$k1}=$prof{$k}{$k1};
					}
					push @fun,$k;
				}
			}
	}
	push @fun,"Other" if (exists $data{Other});
	push @fun,"Unclassified_KO" if (exists $data{Unclassified_KO});
	foreach my $k2(@fun){
		my $k1=(split /\|/,$k2)[-1];
		$k1=~s/\_/ /g;
		print "\t$k1";
	}
	print "\n";
	foreach my $k4(@sam){
		print "$k4";
		foreach my $k3(@fun){
			print "\t$data{$k3}{$k4}";
		}
		print "\n";
	}
}else{
	foreach my $k1(@sam){
		print "\t$k1";
	}
	print "\n";
	my $num=0;
	foreach my $k(sort{$sum{$b}<=>$sum{$a}} keys %sum){
		if ($k=~/Unclassified/){
			foreach my $k6(@sam){
				$data{Unclassified_KO}{$k6}+=$prof{$k}{$k6};
			}
		}else{
			$num+=1;
		    if ($num >$line){
				foreach my $k2(@sam){
					$data{Other}{$k2}+=$prof{$k}{$k2};
				}
			}else{
				print "$k";
				foreach my $k2(@sam){
					print "\t$prof{$k}{$k2}";
				}
				print "\n";
			}
		}
	}
	print "Other" if (exists $data{Other});
	foreach my $k3(@sam){
		print "\t$data{Other}{$k3}" if (exists $data{Other});
	}
	print "\n";
	print "Unclassified_KO" if (exists $data{Unclassified_KO});
	foreach my $k5(@sam){
		print "\t$data{Unclassified_KO}{$k5}" if (exists $data{Unclassified_KO});
	}
	print "\n";
}
