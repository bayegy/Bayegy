#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opt=(size=>10000);
GetOptions(\%opt,"size:s","out:s");
#==========================================================================================
@ARGV || die"Name: samples_draw.pl
Description: script for sample draw
Usage: perl samples_draw.pl <in.table> [-option]
    in.table <file> input frequency data table
    --size <str>    sampling size, default=20000
    --out <str>     outfile prefix, default STDOUT
Example:
  1. perl samples_draw.pl  otu_table.txt -size 10000 > otu_table.1w.txt
  2. perl samples_draw.pl otu_table.txt -size 10000,20000,50000,100000 --out otu_table.txt
Note:
  1. --size cat set mulit value, e.g: '10000,20000,50000,100000'.
  2. while --size set mulit, --out default=ARGV[0], outfile to be: out.v1 out.v2 out.v2\n";
#==========================================================================================
my $intab = shift;
(-s $intab) || die$!;
my @level = split /,/,$opt{size};
(@level > 1) && ($opt{out} ||= $intab);
$opt{out} && ($opt{out}=~/^(\S+)\//) && mkdir"$1";
my $head;
my @percent;
open IN,$intab || die$!;
while(<IN>){
    if($head){
        chomp;
        my @l = /\t/ ? split /\t/ : split;
        push @percent,[@l];
    }else{
        $head = $_;
    }
}
close IN;
for my $d(@level){
    my @out;
    my @simpling;
    my @scal;
    my @acc;
    my @tol;
    for my $i(0..$#percent){
        my @line = @{$percent[$i]};
        for my $j(0..$#line){
            if($j==0 || ($j==$#line && ($line[$j]=~/[^\d+\.e-]/ || $line[$j] eq '-'))){
                $out[$i]->[$j] = $line[$j];
            }else{
                my $temp = $line[$j]*$d;
                $out[$i]->[$j] = int($temp);
                $tol[$j] += $out[$i]->[$j];
                my $rest = $temp - $out[$i]->[$j];
                if($rest){
                    $acc[$j]+=$rest;
                    push @{$scal[$j]},[$i,$acc[$j]];
                }
            }
        }
    }
    simple_draw(\@out,\@acc,\@scal,\@tol,$d);#sub1
    if($opt{out}){
        open OUT,">$opt{out}.$d" || die$!;
        select OUT;
    }
    print $head;
    for (@out){
        print join("\t",@{$_}),"\n";
    }
    $opt{out} && close(OUT);
}
#==========================================================================================
#sub1
sub simple_draw{
    my ($out,$acc,$scal,$tol,$d) = @_;
    ($acc && @$acc) || return(0);
    for my $j(0..$#$acc){
        $acc->[$j] || next;
#        my $num = int($acc->[$j]+0.5);
        my $num = $d - $tol->[$j];       
        for (1..$num){
            my $i = sel_sample($scal->[$j],$acc->[$j]);#sub1.1
            $out->[$i]->[$j]++;
        }
    }
}
#sub1.1
sub sel_sample{
    my ($scal,$i) = @_;
    $i = rand($i);
    for (@$scal){
        ($i <= $_->[1]) && return($_->[0]);
    }
}
