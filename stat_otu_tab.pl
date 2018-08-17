#!/usr/bin/perl -w

#Annotated by Cheng Guo, 2018.05.03. The script is confirmed correct that all relative/absolute abundance is as expected.

use strict;
use FindBin qw($Bin);
use Getopt::Long;
my ($prefix,$outsel,$abs,$sufix,$nomat,$even,$unif,$spestat);

GetOptions("prefix:s"=>\$prefix,"outsel:s"=>\$outsel,"abs"=>\$abs,"sufix:s"=>\$sufix,"nomat"=>\$nomat,
    "even:s"=>\$even,"unif:s"=>\$unif,"spestat:s"=>\$spestat);
@ARGV || die"Name: stat_otu_tab.pl

Usage: perl stat_otu_tab.pl <otu_table.txt> [prefix]
    --prefix <str>      set outfile prefix, default to be ARGV[0] or ARGV[1]
    --abs               output stat at absolute aubndant, default at relative
    --outsel <str>      to select output level, e.g: k,p,o,f,g,s, default output all
    --even <file>       output uniformization OTU matrix, default not output
    --unif <str>        set uniformization Tag number, min/max/med/avg means the min/max/med/avg Tag number, default=med
    --spestat <file>    output species Number stat at {k,p,c,o,f,g,s} level
    outfile contain: 1 prefix.relative.mat, 2 prefix.{k,p,c,o,f,g,s}.relative.mat

Note:
    1. -even not used while -abs or -nomat.
Update:
    1. Others' relative abundance=0 if <0, 2014-10-16,chen,line 133\n";
#====================================================================================

#A typical command used: $ perl ~/pipelines/Github/Queena/stat_otu_tab.pl -unif min feature-table.taxonomy.txt -prefix Relative/otu_table --even Relative/otu_table.even.txt -spestat Relative/classified_stat_relative.xls

my (@spe_name,@tag_num);
my (%matrix,%outsel);
my ($otu_tab) = @ARGV;

$prefix ||= ($ARGV[1] || $otu_tab);

#make directory for Relative and put the otu_table as the prefix for any file generated.
if($prefix =~ /^(\S+)\/(\S+)/){
    (-d $1) || mkdir"$1";
    (-w $1) || ($prefix = $2);
}

#output level of the taxonomy
if($outsel){
    for(split/,/,$outsel){
        $outsel{$_} = 1;
    }
}


#generate the prefix.relative.mat table. 
Tag_stat($otu_tab,\@tag_num);
if(!$nomat){
    open OUT,">$prefix.relative.mat" || die$!;
}

($nomat || $abs) && ($even = 0);
$unif ||= 'med';
my %samp_spe_num;

#generate the %matrix with the max() function correctly but I am not totally follow the script.
open IN,$otu_tab || die$!;
while(<IN>){
    chomp;
    my @l = /\t/ ? split /\t/ : split;
    #print "BBBBBBBBBBBBBBBB", $_, "\n";
    if(!@spe_name){
        /^#OTU/ || next;
        $nomat || (print OUT $_,"\n");
        @spe_name = @l[1..$#l];
        ($spe_name[-1] =~ m/^Tax/i) && (pop @spe_name);
        #print "CCCCCCCCCCC", $spe_name[-1], "\n";
    }elsif(!@tag_num){
        @tag_num = @l[1..$#l];
    }elsif(!/^#/){
        my @ll;
        $abs && (@ll = @l);
        #print "DDDDDDDDDDD", "\t", $abs, "\n";
        if(!($nomat && $abs)){
            for my $i(0..$#tag_num){
                #print "AAAAAAAAAAAAAAAAAAAAAAAAAA", "\t", $i, "\t", $#tag_num, "\t", $l[$i+1], "\t", $tag_num[$i], "\n";
                $l[$i+1] /= $tag_num[$i];
                #print "AAAAAAAAAAAAAAAAAAAAAAAAAA", "\t", $i, "\t", $#tag_num, "\t", $l[$i+1], "\t", $tag_num[$i], "\t", $l[$i+1], "\n";
            }
        }
        $nomat || (print OUT join("\t",@l),"\n");
        $abs && (@l = @ll);
        my $full_tax;
        for (split/;/,$l[-1]){
            m/(\w)__(.+)/ || next;
            $outsel && !$outsel{$1} && next;
            my ($level,$tax) = ($1,$2);
            $full_tax .= "$level\__$tax;";
            for my $i(0..$#tag_num){
                $matrix{$level}{$full_tax}->[$i] += $l[$i+1];
                $spestat && ($samp_spe_num{$level}->[$i] += $l[$i+1]);
            }
            $matrix{$level}{$full_tax}->[$#tag_num+1] = max(@{$matrix{$level}{$full_tax}});
            #print "EEEEEEEEEEEEE", $matrix{$level}{$full_tax}->[$#tag_num+1], "\n";
        }
    }
}
close IN;
$nomat || close(OUT);

#output species Number stat at {k,p,c,o,f,g,s} level
if($spestat){
    open SPT,">$spestat" || die$!;
    my @class = qw(Kingdom Phylum Class Order Family Genus Species);
    my @class_short = qw(k p c o f g s);
    my @class_head;
    for my $i(0 .. $#class_short){
        $samp_spe_num{$class_short[$i]} && (push @class_head,$class[$i]);
    }
    print SPT join("\t","Sample_Name",@class_head),"\n";
    for my $i(0 .. $#spe_name){
        my @class_out;
        push @class_out,$spe_name[$i];
        for my $j(0..$#class_head){
            push @class_out,($samp_spe_num{$class_short[$j]}->[$i] || 0);
        }
        print SPT join("\t",@class_out),"\n";
    }
    close SPT;
}

#Here is calling the sample_draw.pl script to draw samples with different sampling size
if($even){
    my @tnum = sort {$a<=>$b} @tag_num;
    if($unif=~/min/){
        $unif = $tnum[0];
    }elsif($unif=~/max/){
        $unif = $tnum[-1];
    }elsif($unif=~/med/){
        $unif = int(($tnum[int($#tnum/2)]+$tnum[int($#tnum/2+0.5)])/2);
    }elsif($unif=~/avg/){
        $unif = 0;
        for(@tnum){$unif += $_;}
        $unif = int($unif / @tnum);
    }
    system"perl $Bin/samples_draw.pl $prefix.relative.mat -size $unif > $even";
    #print "AAAAAAAAAAAAAAAAAAAAAAAAAAa",$unif, "\n";
}

#check $abs is difined or not
$sufix ||= $abs ? "absolute" : "relative";


for my $level(keys %matrix){
    open OUT,">$prefix.$level.$sufix.mat" || die$!;
    print OUT join("\t","Taxonomy",@spe_name,"Tax_detail"),"\n";
    my @tol_tax;
    for my $full_tax(sort {$matrix{$level}{$b}->[-1]<=>$matrix{$level}{$a}->[-1]} keys %{$matrix{$level}}){
        my @out = @{$matrix{$level}{$full_tax}};
        pop @out;
        my $tax = (split/\w\__/,$full_tax)[-1];
        $tax =~ s/;//;
        print OUT join("\t",$tax,@out,$full_tax),"\n";
        for my $i(0..$#out){
            $tol_tax[$i] += $out[$i];
        }
    }
    if($abs){
        for my $i(0..$#tol_tax){
            $tol_tax[$i] = $tag_num[$i] - $tol_tax[$i];
        }
    }else{
        for(@tol_tax){$_ = 1 - $_;$_=0 if($_<0);} #add for Others' relative abundance =0 if <0, 2014-10-16,chen
    }
    print OUT join("\t","Others",@tol_tax,"Others\n");
    close OUT;
}

#=======================================================================================================
sub Tag_stat{
    my ($otu_tab,$tag_num) = @_;
#    ($otu_tab && -s $otu_tab && `awk 'NR==3' $otu_tab`=~/^#/) && return(0);
    open OTU,$otu_tab || die$!;
    while(<OTU>){
        /^#/ && next;
        my @l = /\t/ ? split /\t/ : split;
        for my $i(1..$#l-1){
            $tag_num->[$i-1] += $l[$i];
        }
    }
    close OTU;
}

sub max{
    my $max = $_[0];
    for (@_){($_ > $max) && ($max = $_);}
    $max;
}

