#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %opt;
GetOptions(\%opt,"row","rank:s","Z","x:f","U","D","BC");
@ARGV || die"Usage: perl cluster.pl <input.mat> > out.tree 
    -row        cluster in row, default in rank
    -x <flo>    to amplify eatch matrix data
    -U          to change zero value into -max_value
    -Z          use Z_value: (X-u)/Sd for cluster distance
    -D          input matrix just be distance matirx
    -BC         use Bray Cutris distance
    -rank       ranks for selected data\n";
#============================================================
my (@name,@data,@sel);
my (%pardis, %sel_name,%uniq_pair);
open IN,shift || die$!;
$opt{row} && <IN>;
while(<IN>){
    chomp;
    my @l = /\t/ ? split /\t/ : split;
    my $tax = shift @l;
    @sel && (@l = @l[@sel]);
    if($opt{row}){
        push @name,$tax;
        push @data,[@l];
    }else{
        if(@name){
            $opt{D} && !(defined $sel_name{$tax}) && next;
            for my $i(0 ..$#l){
                if($opt{D}){
                    for my $i(0..$#name){
                        my $tax2 = $name[$i];
                        $uniq_pair{"$tax $tax2"} && next;
                        $uniq_pair{"$tax $tax2"} = 1;
                        $pardis{$tax}{$tax2} = $l[$i];
                    }
                }else{
                    push @{$data[$i]},$l[$i];
                }
            }
        }else{
            if($opt{rank}){
                if($opt{rank} =~ /(\S+)\.\.(\S+)/){
                    my ($s,$e) = ($1,$2);
                    ($s < 0) && ($s += $#l);
                    ($e < 0) && ($e += $#l);
                    @sel = ($s .. $e);
                }else{
                    @sel = split/,/,$opt{rank};
                }
            }
            @sel && (@l = @l[@sel]);
            @name = @l;
            if($opt{D}){
                for(@l){$sel_name{$_} = 1;}
            }
        }
    }
}
close IN;
if(!$opt{D}){
    if($opt{x}){
        for my $i(0..$#data){
            for my $j(0..$#{$data[$i]}){
                $data[$i]->[$j] *= $opt{x};
            }
        }
    }
    if($opt{U}){
        for my $i(0 .. $#{$data[0]}){
            my $max = 0;
            for my $j(0 .. $#data){
                ($data[$j]->[$i] > $max) && ($max = $data[$j]->[$i]);
            }
            for my $j(0 .. $#data){
                $data[$j]->[$i] ||= -$max;
            }
        }
    }
    if($opt{Z}){
        for my $i(0 .. $#{$data[0]}){
            my ($avg,$sd,$num);
            for my $j(0 .. $#data){
                $num++;
                $avg += $data[$j]->[$i];
                $sd += $data[$j]->[$i] ** 2;
            }
            $avg /= $num;
            $sd = sqrt($sd/$num - $avg**2);
            for my $j(0 .. $#data){
                $data[$j]->[$i] = $sd ? ($data[$j]->[$i] - $avg) / $sd : 0;
            }
        }
    }
    for my $i(0..$#name){
        for my $j(0..$#name){
            $i > $j || last;
            $pardis{$name[$i]}{$name[$j]} = dis_cal($data[$i],$data[$j],$opt{BC});
        }
    }
}
my %sigdis;
for(@name){$sigdis{$_} = 0;}
my $tree;
while(keys %pardis){
    my ($mindis,$s1,$s2);
    for my $k1(keys %pardis){
        for my $k2(keys %{$pardis{$k1}}){
            if(!defined $mindis || $pardis{$k1}{$k2} < $mindis){
                ($mindis, $s1, $s2) = ($pardis{$k1}{$k2}, $k1, $k2);
            }
        }
    }
    my $avgdis = $pardis{$s1}{$s2} / 2;
    my $d1 = $avgdis - $sigdis{$s1};
    my $d2 = $avgdis - $sigdis{$s2};
    my $key = "($s1:$d1,$s2:$d2)";
    $sigdis{$key} = $avgdis;
    delete $sigdis{$s1};
    delete $sigdis{$s2};
    if(@name == 2){
        $tree = $key;
        last;
    }
    my @name2;
    for(@name){
        ($_ ne $s1 && $_ ne $s2) && (push @name2,$_);
    }
    for(@name2){
        my $cd1 = (defined $pardis{$s1} && defined ${$pardis{$s1}}{$_}) ? $pardis{$s1}{$_} : $pardis{$_}{$s1};
        my $cd2 = (defined $pardis{$s2} && defined ${$pardis{$s2}}{$_}) ? $pardis{$s2}{$_} : $pardis{$_}{$s2};
        $pardis{$key}{$_} = ($cd1 + $cd2) / 2;
    }
    push @name2,$key;
    for my $i(0..$#name){
        if($name[$i] eq $s1 || $name[$i] eq $s2){
            delete $pardis{$name[$i]};
            next;
        }
        for my $j(0..$#name){
            $i > $j || last;
            ($name[$j] eq $s1 || $name[$j] eq $s2) && (delete ${$pardis{$name[$i]}}{$name[$j]});
        }
    }
    @name = @name2;
}

print $tree,"\n";

sub dis_cal{
    my ($d1,$d2,$BC) = @_;
    my $dis = 0;
    my $sum = 0;
    for my $i(0..$#$d1){
        if($BC){
            $dis += (abs($d1->[$i]) < abs($d2->[$i])) ? abs($d1->[$i]) : abs($d2->[$i]);
            $sum += abs($d1->[$i]) + abs($d2->[$i]);
        }else{
            $dis += ($d1->[$i] - $d2->[$i]) ** 2;
        }
    }
    if($BC){
        1 - 2 * $dis / $sum;
    }else{
        $dis ** 0.5;
    }
}

