#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/";
use PATHWAY;
(-s "$Bin/Pathway_cfg.txt") || die"error: can't find config at $Bin/bin, $!\n";
BEGIN {
    my( $svg_lib)= get_pathway("$Bin/Pathway_cfg.txt",[qw(SVG_Lib)]);
    unshift @INC, $svg_lib
};
use SVG;
my %opt = (flank_x=>40,flank_y=>20,height=>300,flank_y=>20,lw=>1.5,fsize=>15,type=>0,barlen=>400,
        color1=>'black',color2=>'black',color3=>'blue',#scal_title=>"Divergence, substitutions/site",
        range=>1,tran=>0,rate2=>0.9);
GetOptions(
    \%opt,"flank_x:f","flank_y:i","width:i","height:i","fsize:s","lw:f","end","type:i",
    "color1:s","color2:s","color3:s","turn:s","root","scal_title:s","top","one","color:s",
    "barlen:i","btitle:s","bcolor:s","trantab","flank_b:i","bline","bun:s","bh:f","group:s",
    "range:f","tran:f","noscal","cxy:s","cr:s","opacity:f","inner_fill:i","dot","spcolor","lgranks:s",
    "noname","species:s","spematch:s","symbols:s","branks:s","bsplit","rate:f","nohead","rate2:f"
);
$opt{width} ||= ($opt{type}==4) ? 600 : 400;
#=======================================================================================
@ARGV || die"Name: draw_tree.pl
Description: scription to draw evolution tree
Author: Wenbin Liu
Version: 2.0,  Date: 2014-03-06,  Modified: add function to draw bar diagram 
Usage: perl draw_tree.pl <in.tree.nhx> [bar.stat] > out.tree.svg
    --trantab           to transposition input bar.stat table
    --flank_x <flo>     figure x flank edeg distance, default=50
    --flank_y <flo>     figure y flank edge distance, default=20
    --width <num>       figure width, default=400
    --height <num>      figure height, default=300
    --barlen <num>      bar diagram width, default=400
    --flank_b <num>     bar diagram flank edge distance, default=barlen/4
    --fsize <num>       font size, default=15
    --lw <flo>          line width, default=1.5
    --end               species name alignment to the end
    --type <num>        line type: 0-vertical ,1-oval, 2-slant, 3-smooth, 4-cycle default=0
    --color1 <str>      line color, default=black
    --color2 <str>      species name color, default=black
    --color3 <str>      NHX text color, default=blue
    --color <file>      set species color, form: spename color, default not set
    --group <file>      set group and mark line color according to group
    --bcolor <str>      set bar colors, default=''
    --root              to add the false root branch
    --scal_title <str>  default='Divergence, substitutions/site'
    --btitle <str>      bar diagram title, default not set
    --bline             to add group chain line at bar figure
    --bun <str>         to set bar figure scale unit,number e.g:0.25,4 default auto
    --one               just one scal
    --top               put scal_title to top
    --turn <str>        set species at specify turn
    --species <str>     only specified species for draw
    --symbols <file>    to set symbols\n\n";
#=======================================================================================
my $rgb_txt = "$Bin/rgb_colors.txt";
my ($nhx,$barstat) = @ARGV;
(-s $nhx) && ($nhx = `less $nhx`);
my (%spe_num,%node_yh,%node_xw,%edge,%bardata);
my (@spe,@tree,@branch,@bsym,@end_branch);
if($opt{trun}){
    (-s $opt{turn}) && ($opt{turn} = `less $opt{turn}`);
    $opt{turn} =~ s/;|\s|\&\&NHX:|\[.+?\]//g;
    @spe = ($opt{turn}=~/:/) ? ($opt{turn} =~ /([^(),]+?):/g) : ($opt{turn} =~ /([^(),]+?)/);
}
my @barlens;
my @branks = $opt{branks} ? split/,/,$opt{branks} : ();
my $lgranks = defined $opt{lgranks} ? [(split/,/,$opt{lgranks})] : 0;
my $blen = read_bardata($barstat,\@bsym,\%bardata,$opt{trantab},$opt{bsplit},\@barlens,\@branks,$opt{nohead},$lgranks);
my @bcolor = $opt{bcolor} ? split/,/,$opt{bcolor} :
#    qw(crimson blue lightseagreen orange mediumpurple palegreen lightcoral dodgerblue lawngreen red olive 
#    yellow fuchsia salmon mediumslateblue darkviolet purple sienna tan chocolate skyblue turquoise cadetblue green);
qw(pink orange green cyan blue purple cornflowerblue red darkturquoise sienna bisque blueViolet orangered olive lightseagreen crimson salmon yellow fuchsia mediumslateblue darkviolet tan chocolate skyblue turquoise cadetblue navajowhite slategrey cornflowerblue royalblue dodgerblue deepskyblue mediumaquamarine lawngreen yellowgreen gold saddlebrown indianred deeppink darkred peachpuff);
my %colorh;
my (@colors, @gline);
if($opt{color} && -s $opt{color}){
   %colorh = split/\s+/,`less $opt{color}`;
}elsif($opt{group} && -s $opt{group}){
    my %uniq_group;
    for(`less $opt{group}`){
        chomp;
        my @l = split /\s+/,$_,2;
        if(!$uniq_group{$l[1]}){
            push @gline,$l[1];
            $uniq_group{$l[1]} = $bcolor[$#gline];
        }
        $colorh{$l[0]} = $uniq_group{$l[1]};
    }
}
my $rate = $opt{rate} || 0.8;
my @symbols;
my %sym_colors;
if($opt{symbols}){
    my $s = 0;
    if(-s $opt{symbols}){
        for (`less $opt{symbols}`){
            chomp;
            my @l = /\t/ ? split /\t/ : split;
            push @symbols,$l[0];
            (!$l[1] || $l[1]=~/^[\.\d]+$/) && ($l[1] = $bcolor[$s]);
            $sym_colors{$l[0]} = $l[1];
            $s++;
        }
    }else{
        @symbols = split/,/,$opt{symbols};
        for(@symbols){
            $sym_colors{$_} = $bcolor[$s];
            $s++;
        }
    }
}
my %species_sel;
if($opt{species}){
    if(-s $opt{species}){
#        %species_sel = split/\s+/,`awk '{print \$1,1}' $opt{species}`;
        for(`less $opt{species}`){
            chomp;
            my @l = /\t/ ? split/\t/ : split;
            $l[1] && $sym_colors{$l[1]} && !$colorh{$l[0]} && ($colorh{$l[0]}=$sym_colors{$l[1]});
            $species_sel{$l[0]} = 1;
        }
    }else{
        for(split/,/,$opt{species}){$species_sel{$_} = 1;}
    }
}elsif($opt{spematch}){
    $opt{species} = 1;
}
if(-s $rgb_txt && %colorh){
    my %rgbh = split/\s+/,`less $rgb_txt`;
    for (values %colorh){$rgbh{$_} && ($_=$rgbh{$_});}
}
my ($root_len,$root_B,$no_blen) = get_nhx($nhx,\@spe,\@tree,\@branch,\%spe_num,\%node_yh,\%edge,\%colorh,
    \@colors,\@end_branch,$opt{spematch},$opt{species} ? \%species_sel : 0);#sub1
$no_blen && ($opt{scal_title} = "");
if($opt{bh}){
    $opt{height} = $opt{bh} * @spe;
}elsif($opt{height}/@spe < 8){
    $opt{height} = 8*@spe;
}
my $max_branch = 0;
for(@branch){($_ > $max_branch) && ($max_branch = $_);}
my ($x_unit, $y_unit, $cr1, $cr2);
my $fsize = $opt{fsize};
if($opt{type}==4){
    my $max_fsize = $opt{width}*6/$#branch;
    ($max_fsize < $opt{fsize}) && ($opt{fsize} = $max_fsize);
    my $branch_num = ($opt{range} == 1) ? @branch : $#branch;
    if($opt{cr}){
        ($cr1,$cr2) = split/,/,$opt{cr};
    }else{
        $cr1 = $opt{width} / 2;
    }
    $cr2 ||= $cr1/8;
    ($x_unit, $y_unit) = (($cr1-$cr2) / ($max_branch+$root_len), $opt{range} / $branch_num);
}else{
    ($#branch*$opt{fsize} > $opt{height}) && ($opt{fsize} = $opt{height} / $#branch);
    ($x_unit, $y_unit) = ($opt{width} / ($max_branch+$root_len), $opt{height} / ($#branch+1));
}
my $max_width = 0;
for my $i(0..$#branch){
    my $temp_len = ($opt{end} ? ($opt{type}==4 ? $cr1 : $opt{width}) : ($opt{type}==4 ? $cr2 : 0) + $branch[$i]*$x_unit) +
       ($opt{noname} ? $opt{fsize} : 0.6*$opt{fsize}*(length($spe[$i])+3));
    ($temp_len > $max_width) && ($max_width = $temp_len);
}
my $bar_star = ($opt{type}==4) ? $max_width : $max_width+$opt{flank_x};
$max_width -= $opt{width};
($max_width < $opt{flank_x}) && ($max_width = $opt{flank_x});
$opt{scal_title} && ($opt{scal_title} !~ /\S/) && (delete $opt{scal_title});
my $max_height = $opt{scal_title} ? 3.8*$opt{fsize} : 2.8*$opt{fsize};
($max_height < $opt{flank_y}) && ($max_height = $opt{flank_y});
my $max_len = $opt{flank_x} + $opt{width};
my $pwidth = $opt{width} + $opt{flank_x} + $max_width;
my $pheight = $opt{height} + $opt{flank_y} + $max_height;
$opt{flank_b} ||= $opt{barlen} / 4;
my ($cx,$cy);
if($opt{type} == 4){
    $blen && ($pwidth += $opt{barlen});
    $pwidth = 2*($pwidth - $opt{width} + $cr1);
    if($opt{cxy}){
        ($cx,$cy) = split/,/,$opt{cxy};
    }else{
        $cx = ($cy = $pwidth/2);
    }
    $pheight = $pwidth;
}elsif($blen){
    $pwidth += $opt{barlen}+$opt{flank_b};
}
my $svg = SVG->new(width=> $pwidth,height=> $pheight);

if($opt{type} == 4){
    inner_fill($svg,\@spe,\@colors,$cx,$cy,$cr1,$cr2,$opt{tran},$x_unit,$y_unit,$opt{opacity},\@end_branch,\@branch,$opt{inner_fill},
        $bar_star,$opt{barlen});
    draw_cycle_tree($svg,$cx,$cy,$cr2,$cr1,$opt{root},$root_len,$root_B,$x_unit,$y_unit,\@tree,\@spe,\%node_xw,\%node_yh,$opt{fsize},
        $opt{lw},$opt{end},$opt{color1},$opt{color2},$opt{color3},\@colors,\%edge,$opt{tran},$opt{opacity},$opt{dot},$opt{noname});
    if($blen){
        draw_cycl_bar($svg,\@spe,\%bardata,$cx,$cy,$bar_star,$opt{barlen},$opt{tran},$rate,$opt{rate2},$y_unit,$blen,\@bcolor,\@colors,
            $opt{bline},$opt{spcolor},$opt{bsplit},\@barlens);
        @bsym || (@bsym = @symbols);
        draw_sym($svg,\@bsym,$fsize,$fsize,$pwidth/3,$fsize,\@bcolor,1);
    }elsif(!($no_blen || $opt{noscal})){
        draw_scal($svg,$max_branch,$x_unit,$opt{width},2*$opt{fsize},2*$opt{fsize},$opt{fsize},$opt{lw},$opt{scal_title},1);
    }
}else{
    $opt{top} && ($opt{flank_y} = $max_height);
    draw_tree($svg,$opt{flank_x},$opt{flank_y},$opt{root},$root_len,$root_B,$x_unit,$y_unit,\@tree,\@spe,\%node_xw,
        \%node_yh,$opt{fsize},$max_len,$opt{lw},$opt{end},$opt{type},$opt{color1},$opt{color2},$opt{color3},\@colors,\%edge,
        $opt{inner_fill},$opt{opacity},$opt{noname});
    my ($scal_x,$scal_y) = ($opt{flank_x}, $opt{top} ? 0 : $opt{flank_y}+$opt{height}+$opt{fsize}/3);
    $no_blen ||
    draw_scal($svg,$max_branch,$x_unit,$opt{width},$scal_x,$scal_y,$opt{fsize},$opt{lw},$opt{scal_title},$opt{one});#sub3
    if($blen){
        draw_bar($svg,\@spe,\%bardata,$y_unit,$bar_star,$opt{flank_y},$opt{barlen},$rate,$opt{fsize},$blen,
        \@bcolor,$opt{bline},$opt{bun},$opt{btitle});
        draw_sym($svg,\@bsym,$bar_star+$opt{barlen},$opt{flank_y}+(1-$rate)*$y_unit/2,$opt{flank_b},$opt{fsize},\@bcolor);
    }
    @gline && draw_sym($svg,\@gline,$opt{fsize},$opt{flank_y},$opt{flank_x},$opt{fsize},\@bcolor,2);
}
print $svg->xmlify;

#===================================================================
sub draw_cycl_bar{
    my ($svg,$spe,$bardata,$cx,$cy,$barr,$bwidth,$tran,$rate,$rate2,$y_unit,$blen,$bcolor,$colors,$chain,$spcolor,$barsplit,$blens) = @_;
    my $x_unit = $barsplit ? $bwidth/$barsplit : $bwidth / $blen;
    my @x_units;
    if($barsplit){
        for (@$blens){ push @x_units,$rate2*$x_unit/$_;}
    }
    my $bh = $y_unit * $rate;
    my $by = $tran-$bh/2;
    my (@per_line,@first_line);
    my @line_color = ('black');
    for my $k (0..$#$spe){
#        $bardata->{$spe->[$k]} || die"error: can't find tree-species $k at bardata\n";
        if(!defined $bardata->{$spe->[$k]}){
            $by += $y_unit;
            next;
        }
        my @bd = @{$bardata->{$spe->[$k]}};
        my $br = $barr;
        my @temp_line;
        my $sp_color = $colors->[$k] || 'black';
        ($sp_color =~ /^\d+,\d+,\d+$/) && ($sp_color = "rgb($sp_color)");
        for my $i(0..$#bd){
            my $bL = $bd[$i] * ($x_units[$i] || $x_unit);
            my $bar_col = $spcolor ? $sp_color : $bcolor->[$i];
            my @xy = annulus_path($svg,'d',[$cx,$cy,$br,$br+$bL,$by,$by+$bh,0],'stroke','none','fill',$bar_col);
            $br += $barsplit ? $x_unit : $bL;
            $chain || next;
            $k || ($line_color[$i+1] = $sp_color ? 'black' : $bar_col);
            @temp_line ? (push @temp_line,$xy[1]) : (push @temp_line,@xy);
        }
        $by += $y_unit;
        @temp_line || next;
        if(@per_line){
            for my $i(0..$#per_line){
                $svg->line('x1',$per_line[$i]->[2],'y1',$per_line[$i]->[3],'x2',$temp_line[$i]->[0],'y2',
                    $temp_line[$i]->[1],'stroke',$line_color[$i],'stroke-width',$chain);
            }
        }else{
            @first_line = @temp_line;
        }
        @per_line = @temp_line;
    }
    if(@per_line){
        for my $i(0..$#per_line){
            $svg->line('x1',$per_line[$i]->[2],'y1',$per_line[$i]->[3],'x2',$first_line[$i]->[0],'y2',
                $first_line[$i]->[1],'stroke',$line_color[$i],'stroke-width',$chain);
        }
    }
}

sub inner_fill{
    my ($svg,$spe,$colors,$cx,$cy,$r2,$r1,$tran,$x_unit,$y_unit,$opacity,$end_branch,$branch,$inner_fill,$bar_star,$barlen,$add) = @_;
    ($inner_fill && $colors && @$colors) || return(0);
    $add ||= 0;
    $opacity ||= 0.3;
    my $a0 = $tran - $y_unit / 2;
    my $r10 = $r1;
    my ($a00,$r00,$col0) = ($a0);
    for my $i(0..$#$spe){
        my $col = $colors->[$i];
        ($col =~ /^\d+,\d+,\d+$/) && ($col = "rgb($col)");
        ($inner_fill >= 2) && ($r1 = $r10 + ($branch->[$i] - $end_branch->[$i])*$x_unit);
        if($inner_fill > 2){
            $col0 ||= $col;
            $r00 ||= $r1;
            if($col0 ne $col){
                annulus_path($svg,'d',[$cx,$cy,$r00,$r2,$a00,$a0,$add],'fill-opacity', $opacity,
                'stroke','none','fill',$col0);
                ($inner_fill > 3) &&  annulus_path($svg,'d',[$cx,$cy,$bar_star,$bar_star+$barlen,$a00,$a0,$add],
                    'fill-opacity', $opacity,'stroke','none','fill',$col0);
                ($r00,$a00,$col0) = ($r1,$a0,$col);
            }elsif($r1 < $r00){
                $r00 = $r1;
            }
        }else{
            annulus_path($svg,'d',[$cx,$cy,$r1,$r2,$a0,$a0+$y_unit,$add],'fill-opacity', $opacity,
            'stroke','none','fill',$col);
        }
        $a0 += $y_unit;
    }
    ($inner_fill > 2) && annulus_path($svg,'d',[$cx,$cy,$r00,$r2,$a00,$a0,$add],'fill-opacity', $opacity,
        'stroke','none','fill',$col0);
    ($inner_fill > 3) &&  annulus_path($svg,'d',[$cx,$cy,$bar_star,$bar_star+$barlen,$a00,$a0,$add],
        'fill-opacity', $opacity,'stroke','none','fill',$col0);
}
sub read_bardata{
    my ($barstat,$bsym,$bardata,$tran,$bsplit,$barlens,$branks,$nohead,$lgranks) = @_;
    ($barstat && -s $barstat) || return(0);
    my $blen = 0;
    my @head;
    for(`less $barstat`){
        chomp;
        my @l = /\t/ ? split /\t/ : split;
        $branks && @$branks && (@l = @l[@$branks]);
        my $key = shift @l;
        if($tran){
            if(@head || $nohead){
                for my $i(0..$#head){
                    push @{$bardata->{$head[$i]}},$l[$i];
                }
                push @$bsym,$key;
            }else{
                @head = @l;
            }
        }else{
            (@$bsym || $nohead) ? ($bardata->{$key} = [@l]) : (@$bsym = @l);
        }
    }
    my @minbarlen;
    for my $v(values %{$bardata}){
        my $tem_len = sum(@$v);
        ($tem_len > $blen) && ($blen = $tem_len);
        if($opt{bsplit}){
            if(@$barlens){
                for my $i(0..$#$v){
                    ($v->[$i] > $barlens->[$i]) && ($barlens->[$i] = $v->[$i]);
                }
                if($lgranks && @$lgranks){
                    for my $i(@$lgranks){
                        ($v->[$i] < $minbarlen[$i]) && ($minbarlen[$i] = $v->[$i]);
                    }
                }
            }else{
                @minbarlen = (@$barlens = @$v);
                $opt{bsplit} = @$v;
            }
        }
    }
    if($lgranks && @$lgranks){
        for my $i(@$lgranks){
            my $min = log($minbarlen[$i]);
            my $max = log($barlens->[$i]) - $min;
            $barlens->[$i] = 1;
            for my $v(values %{$bardata}){
                $v->[$i] = (log($v->[$i])-$min)/$max;
            }
        }
    }       
    $blen;
}
sub draw_sym{
    my ($svg,$bsym,$flank_x,$flank_y,$bwidth,$fsize,$color,$cycle) = @_;
    my $nlen = 0;
    $cycle ||= 0;
    for (@$bsym){
        my $temp_len = length;
        ($temp_len > $nlen) && ($nlen = $temp_len);
    }
    my $may_size = $bwidth / (0.6*$nlen+3);
    ($may_size < $fsize) && ($fsize = $may_size);
    my $sx = $flank_x + $fsize;
    my $sy = $flank_y;
    for my $i(0..$#$bsym){
        if($cycle==2){
            $svg->line('x1',$sx-$fsize/2,'y1',$sy+$fsize/2,'x2',$sx+$fsize,'y2',$sy+$fsize/2,'stroke',$color->[$i]);
        }elsif($cycle==1){
            $svg->circle('cx',$sx+$fsize/2,'cy',$sy+$fsize/2,'r',$fsize/2,'stroke','none','fill',$color->[$i]);
        }else{
            $svg->rect('x',$sx,'y',$sy,'width',$fsize,'height',$fsize,'stroke','black','fill',$color->[$i]);
        }
        $svg->text('x',$sx+1.5*$fsize,'y',$sy+0.85*$fsize,'stroke','none','fill','black','-cdata',$bsym->[$i],
                'font-size', $fsize, 'text-anchor','start','font-family','Arial');
        $sy += 1.5*$fsize;
    }
}
sub draw_bar{
    my ($svg,$spe,$bardata,$y_unit,$flank_x,$flank_y,$barlen,$rate,$fsize,$blen,$bcolor,$chain,$bun,$btitle) = @_;
    my ($bs_min,$bs_num) = $bun ? split/,/,$bun : axis_split($blen);
    my $x_unit = $barlen/$bs_min/$bs_num;
    my $by = $flank_y + (1 - $rate)*$y_unit/2;
    my $bh = $y_unit * $rate;
    my @ln;
    ($rate < 1) || ($chain = 0);
    for my $k(@$spe){
        my $bx = $flank_x;
        $bardata->{$k} || die"error: can't find tree-species $k at bardata\n";
        my @bd = @{$bardata->{$k}};
        for my $i(0..$#bd){
            my $rect_w = $bd[$i] * $x_unit;
            simpnum($bx,$by,$rect_w,$bh,2);
            if($rect_w > 0){
                $svg->rect('x',$bx,'y',$by,'width',$rect_w,'height',$bh,'stroke','black','fill',$bcolor->[$i]);#,'stroke-width',1);
                $bx += $rect_w;
            }
            $chain || next;
            $ln[$i] && simpnum($ln[$i]->[0],$ln[$i]->[1],2);
            $ln[$i] && 
                $svg->line('x1',$ln[$i]->[0],'y1',$ln[$i]->[1],'x2',$bx,'y2',$by,'stroke','black','stroke-width',$chain);
            $ln[$i] = [$bx,$by+$bh];
        }
        $by += $y_unit;
    }
    simpnum($flank_x,$by,$barlen,$fsize,2);
    $svg->line('x1',$flank_x,'y1',$by,'x2',$flank_x+$barlen,'y2',$by,'stroke','black','stroke-width',1);
    my $bx = $flank_x;
    for my $num(0..$bs_num){
        simpnum($bx,2);
        $svg->line('x1',$bx,'y1',$by,'x2',$bx,'y2',$by+$fsize/2,'stroke','black','stroke-width',1);
        $svg->text('x',$bx,'y',$by+1.5*$fsize,'stroke','none','fill','black','-cdata',$num*$bs_min,
                'font-size', $fsize, 'text-anchor','middle','font-family','Arial');
        $bx += $barlen/$bs_num;
    }
    if($btitle){
        simpnum($barlen,2);
        $svg->text('x',$flank_x+$barlen/2,'y',$by+2.5*$fsize,'stroke','none','fill','black','-cdata',$btitle,
            'font-size', $fsize,'text-anchor','middle','font-family','Arial');
    }
}
sub simpnum{
    my $num = pop;
    my $lim = 10**$num;
    for(@_){
        ($_ > 0.5/$lim) && ($_ = int($lim*$_+0.5)/$lim);
    }
}
sub sum{
    my $sum = 0;
    for(@_){$sum += $_;}
    $sum;
}
#sub1
#===========
sub get_nhx{
#===========
    my ($nhx,$spe,$tree,$branch,$spe_num,$node_yh,$edge,$colorh,$colors,$end_branch,$spematch,$spe_sel) = @_;
    $nhx =~ s/\'[\.\d]+?:/:/g;
    $nhx =~ s/;|\s|\&\&NHX:|'//g;
    my $i = "0000001";
    my $j = 0;
    my $no_blen = ($nhx =~ /:\d/) ? 0 : 1;
    if(!@$spe){
        if($no_blen){
            for($nhx =~ /([^(),]+?)/g){
                /^\[|^:|\]$/ && next;
                s/\:.+//g;
                push @$spe,$_;
            }
        }else{
            $nhx =~ s/\)[\.\d]+/\)/g;
            @$spe = ($nhx =~ /([^(),]+?):/g);
        }
    }
    my %del_spe;
    if($spe_sel && ($spematch || %$spe_sel)){
        my @new_spe;
        for (@$spe){
            if($spematch){
                /$spematch/ ? ($spe_sel->{$_}=1, push @new_spe,$_) : ($del_spe{$_} = 1);
            }else{
                $spe_sel->{$_} ? (push @new_spe,$_) : ($del_spe{$_} = 1);
            }
        }
        @$spe = @new_spe;
    }
    for(@$spe){
        $spe_num->{$_} = $i;
        $node_yh->{$i} = $j;
        push @{$edge->{$i}},$j;
        $colorh->{$_} && ($colors->[$j] = $colorh->{$_});
        $i++;$j++;
    }
    $i = 0;
    my $di = 0;
    my %single_node;
    my %sel_node = $spe_sel ? %$spe_sel : ();
    my %node_change;
    while($nhx && $nhx =~ /\([^()]+\)/){
        my @node = ($nhx =~ /(\([^()]+\))/g);
        my @sub_tree;
        my $dj = 0;
        for my $j(0..$#node){
            my $old_node = $node[$j];
            $old_node =~ s/([()\[\]\.])/\\$1/g;
            my $new_node0 = "$di->$dj";
            my $new_node = "$i->$j";
            $node_change{$new_node} = $new_node0;
            $node[$j] =~ s/^\(|\)$//g;
            my ($node_y,$node_n);
            for (split/,/,$node[$j]){
                /\S/ || next;
                my @T = split(/:/,$_,2);
                if($del_spe{$T[0]}){
                    delete $del_spe{$T[0]};
                    next;
                }
                if($spe_sel && (defined $spe_num->{$T[0]} || $T[0]!~/->/)){
                    $spe_sel->{$T[0]} ? delete $sel_node{$T[0]} : next;
                }
                $T[1] ||= $no_blen;
                if($T[1]=~/(.+?)\[(.+)\]/){
                    @T[1,2] = ($1,$2);
                }else{
                    $T[2] = 0;
                }
                if($spe_num->{$T[0]}){
                    $T[0] = $spe_num->{$T[0]};
                    $end_branch->[$node_yh->{$T[0]}] = $T[1];
                }
                (defined $node_yh->{$T[0]}) || die"$T[0]";
                $node_y += ($T[3] = $node_yh->{$T[0]});
                $T[4] = $new_node;
                for my $k(@{$edge->{$T[0]}}){
                    $branch->[$k] += $T[1]; # species $k total branch length
                    push @{$edge->{$new_node}},$k; # node_name -> [species contain]
                }
                if($single_node{$T[0]}){
                    $T[1] += $single_node{$T[0]}->[1];
                    ($single_node{$T[0]}->[0] !~ /->/) && 
                        ($end_branch->[$node_yh->{$single_node{$T[0]}->[0]}] = $T[1]);
                    $single_node{$T[0]}->[2] && ($T[2] = $single_node{$T[0]}->[2]);
                    $T[0] = $single_node{$T[0]}->[0];
                }
                push @{$sub_tree[$dj]},[@T]; # node_nmae branch_length &&NHX y-corrdinate
                $node_n++;
            }
            if($node_n){
                $nhx =~ s/$old_node([\.\d]+):([\.\d]+)/$new_node:$2\[B=$1\]/ ||
                $nhx =~ s/$old_node/$new_node/;
                $node_yh->{$new_node} = $node_y / $node_n; # node_name -> y-axis coordinate
                if($node_n == 1){
                    my $end_sub = pop @sub_tree;
                    $single_node{$new_node} = [@{$end_sub->[0]}];
                }else{
                    $dj++;
                }
            }else{
                $nhx =~ s/$old_node([^\),]+),// || $nhx =~ s/,$old_node([^\)\,)]*)// ||
                    $nhx =~ s/\($old_node\)// || $nhx =~ s/$old_node//;
            }
        }
        if(@sub_tree){
            push @$tree,[@sub_tree];
            $di++;
        }
        $i++;
        ($spe_sel && !%sel_node && $nhx!~/->\d+.+\d+->/) && ($nhx = "");
    }
    if($spe_sel && %$spe_sel){
        my (%new_edge,%new_node_yh);
        my $top_node;
        for my $i (0 .. $#$tree){
            for my $lev (@{$tree->[$i]}){
                for my $p (@$lev){
                    my $node_name = $p->[0];
                    my $new_name = $node_change{$node_name} || $node_name;
                    $p->[0] = $new_name;
                    $top_node = $p->[4];
                    $p->[4] = $node_change{$p->[4]};
                    $new_edge{$new_name} = [@{$edge->{$node_name}}];
                    $new_node_yh{$new_name} = $node_yh->{$node_name};
                }
            }
        }
        my $new_top_node = $node_change{$top_node};
        $new_edge{$new_top_node} = [@{$edge->{$top_node}}];
        $new_node_yh{$new_top_node} = $node_yh->{$top_node};
        %{$edge} = %new_edge;
        %{$node_yh} = %new_node_yh;
    }
    my ($root_len,$root_B) = ($nhx && $nhx=~/:(.+?)\[(.+)\]$/) ? ($1, $2) : (0, 0);
    ($root_len,$root_B,$no_blen);
}
sub draw_cycle_tree{
    my ($svg,$cx,$cy,$r1,$r2,$root,$root_len,$root_B,$x_unit,$y_unit,$tree,$spe,$node_xw,$node_yh,$fsize,
        $lw,$to_end,$color1,$color2,$color3,$colors,$edge,$tran,$opacity,$dot,$noname) = @_;
    ($tree && @$tree) || return(0);
    $tran ||= 0;
    my $add = 0;
    $lw ||= 1;
    $color1 ||= 'black'; #line
    $color2 ||= 'black'; #species name
    $color3 ||= 'blue'; #Boot
    my ($r0,$a0) = ($r1,$tran+$node_yh->{"$#$tree\->0"}*$y_unit);
    if($root_len){
        ($r0,$a0) = arc_path($svg,'d',[$cx,$cy,$r1,$a0,-1,$add,$root_len*$x_unit,0,$lw,$root_B,$fsize,$color3],
            'stroke',$color1,'stroke-width',$lw);
    }elsif($root){
        ($r0,$a0) = arc_path($svg,'d',[$cx,$cy,0.65*$r1,$a0,-1,$add,0.35*$r1,0,$lw,$root_B,$fsize,$color3],
            'stroke',$color1,'stroke-width',$lw);
    }
    $node_xw->{"$#$tree\->0"} = [$r0,$a0];
    for my $i(reverse (0..$#$tree)){
        for my $j(0..$#{$tree->[$i]}){
            defined $node_xw->{"$i->$j"} || die"error: x $i->$j\n";
            defined $node_yh->{"$i->$j"} || die"error: y $i->$j\n";
            ($r0,$a0) = @{$node_xw->{"$i\->$j"}};
            for my $branch(@{$tree->[$i]->[$j]}){
                my ($spe_name,$Bvalue,$L2);
                my $L = $branch->[1]*$x_unit;
                if($branch->[0]=~/^(\d+)$/ && $spe->[$1-1]){
                    $spe_name = $spe->[$1-1];
                    if($to_end && $r2-$r0-$L){
                        $L2 = $r2-$r0-$L;
                    }
                }
                if($branch->[2] && $branch->[2]=~/B=([\d\.]+)/){
                    $Bvalue = $1;
                }
                my ($line_color,$spe_color) = ($color1, $color2);
                if($colors && @$colors && $edge->{$branch->[0]}){
                    $line_color = ($spe_color = branch_color($edge->{$branch->[0]},$colors));#sub2.1
                }
                my ($r1,$a1) = arc_path($svg,'d',[$cx,$cy,$r0,$a0,$tran+$node_yh->{$branch->[0]}*$y_unit,$add,
                    $L,$L2,$lw,$Bvalue,$fsize,$color3,$spe_name,$line_color,$spe_color,$opacity,$dot,$noname],'stroke',$line_color,
                    'fill','none','stroke-width',$lw);
                $node_xw->{$branch->[0]} = [$r1,$a1];
            }
        }
    }
}
sub arc_path{
    my $svg = shift;
    my %path = @_;
    my ($cx,$cy,$r,$a1,$a2,$add,$L,$L2,$lw,$Bvalue,$fsize,$color3,$name,$color1,$color2,$opacity,$dot,$noname) = @{$path{d}};
    my @xy = (cycle_pos($cx,$cy,$r,$a1,$add), cycle_pos($cx,$cy,$r,$a2,$add));
    if($a2 == -1){
        delete $path{d};
        @xy[2,3] = cycle_pos($cx,$cy,$r+$L,$a1,$add);
        simpnum(@xy[0..3],2);
        $svg->line('x1',$xy[0],'y1',$xy[1],'x2',$xy[2],'y2',$xy[3],%path);
        path_name($svg,$Bvalue,$a1,$xy[2],$xy[3],$fsize,$color3,'end');
        return($r+$L,$a1-$add);
    }
    my $large = ($a2-$a1>0.5) ? 1 : 0;
    my $sweep = ($a2>$a1) ? 1 : 0;
    if($L){
        push @xy,(cycle_pos($cx,$cy,$r+$L,$a2,$add));
        $path{d} = "M @xy[0,1] A $r $r, 0, $large, $sweep, @xy[2,3] L @xy[4,5]";
    }else{
        $path{d} = "M @xy[0,1] A $r $r, 0, $large, $sweep, @xy[2,3]";
        $L = 0;
    }
    if($L2){
        $L += $L2;
        push @xy, (cycle_pos($cx,$cy,$r+$L,$a2,$add));
        simpnum(@xy[-4,-3,-2,-1],2);
        $svg->line('x1',$xy[-2],'y1',$xy[-1],'x2',$xy[-4],'y2',$xy[-3],'stroke',$color1,'stroke-width',$lw/2,
                'stroke-dasharray',"3 2");
    }
    $svg->path(%path);
    $noname || path_name($svg,$name,$a2,$xy[-2],$xy[-1],$fsize,$color2,$dot);
    path_name($svg,$Bvalue,$a1,$xy[0],$xy[1],$fsize,$color3,'end');
    $L && $opacity && annulus_path($svg,'d',[$cx,$cy,$r,$r+$L,$a1,$a2,$add],'fill-opacity', $opacity,'stroke','none','fill',$color1);
    ($r+$L,$a2-$add);
}
sub path_name{
    my ($svg,$name,$a2,$cx,$cy,$fsize,$color,$dot,$text_anchor) = @_;
    $name || return(0);
    $color ||= 'black';
    $text_anchor ||= 'start';
    my $rra = 360*$a2 - 90;
    my $rraa = 2*3.1415926*$a2;
    simpnum($cx,$cy,2);
    $dot && $svg->circle('cx',$cx,'cy',$cy,'r',$fsize/3,"fill",$color,"stroke",'none');
    my ($tx,$ty) = ($text_anchor ne 'end') ? ($cx+$fsize*sin($rraa)/2, $cy-$fsize*cos($rraa)/2) :
        ($cx-$fsize*(sin($rraa)+cos($rraa))/2, $cy-$fsize*(sin($rraa)-cos($rraa))/2);
    my $group = $svg->group("transform","rotate($rra,$tx,$ty)",'stroke','none', 'fill',$color,'font-family','Arial');
        $group->text('x',$tx,'y',$ty+$fsize/3,'-cdata',$name,'font-size',$fsize,'text-anchor',$text_anchor);
}
#sub
#==============#
sub annulus_path
#==============#
{
    my $svg = shift;
    my %path = @_;
    my ($cx,$cy,$r1,$r2,$a1,$a2,$add) = @{$path{d}};
    $add ||= 0;
    my $large = ($a2-$a1>0.5) ? 1 : 0;
    my $sweep = ($a2>$a1) ? 1 : 0;
    my $invsweep = 1-$sweep; 
    my @xy = (cycle_pos($cx,$cy,$r1,$a1,$add), cycle_pos($cx,$cy,$r1,$a2,$add));
    if($r2){
        push @xy,(cycle_pos($cx,$cy,$r2,$a2,$add), cycle_pos($cx,$cy,$r2,$a1,$add));
        $path{d} = "M @xy[0,1] A $r1 $r1, 0, $large, $sweep, @xy[2,3] L @xy[4,5] A $r2 $r2, 0, $large, $invsweep, @xy[6,7] Z";
    }else{
        $path{d} = "M @xy[0,1] A $r1 $r1, 0, $large, $sweep, @xy[2,3] L $cx $cy Z";
    }
    $svg->path(%path);
    $r2 ? ([@xy[0..3]],[@xy[6,7,4,5]]) : @xy;
}
sub cycle_pos{
    my ($cx,$cy,$r,$a,$add) = @_;
    my $pi = 3.14159265358979;
    $add ||= 0;
    my $aa = 2*$pi*($a+$add);
    my $x = sprintf("%.5f",$cx + $r*sin($aa));
    my $y = sprintf("%.5f",$cy - $r*cos($aa));
    foreach($x,$y){s/\.0+$//;}
    ($x, $y);
}


#sub2
sub draw_tree{
    my ($svg,$flank_x,$flank_y,$root,$root_len,$root_B,$x_unit,$y_unit,$tree,$spe,$node_xw,$node_yh,$fsize,
        $max_len,$lw,$to_end,$type,$color1,$color2,$color3,$colors,$edge,$inner_fill,$opacity,$noname) = @_;
    ($tree && @$tree) || return(0);
    $flank_y += $y_unit/2;
    my $ry = $flank_y + $node_yh->{"$#$tree->0"}*$y_unit;
    if($root_len){
        $svg->line('x1',$flank_x,'y1',$ry,'x2',($flank_x += $root_len*$x_unit),'y2',$ry,'stroke',$color1,'stroke-width',$lw);
    }elsif($root){
        $svg->line('x1',$flank_x/2,'y1',$ry,'x2',$flank_x,'y2',$ry,'stroke',$color1,'stroke-width',$lw);
    }
    if($root_B){
        my $By = $flank_y + $node_yh->{"$#$tree\->0"} * $y_unit;
        $svg->text('x',$flank_x-$fsize/2,'y',$By-$fsize/6,'stroke','none','fill',$color3,'-cdata', $root_B,
                'font-size', $fsize,'text-anchor','end','font-family','Arial');
    }
    $node_xw->{"$#$tree\->0"} = $flank_x;
    for my $i(reverse (0..$#$tree)){
        for my $j(0..$#{$tree->[$i]}){
            defined $node_xw->{"$i->$j"} || die"error: x $i->$j\n";
            defined $node_yh->{"$i->$j"} || die"error: y $i->$j\n";
            my ($root_x,$root_y) = ($node_xw->{"$i->$j"}, $node_yh->{"$i->$j"});
            for my $branch(@{$tree->[$i]->[$j]}){
                draw_branch($svg,$flank_y,$root_x,$root_y,$x_unit,$y_unit,$branch,$spe,$node_xw,$fsize,
                $max_len,$lw,$to_end,$type,$color1,$color2,$color3,$colors,$edge,$inner_fill,$opacity,$noname);#sub2.2
            }
        }
    }
}
#sub2.1
#================
sub branch_color{
#================
    my ($edge,$colors) = @_;
    my %colh;
    for(@$edge){$colh{$_}++;}
    my @k = keys %colh;
    (@k == 1) && return("rgb($colors->[$k[0]])");
    my $tol = @$edge;
    my @rgb;
    for (@k){
        $colors->[$_] || die"$_";
        $colors->[$_]=~/(\d+),(\d+),(\d+)/;
        $rgb[0] += $1 * $colh{$_} / $tol;
        $rgb[1] += $2 * $colh{$_} / $tol;
        $rgb[2] += $3 * $colh{$_} / $tol;
    }
    for(@rgb){$_ = int($_+0.5);}
    return("rgb(".join(",",@rgb).")");
}
#sub2.2
#===============
sub draw_branch{
#===============
    my ($svg,$flank_y,$root_x,$root_y,$x_unit,$y_unit,$branch,$spe,$node_xw,$fsize,
        $max_len,$lw,$to_end,$type,$color1,$color2,$color3,$colors,$edge,$inner_fill,$opacity,$noname) = @_;
    $lw ||= 1;
    $opacity ||= 0.3;
    $color1 ||= 'black'; # line color
    $color2 ||= 'black'; # name color
    $color3 ||= 'blue'; # boot value color
    $root_y = $flank_y + $root_y * $y_unit;
#   branch -> [node_name branch_length NHX y-corrdinate]
    if($colors && @$colors && $edge->{$branch->[0]}){
        $color1 = ($color2 = branch_color($edge->{$branch->[0]},$colors));#sub2.1
    }
    my ($bx, $by) = ($root_x+$branch->[1]*$x_unit, $flank_y+$branch->[3]*$y_unit);
    my $sweep = 1;
    if($type==2){
        simpnum($root_x,$root_y,$bx,$by,2);
        $svg->line('x1',$root_x,'y1',$root_y,'x2',$bx,'y2',$by,'stroke',$color1,'stroke-width',$lw);
    }elsif($type==1 || $type==3){
        my $ra = $bx - $root_x;
        my $rb = $root_y - $by;
        if($rb < 0){
            $rb = -$rb;
            $sweep = 0;
        }
        if($type==3){
            my $cyclen = $ra > $rb ? $rb/6 : $ra/6;
            ($cyclen > $fsize) && ($cyclen = $fsize);
            my $rx1 = $root_x + ($ra ? ($bx-$root_x)*$cyclen/$ra : 0);
            my $ry1 = $by + ($root_y-$by)*$cyclen/$rb;
            simpnum($root_x,$root_y,$ry1,$rx1,$bx,$by,$cyclen,2);
            $svg->path(d=>"M $root_x $root_y L $root_x $ry1 A $cyclen $cyclen 0 0 $sweep $rx1,$by L $bx $by",
                style=>{stroke=>$color1,fill=>'none','stroke-width'=>$lw});
        }else{
            simpnum($root_x,$root_y,$bx,$by,2);
            $svg->path(d=>"M $root_x $root_y A $ra,$rb 0 0 $sweep $bx,$by",
                style=>{stroke=>$color1,fill=>'none','stroke-width'=>$lw});
        }
    }else{
#        $svg->line('x1',$root_x,'y1',$root_y,'x2',$root_x,'y2',$by,'stroke',$color1,'stroke-width',$lw);
#        $svg->line('x1',$root_x,'y1',$by,'x2',$bx,'y2',$by,'stroke',$color1,'stroke-width',$lw);
        simpnum($root_x,$root_y,$bx,$by,2);
        $svg->path(d=>"M $root_x $root_y V $by H $bx",style=>{'stroke',$color1,fill=>'none','stroke-width',$lw});
    }
    if($branch->[0]=~/^(\d+)$/ && $spe->[$1-1]){
        my $spe_name = $spe->[$1-1];
        if($to_end && $max_len - $bx){
            simpnum($max_len,$bx,$by,2);
            $svg->line('x1',$bx,'y1',$by,'x2',$max_len,'y2',$by,'stroke',$color1,'stroke-width',$lw/2,
                    'stroke-dasharray',"3 2");
            $bx = $max_len;
        }
        $noname ||
        $svg->text('x',$bx+$fsize/2,'y',$by+$fsize/3,'stroke','none','fill',$color2,'-cdata', $spe_name,
                'font-size', $fsize,'text-anchor','start','font-family','Arial');
        if($inner_fill && ($inner_fill == 2)){
            my ($by1,$by2) = ($by-$y_unit/2, $by+$y_unit/2);
            simpnum($root_x,$bx,$by1,$by2,2);
            $svg->path(d=>"M $root_x $by1 V $by2 H $bx V $by1",'stroke','none',fill=>$color1,'fill-opacity',$opacity);
        }
    }else{
        $node_xw->{$branch->[0]} = $bx;
#        $inner_fill && $inner_fill==2 && same_arr($edge->{$branch->[0]},$colors) &&
#            $svg->path(d=>"M $root_x $root_y V $by H $bx V $root_y",'stroke','none',fill=>$color1,'fill-opacity',$opacity);
    }
    $inner_fill && ($inner_fill == 1) && simpnum($root_x,$root_y,$bx,$by,2);
    $inner_fill && ($inner_fill == 1) && 
        $svg->path(d=>"M $root_x $root_y V $by H $bx V $root_y",'stroke','none',fill=>$color1,'fill-opacity',$opacity);
    if($branch->[2] && $branch->[2]=~/B=([\d\.]+)/){
        my $By = ($type && $sweep == 0) ? $by+5*$fsize/6 : $by-$fsize/6;
        my $Bx = $bx-$fsize/2;
        simpnum($Bx,$By,$fsize,2);
        $svg->text('x',$Bx,'y',$By,'stroke','none','fill',$color3,'-cdata', $1,
                'font-size', $fsize,'text-anchor','end','font-family','Arial');
    }
}
sub same_arr{
    my ($arr,$color) = @_;
    my $fa = $color->[$arr->[0]];
    for (@$arr){
        ($fa ne $color->[$_]) && return(0);
    }
    1;
}
#sub3
#=============#
sub draw_scal{
#=============#
    my ($svg,$max_branch,$x_unit,$width,$scal_x,$scal_y,$fsize,$lw,$title,$one) = @_;
    my ($s_unit,$num) = axis_split($max_branch);
    ($num > 6) && ($num/=2, $s_unit*=2);
    my $flen = 0;
    for(0..$num){
        my $temp_scal = length($_*$s_unit);
        ($temp_scal > $flen) && ($flen = $temp_scal);
    }
    my $max_fsize = $width*2/($num-1)/($flen+2);
    ($max_fsize < $fsize) && ($fsize = $max_fsize);
    if($title){
        $scal_y += $fsize;
        my $tsize = 0.75*2*$width/length($title);
        ($tsize > $fsize) && ($tsize = $fsize);
        $svg->text('x',$scal_x,'y',$scal_y,'stroke','none','fill','black','-cdata',$title,
            'font-size', $tsize,'text-anchor','start','font-family','Arial');
    }
    $scal_y += 5*$fsize/5;
    my $ss_unit = $s_unit / 10;
    my $sx_unit = $ss_unit * $x_unit;
    my $x_end = $one ? $scal_x+10*$sx_unit : $scal_x+$width;
    $svg->line('x1',$scal_x,'y1',$scal_y,'x2',$x_end,'y2',$scal_y,'stroke','black','stroke-width',$lw);
    my $scal = 0;
    my $sx = $scal_x;
    my $i = 0;
    while($scal <= $max_branch){
        my $sy = ($i % 5) ? $scal_y-$fsize/4 : $scal_y-$fsize/2;
        $svg->line('x1',$sx,'y1',$scal_y,'x2',$sx,'y2',$sy,'stroke','black','stroke-width',$lw);
        if($i % 10 == 0){
            $svg->text('x',$sx,'y',$scal_y+$fsize,'stroke','none','fill','black','-cdata',$scal,
                'font-size', $fsize,'text-anchor','middle','font-family','Arial');
            $one && ($i == 10) && last;
        }
        $i++;
        $scal = $ss_unit * $i;
        $sx += $sx_unit;
    }
}

#sub3.1
#=============#
sub axis_split
#=============#
{
    #useage: axis_spli(a[,b,c]),a is the max value in the axis scale polt
    #b is the ratio of the max value to the length of the Y axis.
    #c is for precision, often use 2,4,8,16
    my ($maxv,$rate,$preci) = @_;
    $rate ||= 0.9;
    die"the ratio must between 0.5 to 0.98, if not you should revise your figure.\n"  if ($rate > 0.98 || $rate < 0.5);
    $preci ||= 2;$preci = 2** $preci;
    my ($mbs,$mag) = split/\.?\d?e/,sprintf("%1.1e", $_[0]);# the MSB of the max value in the plot
    $mag =~ s/^\+//;                                        # the order of magnitude of the max value in the plot.
    $mag =~ s/^0*//;
    $mag ||= 0;
    my $k = $rate / (1 - $rate) / $preci;              # the middle value used to caclutate $min_value-
                       # -you can also change preci into 2 or 1, the y-scal will become more precision
    my $min_value;     # the min value show in y axis
    foreach(2,1,0.5,0.25,0.125,0.1,0.05){
        $min_value = $_;
        ($mbs >= $_ * $k) && last;
    }
    $min_value = $min_value * 10**$mag;
    my $value_number = int($maxv / $min_value );  # the number of value show in y axis
    ($maxv > $value_number * $min_value) && ($value_number++);
    while(($value_number-1)*$min_value >= $maxv){$value_number--;}
    ($min_value, $value_number);
}
