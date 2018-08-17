#!/usr/bin/perl  -w
use FindBin qw($Bin);
use lib "$Bin";
use PATHWAY;
(-s "$Bin/Pathway_cfg.txt") || die"error: can't find config at $Bin/, $!\n";
BEGIN { 
    my( $svg_lib)= get_pathway("$Bin/Pathway_cfg.txt",[qw(SVG_Lib)]);
    unshift @INC, $svg_lib 
};
use SVG;
=head1 Name

bar_diagram.pl

=head1 Description

draw the bar diagram of SVG

=head1 Version

 Author: Wenbin Liu, liuwenbin@genomics.org.cn
 Version: 1.0,  Date: 2010-04-24
 Version: 2.0,  Date: 2010-11-27  Update: add option --sig_pos,dge,right,nosign,blank
 Version: 2.1,  Date: 2011-01-17  Update: add option --tranxy,preci
 Version: 2.2,  Date: 2011-02-10  Do some improve in the advice of Liang Xinming
 Version: 2.3,  Date: 2011-07-25  Update: add option -interval,stroke
 Version: 2.4,  Date: 2013-04-08  Update: add option -end_mask --ordercols
                                   inflie also can be a list of number instead of a file like: "x1 x2 x3 .. xn;y1 y2 y3 .. yn"

=head1 Usage

 perl bar_diagram.pl  <infile>  [--Option --help] >out.svg
 infile             input file store the data to draw diagram
 --Option:
 1 about paper size
 --width <num>      the length of the x axis, default=400
 --height <num>     the length of the y axis, default=300
 --barw <num>       the width of a bar, failure when --width set
 --flank_y <num>    distance from y axis to edge of drawing paper, default=height/3
 --flank_x <num>    distance from x axis to edge of drawing paper, default=flank_y
 
 2 about symbol
 --single           the inflie just be date for draw, not contaion sign or xlab scale
 --nohead <num>     the head specify rows not descirbe data, x-area date get from first rank
 --nosign           do not draw the symbol
 --table            infile in table form
 --ranks <str>      only show specify ranks when -table, ranks star from 0, default show all ranks
 --row <num>        max row the get at infile(head not concluded), default use all row.
 --symbol <str>     symbol text split by ',', default get form head info
 --size_sig <num>   figure symbol sign text font-size, default=0.05*height
 --sig_pos <str>    the position of the symbol star coordinate x,y, default auto
 --right            put the sign at the right of figure, or put at head
 --column <num>     the column number of symbol when symbol at head, default=2;
 
 3 about title and scale
 --h_title <str>    the head title, default no title
 --y_title <str>    the title of ylab, default no title
 --x_title <str>    the tilte of xlab, default no title
 --h_pos <str>      the middle position of the head title, default auto
 --size_ht <num>    the head title font-size, default=0.09*height
 --size_yt <num>    y-area title font-size, default=0.08*height
 --size_xt <num>    xlab title font-size, default = size_yt
 --size_ys <num>    y-area scale font-size, default=0.05*height
 --size_xs <num>    x-area text font-size, default=0.035*width
 --edge             put the x-area text at the edge of two group, default in middle
 --interval         to show x-area text interval
 --inter_num <num>  interval number, default=1
 --micro_scale      show micro scale
 --preci <num>      y-area scale precision option, select from 1..4, default=3
 --y_mun <str>      specified y-area scale unit,number, default not set
 --show_data        show the data number at figure
 --end_mask <num>   the end number of rank to be mask text for --show_data
 --rotate <str>     to rotate x-area scale text at specify angle
 --rotate2 <str>    to rotate data number  text at specify angle
 --log <flo>        to make data=log(data)/log(specify number), defualt not use
 --vice             draw vice y-axis bar, use only tow group
 --y_title2         the title of vice ylab, default no title
 --cluster          to cluster x-axis
 --cluhor           not to rotate y_title while -cluster -tranxy

 
 4 about bar
 --style <1|2|3|4>  bar style: bars in a group 1 heap up, 2 abreast, 3 mixed; default=1
 --pair <flo>       set overlap rate to put the nearly group of bars overlap
 --tranxy           to exchange x-y area, but x y title will not change
 --blank <fig>      the space between two group, 1 is width of one bar, default=1
 --colors <str>     the colours of the bar fill, splited by ',', default auto
 --ordercols <str>  just use sortorder colors, not --colors
 --stroke <str>     colors for bar stroke, default equal to fill color
 --help             output help information to screen

=head1 Notice

 1 The from of infile to be: first line x-axis name group, second line graphic symbol,
   other line eatch data of the x-axis group turn to the symbol.
 2 The symbol could have blank in second line of infile must split by "\t" or ';'
 3 when --tranxy, width and height set will exchange
 4 When width and barw not set but budget barw < 10, than barw will set to be 10
 5 The out svgfile must end by suffix ".svg".

=head1 Example

 perl bar_diagram.pl   orthologs.statdate > out.svg
 perl bar_diagram.pl "1 2 3 4 5;1 4 9 16 26" > out.svg

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin);
my ($height,$width, $flank_x,$flank_y, $size_x,$micro_scale,
    $size_xt, $size_scale, $size_sig, $tranxy,$size_yt,
    $ylab_title,$xlab_title, $style,$color,$sig_pos,$minus_xy,
    $blank,$edge,$nosign, $column, $right,$h_title, $no_xscale,
    $preci,$help,$h,$barw,$y_mun,$h_pos,$size_ht,$row,
    $nohead,$Symbol,$interval,$stroke,$inter_num,$table,$cluhor,
    $ranks,$grid,$rotate,$single,$log,$show_data,$rotate2,$pair,
    $vice,$y_title2,$y_mun2,$add_xblock,$add_yblock,$scale_extend,
    $end_mask,$ordercols,$opacity,$zoom,$cluster,$grid2,$aa_flank_x,
    $endmask,$textup,$width_lim,$colsel,$avg,$scalrate,$rev_sym,$sr
);
my @add_line;
GetOptions(
           "height:i"   => \$height,        "vice"       => \$vice,
           "width:i"    => \$width,         "y_title2:s" => \$y_title2,
           "flank_x:i"  => \$flank_x,       "interval"   => \$interval,
           "flank_y:i"  => \$flank_y,       "colors:s"   => \$color,
           "barw:i"     => \$barw,          "stroke:s"   => \$stroke,
           "inter_num:i"=> \$inter_num,     "right"      => \$right,
           "nohead:i"	=> \$nohead,        "nosign"     => \$nosign,
           "symbol:s"	=> \$Symbol,        "blank:f"    => \$blank,
           "size_xs:f"  => \$size_x,        "column:i"   => \$column,
           "size_xt:f"  => \$size_xt,       "sig_pos:s"  => \$sig_pos,
           "size_yt:f"  => \$size_yt,       "h_pos:s"    => \$h_pos,
           "tranxy"     => \$tranxy,        "size_ht:f"  => \$size_ht,
           "preci:i"    => \$preci,         "h_title:s"  => \$h_title,
           "y_mun:s"	=> \$y_mun,         "rotate:i"   => \$rotate,
           "size_ys:f"  => \$size_scale,    "rotate2:i"  => \$rotate2,
           "size_sig:f" => \$size_sig,      "help"       => \$help,
           "y_title:s"  => \$ylab_title,    "h"          => \$h,
           "x_title:s"  => \$xlab_title,    "table"      => \$table,
           "style:i"    => \$style,         "ranks:s"    => \$ranks,
           "edge"       => \$edge,          "row:i"      => \$row,
           "grid"       => \$grid,          "y_mun2:s"    => \$y_mun2,
           "single"     => \$single,        "add_xblock:i" => \$add_xblock,
           "log:f"      => \$log,           "add_yblock:i" => \$add_yblock,
           "show_data"  => \$show_data,     "scale_extend" => \$scale_extend,
           "micro_scale"=> \$micro_scale,   "minus_xy:s"    => \$minus_xy,
           "no_xscale"  => \$no_xscale,     "end_mask:i"    => \$end_mask,
           "ordercols:s"=> \$ordercols,     "opacity:f"     => \$opacity,
           "zoom:f"     => \$zoom,          "cluster"       => \$cluster,
           "grid2"      => \$grid2,         "add_line:s"    => \@add_line,
           "aa_flank_x:i"=> \$aa_flank_x,   "pair:f"        => \$pair,
           "cluhor"     => \$cluhor,        "endmask:i"     => \$endmask,
           "textup"   => \$textup,          "width_lim:i"   => \$width_lim,
           "colsel:s" => \$colsel,          "avg"           => \$avg,
           "scalrate:f" => \$scalrate,      "rev_sym"       => \$rev_sym,
           "sr:f"       => \$sr
);
$help && (die `pod2text $0`);
if (@ARGV != 1 || $h){
    die "Name: bar_diagram.pl
Author: Wenbin Liu, liuwenbin\@genomics.org.cn
Version: 2.3,  Date: 2011-07-25
Usage: perl bar_diagram.pl  <infile>  [--Option --help] >out.svg
 <infile>          the file store the data in the diagram
 -nosign           do not draw the symbol
 -nohead           the head 2rows not descirbe data, x-area date get from first rank
 -table            infile in table form
 -symbol <str>     symbol text split by ',', default get form head info
 -sig_pos <str>    the symbol star area 'x,y', default auto
 -h_title <stre>   the head title, default no title
 -y_title <str>    the title of ylab, default no title
 -x_title <str>    the tilte of xlab, default no title
 -style <1|2|3|4>  bar sytle: bars in a group 1 heap up, 2 abreast, 3 mixed; default=1
 -tranxy           to exchange x-y area, but x y title will not change
 -blank <fig>      the space between two group, 1 is width of one bar, default=1
 -colors <str>     the colours of the bar, split by ',', default auto
 -stroke <str>     colors for bar stroke, default equal to fill color
 -edge             put the x-area text at the edge of two group, default in middle
 -right            put the sign at the right of figure, or put at head
 -column <num>     the column number of symbol when symbol at head, default=2;
 -preci <num>      y-area scale precision option, select from 1..4, default=3
 -sr <flo>         the y-area Zoom Scale, default not set
 -h                output brief help information to screen
Note: You can use --help to get detail help information\n";
}

#************************************************************#
#
#                          MAIN                              #
#
#************************************************************#
my $infile = shift;
my @sel = $ranks ? split/,/,$ranks : ();
my (@x_group,@symbol,$sample_name);
my (@y_group,@ys_group,@endmask,@ordercol);
my ($tol_amin,$max_y,$max_y2) = (0) x 3;
my $del_lim = $nohead ? 2 : 1;
$ordercols && (@ordercol = split/[,\s]+/,$ordercols);
$scalrate ||= 0.95;
$style ||= 1;
($style =~ /\b[1234]\b/) || die"Error: --style noly can be 1 or 2 or 3\n";
($style == 3) && ($opacity ||= 0.5);
(defined $pair) && ($opacity ||= 0.8);
if(-s $infile){ ####### =>
open IN, $infile || die$!;
if($single){
#    $nosign = 1;
}elsif($nohead){
   for(1..$nohead){
       <IN>;
   }
}else{
	chomp(my $fline = <IN>);
	@x_group = ($fline=~/;/) ? (split/;/,$fline) : ($fline=~/\t/) ? (split/\t/,$fline) : (split/\s+/,$fline);
    if($table){
        $sample_name = shift @x_group;
        @symbol = @x_group;
        $ranks && (@symbol = @symbol[@sel]);
        @x_group = ();
        $nohead = 1;
    }else{
	    chomp(my $sline = <IN>);
	    ($sline =~ /\S/) ? (@symbol = split /;|\t/, $sline): ($nosign = 1);
    }
}
$inter_num ||= ($interval ? 2 : 1);
($inter_num > 1) && ($interval = 1);
$del_lim = $nohead ? 2 : 1;
my $line_num = 0;
if($row && $row<0){
    my $fl = (split/\s+/,`wc -l $infile`)[0];
    $fl -= $table ? 1 : defined $nohead ? $nohead : 2;
    $row += $fl;
}
while (<IN>){
	/\d/ || next;
	chomp;
	my @l = (/\t/) ? split /\t+/ : split /\s+/;
    (@l < $del_lim) && next;
    $ranks && (@l = @l[@sel]);
    my $head_tax;
    ($nohead || $single) && ($head_tax = shift @l);
    $line_num && ($line_num > @l) && next;
    ($l[-1] =~ /^\d/) || next;
    $line_num = @l;
	$head_tax && (push @x_group,$head_tax);
    $end_mask && (push @endmask,[splice(@l,-$end_mask)]);
    (defined $endmask) && (push @endmask,[$l[$endmask]]);
	$tol_amin++;
    push @ys_group,[@l];
    if($log){
        foreach(@l){$_ &&= log($_)/log($log);}
    }elsif($zoom){
        foreach(@l){$_ &&= $_ * $zoom;}
    }
    push @y_group,[@l];
    if("@l" =~ /[\/\|\\:]/){
        foreach(@l){$_ = (split/[\/\|\\:]/)[0];}
    }
    if(@l>=2 && $vice){
        $style = 2;
        ($max_y < $l[0]) && ($max_y = $l[0]);
        ($max_y2 < $l[-1]) && ($max_y2 = $l[-1]);
    }else{
        $max_y = ($style==3) ? max($max_y,style_len(@l)) : ($style == 1) ? max($max_y, sum(@l)) : max($max_y, @l);
    }
    $row && ($tol_amin == $row) && last;
}
close IN;
(defined $endmask) && ($end_mask = 1);
}else{  ##### ==
    my ($inxg, $inyg) = split/;/,$infile;
    @x_group = split /[,\s]+/,$inxg;
    for(split /[,\s]+/,$inyg){
        /\d/ || next;
        my @xx = split /[\/\|\\:]/;
        push @y_group,[@xx];
        $max_y = ($style==3) ? max($max_y,style_len(@xx)) : ($style == 1) ? max($max_y, sum(@xx)) : max($max_y, @xx);
    }
    $tol_amin = @x_group;
} ##### <==
sub style_len{
    my $len = 0;
    my $add = 1;
    for(@_){
        $len += $add * $_;
        $add *= -1;
    }
    $len;
}
@y_group || exit;
$Symbol && (@symbol = split/,/,$Symbol);
my @avg;
if($avg){
    for my $iy(@y_group){
        for my $i(0..$#$iy){
            $avg[$i] += ($iy->[$i] ||= 0);
        }
    }
    my $y_num = @y_group;
    for (@avg){
        $_ = /\./ ? sprintf("%.2f",$_/$y_num) : int($_/$y_num+0.5);
    }
    for my $i(0..$#symbol){
        $symbol[$i] .= "(avg:$avg[$i])";
    }
}

@symbol || ($nosign = 1);
$blank  ||= (defined $pair) ? 0.3 : 1;
my $budget_bar_num = ((($style==2) ? @{$y_group[0]} : 1) + $blank) * $tol_amin;
($budget_bar_num > 80) && ($barw ||= 8);
($budget_bar_num > 50) && ($barw ||= 10);

##======        Global Variable       ======##
$width   ||= ($barw ? $barw*$budget_bar_num : $tranxy ? 300 : 600);
($width_lim && $width_lim < $width) && ($width = $width_lim);
$height  ||= $tranxy ? 450 :300;
$flank_y ||= $height / 3;
$flank_x ||= $flank_y;
$tranxy && ($flank_x *= 1.2);
$size_scale ||= max(0.03 * $height,15);    #################
$size_yt    ||= max(0.04 * $height,15);
$size_xt    ||= $size_yt;
$size_ht	||= max(0.05 * $height,15);
$tranxy && (($width,$height) = ($height,$width));
my $pwidth  = $width + $flank_x * 1.5;     # Calculate the width of the paper
my $pheight = $height + $flank_y * 1.6;    # Calculate the height of the paper
(!$tranxy && $rotate) && ($pheight += 0.4*$flank_y);
$add_xblock && ($pwidth += $add_xblock);
$add_yblock && ($pheight += $add_yblock);
($vice || $right) && ($pwidth += $flank_x/2);
if($minus_xy && $minus_xy =~ /(\S+),(\S+)/){
    $flank_y -= $2;
    $pwidth -= $1;
}
$vice && $tranxy && ($right = 1);
$size_x ||= $tranxy ? $size_scale : 0.03 * $width;  ###### The size of the x  axis title
my $xpolt_unit = ($tranxy ? $height : $width) / $tol_amin;    # The min distance in the x ordinate while drawing
($size_x > 0.8*$xpolt_unit) && ($size_x = 0.8*$xpolt_unit);
my $may_size_x;
my (@clusters,@xgroup_len);
my ($longest_x_str,$clulen) = longest_str(\@x_group,\@clusters,\@xgroup_len,$cluster,$cluhor);
my $clus_len = @clusters ? 4 : 0;
if($tranxy){
    $may_size_x = 2*$flank_x/($longest_x_str + 2 + $clus_len);
    if($may_size_x < $size_x){
       if($scale_extend){
           $pwidth += $size_x * ($longest_x_str+2+$clus_len)/2 - $flank_x;
           $flank_x = $size_x * ($longest_x_str+2+$clus_len)/2 ;
       }else{
           $size_x = $may_size_x;
       }
    }
}elsif($rotate){
    my $str_len = abs(($longest_x_str+$clus_len) * sin(3.1415926*abs($rotate)/180) )+2;
#    $may_size_x = 2*($pheight-$flank_y-$height)/$str_len;
#    ($may_size_x < $size_x) && ($size_x = $may_size_x);
    my $may_flank_y = 0.5*$str_len*$size_x;
    $xlab_title && ($may_flank_y += 2*$size_xt);
    if($pheight-$flank_y-$height < $may_flank_y){
        $pheight += $may_flank_y - ($pheight-$flank_y-$height);
    }
    $str_len = $size_x * abs( ($longest_x_str+$clus_len) * cos(3.1415926*abs($rotate)/180) )/2;
    if($str_len > $flank_x){
        my $may_flank_x = 0;
        if(!$clus_len){
            foreach my $xg(0..$#x_group){
                my $may_flank_x_temp = $size_x * abs( (length($x_group[$xg])+$clus_len) * cos(3.1415926*$rotate/180) )/2;
                if($rotate > 0){
                    $may_flank_x_temp -= $width*($#x_group-$xg)/@x_group;
                }else{
                    $may_flank_x_temp -= $width*$xg/@x_group;
                }
                ($may_flank_x_temp > $may_flank_x) && ($may_flank_x = $may_flank_x_temp);
            }
        }else{
            $may_flank_x = $size_x * (abs($longest_x_str * cos(3.1415926*$rotate/180) ) + 2+$clus_len) /2;
        }
        if($may_flank_x > $flank_x){
            $pwidth += ($may_flank_x - $flank_x);
            ($rotate > 0) || ($flank_x = $may_flank_x);
        }
    }
}else{
	$may_size_x = ($interval ? 1.5 : 2) * $width / (@x_group * $longest_x_str);
	($may_size_x < $size_x) && ($size_x = $may_size_x);
}
($size_scale < $size_x) && ($size_x = $size_scale);
$size_sig   ||= ($tol_amin > 5 ? 0.04 : 0.05) * $height;
$preci      ||= 3;
my @colors =
#  qw(darkblue cornflowerblue deeppink bisque lightgreen lightblue green yellow orange dodgerblue black);
#qw(crimson blue lightseagreen orange mediumpurple palegreen lightcoral dodgerblue lawngreen olive
#    yellow fuchsia salmon mediumslateblue darkviolet purple sienna  black  tan chocolate skyblue turquoise cadetblue);
qw(pink orange green cyan blue purple cornflowerblue red darkturquoise sienna bisque blueViolet orangered olive lightseagreen crimson salmon yellow fuchsia mediumslateblue darkviolet tan chocolate skyblue turquoise cadetblue navajowhite slategrey cornflowerblue royalblue dodgerblue deepskyblue mediumaquamarine lawngreen yellowgreen gold saddlebrown indianred deeppink darkred peachpuff);
if($colsel && $colsel=~/(\d+),(\d+)/){
    @colors = @colors[$1 .. $2];
}
if($color){
    if(-s $color){
        unshift @colors,split/\s+/,`less $color`;
    }else{
    	@colors = ((split/[,\s]+/,$color),@colors);
    }
}
if(@{$y_group[0]} > @colors && -s "$Bin/rgb_colors.txt"){
    unshift @colors ,split/\s+/,`awk '(\$1!=\"white\"){print \$1}' $Bin/rgb_colors.txt`;
}

#==

#=========================================================#
#                 Drawing The Bar graphs                  #
#=========================================================#

##=======   Draw  the X Y axis   =======##
##
if($vice && $xlab_title && $y_title2 && $tranxy){
    ($max_y, $max_y2) = ($max_y2, $max_y);
    ($xlab_title, $y_title2) = ($y_title2, $xlab_title);
}
my ($yvalue_min, $yvalue_number) = $y_mun ? (split/,/,$y_mun) : &axis($max_y, $scalrate, $preci,$sr);
# Use &axis to caclutate the min value and the number of value that show in the y axis
# Caclutate the space the value in the y axis take in the paper,one number's length is likely to bo 7;
##=======   Draw  the  Y axis   =======##
my $y_max = caculate_ymax($yvalue_number,$yvalue_min,$log);#sub5.1
my $yused_size = $size_scale * $y_max / 2;
if($yused_size + 2*$size_yt > $flank_x){
    $pwidth += ($yused_size + 2*$size_yt - $flank_x);
    $flank_x = $yused_size + 2*$size_yt;
}
if($aa_flank_x){
    $flank_x += $aa_flank_x;
    $pwidth += $aa_flank_x;
}
$pheight -= $flank_y / 2;
$flank_y /= 2;
my @signs = @symbol;
#$size_sig ||= 0.032 * $width;    ########################################
$column ||= 2;
my $splitn = int(@signs / $column) || 1;
(@signs > $column * $splitn) && ($splitn++);
my ($sig_x, $sig_y,$sig_y2,$size_sig_max);
my @sig_long;
if ($right){
    my $rig_str_len = longest_str(\@signs);
    $size_sig_max = ($pwidth-$width-$flank_x) / $rig_str_len;
    if($size_sig_max < 10){
        $pwidth += (10-$size_sig_max)*$rig_str_len;
        $size_sig_max = 10;
    }
    ($size_sig_max < $size_sig) && ($size_sig = $size_sig_max);
    ($sig_x, $sig_y) =($flank_x + $width + 1.2 * $size_sig, $flank_y + ($height - 1.2 * $size_sig * ($#signs + 1))/4);
    $sig_y2 = $flank_y + 1.2 * $size_sig;
}else{
    ($sig_x, $sig_y) = ($flank_x + $xpolt_unit / 4, $flank_y + 5 - 1.2 * $size_sig * ($splitn+1));
    ($sig_y < 5) && ($size_sig = ($flank_y-10)/1.2/($splitn+1), $sig_y=5);
    $sig_y2 = $sig_y;
    my $s = 0;
    foreach (1 .. $column){
        my $e = $s + $splitn - 1;
        ($e > $#signs) && ($e = $#signs);
        push @sig_long, (longest_str([@signs[$s .. $e]]) + 5);
        $s += $splitn;
    }
}

##======   Creat A New drawing paper  ======##
my $w1 = "http://www.w3.org/2000/svg";
my $w2 = "http://www.w3.org/1999/xlink";
my $svg = SVG->new(width=> $pwidth,height=> $pheight,xmlns => $w1,"xmlns:xlink" => $w2);
#==


my $ypolt_unit = draw_y_axis($flank_x,$flank_y,$width,$height,$size_scale,$yvalue_min,$yvalue_number,$grid,$micro_scale,$tranxy,$log,0);#sub5
my ($yvalue_min2,$yvalue_number2,$ypolt_unit2,$y_max2);
if($vice){
    ($yvalue_min2, $yvalue_number2) = $y_mun2 ? (split/,/,$y_mun2) : &axis($max_y2, $scalrate, $preci, $sr);
    $y_max2 = caculate_ymax($yvalue_number2,$yvalue_min2,$log);#sub5.1
    $ypolt_unit2 = draw_y_axis($flank_x,$flank_y,$width,$height,$size_scale,$yvalue_min2,$yvalue_number2,$grid,$micro_scale,$tranxy,$log,1);
}
##===== Draw the X Y axis title ======##
my $xused_size = $longest_x_str * $size_x / 2 + $size_x;
if($cluster){
    $xused_size += (1+$clulen/2) * $size_x;
}elsif($rotate){
    $xused_size = $size_x * abs( ($longest_x_str+$clus_len) * sin(3.1415926*abs($rotate)/180) )*0.56 + $size_x;
}else{
    $xused_size = $longest_x_str * $size_x / 2 + $size_x;
}
$tranxy && (($size_yt, $size_xt, $yused_size) = ($size_xt, $size_yt, $xused_size));
draw_xtitle($xlab_title,$size_xt,$flank_x,$flank_y,$width,$height,$rotate ? $xused_size : 2*$size_x,$tranxy);#sub6
draw_ytitle($ylab_title,$size_yt,$flank_x,$flank_y,$width,$height,$yused_size);#sub7
if($vice && $y_title2){
    $tranxy ? draw_xtitle($y_title2,$size_xt,$flank_x,$flank_y,$width,$height,2*$size_scale,1,1) :#sub6
        draw_ytitle($y_title2,$size_yt,$flank_x,$flank_y,$width,$height,$size_scale*$y_max2/2,1);
}
($vice && $xlab_title && $y_title2 && $tranxy) && (($ypolt_unit,$ypolt_unit2) = ($ypolt_unit2, $ypolt_unit));
#sub6
#===============
sub draw_xtitle{
#===============
    my ($xlab_title,$size_xt,$flank_x,$flank_y,$width,$height,$xused,$tranxy,$vice) = @_;
    ($xlab_title) || return(0);
    my $may_size_xt = 2*$width/length($xlab_title);
    ($size_xt > $may_size_xt) && ($size_xt = $may_size_xt);
    my $ttxx = $flank_x + $width / 2;
    my $ttxy = $vice ? $flank_y - $xused :
        $height +$flank_y + ($tranxy ? 2*$size_scale : $xused) + $size_xt;
    $svg->text('x',$ttxx,'y',$ttxy,'stroke','none','fill','black','-cdata',$xlab_title, 'font-size', $size_xt,
            'text-anchor', 'middle','font-family','Arial');
}
#sub7
#================
sub draw_ytitle{
#================
    my ($ylab_title,$size_yt,$flank_x,$flank_y,$width,$height,$yused_size,$vice) = @_;
    ($ylab_title) || return(0);
    my $may_size_yt = 2*$height/length($ylab_title);
    ($size_yt > $may_size_yt) && ($size_yt = $may_size_yt);
    my $r = $vice ? 90 : -90;
    my ($ytitle_x, $ytitle_y) = ($vice ? $flank_x+$width+$size_yt+$yused_size : $flank_x-$size_yt-$yused_size,$flank_y+$height/2);
#            $flank_y + 0.5 * $height - $size_yt * length($ylab_title) / 4);
    my $g = $svg->group("transform" => "rotate($r,$ytitle_x,$ytitle_y)");
#    ($ytitle_x, $ytitle_y) = $vice ? ($ytitle_y, -$ytitle_x) : (-$ytitle_y, $ytitle_x);
    $g->text('x',$ytitle_x,'y',$ytitle_y,'stroke','none', 'fill','black','-cdata',$ylab_title, 'font-size', $size_yt,'font-family','Arial','text-anchor','middle');
#            'font-family','Arial','text-anchor', 'end');
}

##=======   darw the bar diagram  =======##
my $sym_num = @{$y_group[0]};
my $split_unit_x = ($style =~ /[134]/) ? $xpolt_unit / (1 + $blank)
  : $xpolt_unit / ($sym_num + $blank);
my $x1 = ($tranxy ? $flank_y : $flank_x) + $blank * $split_unit_x / 2;
my ($y_0, $x_0);
$y_0 = $tranxy ? $flank_x : ($height + $flank_y);
my ($xtitle_x, $xtitle_y) = $tranxy ? 
($flank_x - $size_x , ($edge ? $flank_y : $flank_y + 0.5 * $xpolt_unit))
  : ($edge ? $flank_x : ($flank_x + 0.5 * $xpolt_unit),$flank_y + $height + $size_x);
my ($texanchor, $texx, $texy, $texh, $texw,$dotl_lim) = $tranxy ? ('end', 'y', 'x', 'width', 'height', $width/6)
  : ('middle', 'x', 'y', 'height', 'width', $height/6);
my %key_pos;
my @kxy;
if($cluster && @clusters){
    for (@clusters){
        $key_pos{$_->[1]} = 1;
        $key_pos{$_->[2]} = 2;
    }
}
my $x00 = $x1;
my $data_size = 0;
if($show_data){
    for my $d(@y_group){
        for (@$d){
            my $d_len = length;
            ($d_len > $data_size) && ($data_size = $d_len);
        }
    }
    $data_size = 2*( $style =~ /[134]/ ? $xpolt_unit : $split_unit_x)/($data_size+2);
}
foreach my $i (0 .. $tol_amin - 1){
    if(defined $pair){
        $x1 = $x00 + ($i % 2 ? -1 : 1) * ($blank/2 + $pair) * $xpolt_unit;
    }
    $x_0 = $x1;
    my $xline = ($tranxy ? $flank_y : $flank_x) + $xpolt_unit * ($i + 1);
    if(!$no_xscale){
        if($tranxy){
            $svg->line('x1', $flank_x, 'y1', $xline, 'x2', $flank_x - $size_x / 5,
                    'y2', $xline, 'stroke', 'black', 'stroke-width', 2);
            $grid2 && $svg->line('x1', $flank_x, 'y1', $xline, 'x2', $flank_x + $width,
                    'y2', $xline, 'stroke', 'black', 'stroke-width',1);
        }else{
            $svg->line('x1', $xline, 'y1', $y_0, 'x2', $xline, 'y2',
                   $y_0 + $size_x / 5,'stroke', 'black', 'stroke-width', 2);
            $grid2 && $svg->line('x1', $xline, 'y1', $y_0, 'x2', $xline, 'y2',
                    $y_0 - $height, 'stroke', 'black', 'stroke-width', 1);
        }
    }
    my $xalia_title;
    my ($t_x, $t_y) = ($xtitle_x,$xtitle_y + ($tranxy ? $size_x / 3 : 0));
    my $t_y2 = $t_y;
    $rotate && ($t_y2 -= $size_x/3);
    my $group = $rotate ? $svg->group("transform","rotate($rotate,$t_x,$t_y2)",'stroke','none', 'fill','black','font-size',
        $size_x,'font-family','Arial','text-anchor', $rotate > 0 ? 'start' : 'end') :
    $svg->group('stroke','none', 'fill','black','font-size',$size_x,'font-family','Arial','text-anchor', $texanchor);

    if($interval && !($i % $inter_num)){
    	$xalia_title = $x_group[$i/$inter_num];
    	$group->text('x',$t_x,'y',$t_y,'-cdata',$xalia_title);
    }elsif(!$interval){
    	$xalia_title = $x_group[$i];
    	$group->text('x',$t_x,'y',$t_y,'-cdata',$xalia_title);
    }
    if($key_pos{$i}){
        if($tranxy || $rotate){
            my $raf = $rotate ? 3.1415926 * $rotate / 180 : 0; 
            @kxy[0..3] = ($t_x - $size_x * ($xgroup_len[$i]+0.3) * cos($raf) * 0.6,
                        $t_y - $size_x * ($xgroup_len[$i]+0.3) * sin($raf) * 0.6,
                        $t_x - $size_x * ($longest_x_str + 0.3) * cos($raf) * 0.6,
                        $t_y - $size_x * ($longest_x_str + 0.3) * sin($raf) * 0.6);
            if($tranxy && !$rotate){
                $kxy[1] -= $size_x/3;$kxy[3] -= $size_x/3;
            }
        }else{
            @kxy[0..3] = ($t_x, $t_y+$size_x/2,$t_x,$t_y+$size_x);
        }
        $svg->line('x1',$kxy[0],'y1',$kxy[1],'x2',$kxy[2],'y2',$kxy[3],'stroke','black','stroke-width',1);
        if($key_pos{$i} == 2 && $kxy[5]){
            $svg->line('x1',$kxy[4],'y1',$kxy[5],'x2',$kxy[2],'y2',$kxy[3],'stroke','black','stroke-width',1);
            my ($kcx,$kcy) = ( ($kxy[2]+$kxy[4])/2, ($kxy[3]+$kxy[5])/2);
            if($tranxy){
                $kcy += $size_x/3;
                if($cluhor){
                    $kcx -= $size_x/2;
                    $group = $svg->group('stroke','none', 'fill','black','font-size',$size_x,'text-anchor','end');
                }else{
                    $kcx -= $size_x;
                    $group = $svg->group("transform","rotate(-90,$kcx,$kcy)",'stroke','none', 'fill','black','font-size',$size_x,'text-anchor','middle');
                }
            }else{
                $kcy += $size_x;
                $group = $svg->group('stroke','none', 'fill','black','font-size',$size_x,'text-anchor','middle');
            }
            my $kc_title = shift @clusters;
            $group->text('x',$kcx,'y',$kcy,'-cdata',$kc_title->[0],'font-family','Arial');
        }
        @kxy[4,5] = @kxy[2,3];
    }
    my $y1 = $y_0;
    my $t = 0;
    my $per_ty;
    my @tol_text;
    foreach my $j (0 .. $sym_num - 1){
        my $yu = ($vice && $j==$sym_num-1) ? $ypolt_unit2 : $ypolt_unit;
        my ($yh,$yh1,$yh2,$yh3,$yh4);
        my $dot_end = ($y_group[$i]->[$j] =~ s/\.$//) ? 1 : 0;
        my $opac_p = ($y_group[$i]->[$j] =~ s/[\/\|\\:]p/:/) ? 1 : 0;
        if($y_group[$i]->[$j]=~/(\S+)[\/\|\\:](\S+)/){
            ($yh,$yh2,$yh3) = ($yu*$1, $yu*$2, 1);
            $opac_p && (($yh3,$yh4) = (0, $2));
        }else{
            $yh = $yu * (${$y_group[$i]}[$j] || 0);
        }
        $yh1 = $y1;
        if($style == 3 && $j % 2){
            my $jo = ($j+3) % 4 ? 1 : -1;
            $tranxy && ($jo *= -1);
            $tranxy ? ($y1 -= $yh, $x1 -= $jo*$xpolt_unit/10) : ($y1 += $yh, $x1 += $jo*$xpolt_unit/10);
            $t++;
            next;
        }
        $j -= $t;
        my $stroke_color = ($stroke || $ordercol[$i*$sym_num+$j] || $colors[$j]);
        $tranxy || ($y1 = ($style =~ /[13]/) ? $y1 - $yh : ($y_0 - $yh));
        if($yh3){
            $tranxy || ($yh1 = ($style == 1) ? $yh1 - $yh2 : ($y_0 - $yh2));
            $yh && $svg->rect($texx,$x1,$texy, $y1,$texw,$split_unit_x,$texh,$yh,'stroke',$stroke_color, 'fill', 'none');
            $yh2 && $svg->rect($texx,$x1,$texy, $yh1,$texw,$split_unit_x,$texh,$yh2,'stroke', 'none', 'fill', $stroke_color);
        }else{
            my @svg_rect = ('stroke', $stroke_color, 'fill', $ordercol[$i*$sym_num+$j] || $colors[$j]);
            if($opac_p){
                push @svg_rect,('fill-opacity',$yh4||0);
            }elsif($opacity){
                push @svg_rect,('fill-opacity', $opacity);
            }
            $yh && $svg->rect($texx,$x1,$texy, $y1,$texw,$split_unit_x,$texh,$yh,@svg_rect);
        }
        if($dot_end){
            my ($dot_x,$dot_y,$dot_w,$dot_h,$dot_t,$dot_l);
            if($tranxy){
                $dot_l = (1.1*$flank_x+$width-$y1-$yh)/3;
                ($dot_l > $dotl_lim) && ($dot_l = $dotl_lim);
                ($dot_x,$dot_y,$dot_w,$dot_h,$dot_t) = ($y1+$yh+$dot_l/4,$x1,$dot_l,$split_unit_x,2);
            }else{
                $dot_l = ($y1-$flank_y)/3;
                ($dot_l < $flank_y/5) && ($dot_l = $flank_y/5);
                ($dot_l > $dotl_lim) && ($dot_l = $dotl_lim);
                ($dot_x,$dot_y,$dot_w,$dot_h,$dot_t) = ($x1,$y1-5*$dot_l/4,$split_unit_x,$dot_l,0);
            }
            arrow($svg,$dot_x,$dot_y,$dot_w,$dot_h,$dot_t,1.5,$stroke_color);
        }
        if($show_data || $end_mask){
            my $size_x2 = $size_x * 0.6;
            if($rotate2){
                if($tranxy && $size_x2 > $yh){
                    $size_x2 = $yh;
                }elsif($size_x2 > $split_unit_x){
                    $size_x2 = $split_unit_x;
                }
            }elsif($data_size && $size_x2 > $data_size){
                $size_x2 = $data_size;
            }
            if($tranxy){
                $t_x = $y1 + $yh + $size_x2/6;
                $t_y = $x1 + $split_unit_x/2;
            }else{
                $t_y = $y1-$size_x2/6;
                $t_x = $x1 + $split_unit_x/2;
            }
            $rotate2 && ($t_y2 = $t_y - $size_x2/3);

            my @gp = $rotate2 ? ('stroke','none', 'fill','black','font-size',$size_x2,'font-family','Arial','text-anchor', 'start') :
                ('stroke','none', 'fill','black','font-size',$size_x2,'font-family','Arial','text-anchor', 'middle');
            my $text = $end_mask ? $endmask[$i]->[$j] : $ys_group[$i]->[$j];
            if($textup && !$tranxy){
                push @tol_text,[$t_x,$t_y,$text,[@gp],$size_x2];
            }else{
                $rotate2 && (push @gp,("transform","rotate($rotate2,$t_x,$t_y2)"));
                my $gp = $svg->group(@gp);
                $gp->text('x',$t_x,'y',$t_y,'-cdata',$text);
            }
        }
        if($style =~ /[134]/){
            $tranxy  && ($y1 += $yh);
        }else{
            $x1 += $split_unit_x;
        }
    }
    if(@tol_text){
        @tol_text = sort {$b->[1]<=>$a->[1]} @tol_text;
        for my $txy(@tol_text){
            my $size_x2 = $txy->[4];
            if($per_ty && $per_ty-$size_x2 < $txy->[1]){
                $txy->[1] = $per_ty-$size_x2;
            }
            $per_ty = $txy->[1];
            $t_y2 = $txy->[1]-$size_x2/3;
            $rotate2 && (push @{$txy->[3]},("transform","rotate($rotate2,$txy->[0],$t_y2)"));
            my $gp = $svg->group(@{$txy->[3]});
            $gp->text('x',$txy->[0],'y',$txy->[1],'-cdata',$txy->[2]);
        }
    }
    $x1 = $x_0 + $xpolt_unit;
    $x00 += $xpolt_unit;
    $tranxy ? ($xtitle_y += $xpolt_unit) : ($xtitle_x += $xpolt_unit);
}

$tranxy && $svg->line('x1', $flank_x, 'y1', $flank_y, 'x2', $flank_x - +$size_x / 5,
        'y2', $flank_y, 'stroke', 'black', 'stroke-width', 2);
$edge && ($xlab_title = "$x_group[-1]",
    $svg->text('x',$xtitle_x,'y',$xtitle_y + ($tranxy ? $size_x / 3 : 0),'font-family','Arial','stroke','none',
    'fill','black','-cdata',$xlab_title,'font-size',$size_x,'text-anchor', $texanchor));
#if(!($tranxy && $grid)){
    $svg->line('x1',$flank_x,'y1',$height + $flank_y,'x2',$width + $flank_x,'y2',$height + $flank_y,'stroke','black','stroke-width', 2);    # x axis
    $svg->line('x1',$flank_x,'y1',$flank_y,'x2',$flank_x,'y2',$height + $flank_y + ($tranxy ? 0 : ($size_x / 5)),
            'stroke','black','stroke-width', 2);    # y axis
#}

##=======   Draw  the SIGNS  =======##
##
my ($h_posx,$h_posy) = $h_pos ? (split/,/,$h_pos) : ($flank_x + $width/2, $flank_y - 1.5 * $size_ht);
if($nosign){
    $h_posy = $flank_y - 0.5 * $size_ht;
    goto END;
}
if ($sig_pos){
    $sig_pos =~ s/\s+//g;
    ($sig_x, $sig_y) = split /,/, $sig_pos;
}
if (!$right){
    $size_sig_max = 2 * ($pwidth - $sig_x) / max(@sig_long);
    ($size_sig_max < $size_sig) && ($size_sig = $size_sig_max);
}
($sig_y2 > $flank_y) && ($h_posy = $sig_y2 - 1.5*$size_ht);
my $y_top = $sig_y;
my $a_height = 0.7 * $size_sig;    # size length of squear to sign
my $i = 0;
for (0 .. $#signs){
    $rev_sym && ($_ = $#signs - $_);
    if (!$right && $_ && !($_ % $splitn)) {
        $sig_x = $sig_x + $size_sig * $sig_long[$i] / 2;
        $sig_y = $y_top;
        $i++;
    }
    $sig_y += 1.2 * $size_sig;
    my $signs_title = ($signs[$_] || 'Null');
    my $stroke_color = ($stroke || $colors[$_]);
    my @svg_rect = ('stroke', $stroke_color, 'fill',   $colors[$_]);
    $opacity && (push @svg_rect,('fill-opacity', $opacity));
    $svg->rect('x',$sig_x,'y',$sig_y - $a_height,'width', $a_height,'height', $a_height,@svg_rect);
    $svg->text('x', $sig_x + 1.5 * $a_height, 'y', $sig_y,'stroke', 'none','fill','black',
               'font-family','Arial','-cdata', $signs_title,'font-size', $size_sig);
}

END:{
	if($h_title){
        my $may_size_ht = 2*$width/length($h_title);
        ($may_size_ht < $size_ht) && ($size_ht = $may_size_ht);
        $svg->text('x', $h_posx, 'y', $h_posy,'stroke','none','fill','black',
               '-cdata', $h_title,'font-size', $size_ht,'text-anchor', 'middle','font-family','Arial');
    }
for my $al(@add_line){
    my ($ay,$alen,$up_text,$down_text,$acolor,$bcolor) = split/,/,$al;
    $alen ||= 0;
    $acolor ||= 'black';
    $bcolor ||= $acolor;
    $ay = $flank_y+$height-$ay*$ypolt_unit;
    $svg->line('x1',$flank_x-10,'y1',$ay,'x2',$flank_x+$width+$alen,'y2',$ay,'stroke',$acolor,'stroke-width',1);
    $up_text && $svg->text('x',$flank_x+$width,'y',$ay-$size_scale/6,'stroke','none','fill',$bcolor,
            '-cdata',$up_text,'font-family','Arial', 'font-size',$size_scale,'text-anchor','middle');
    $down_text && $svg->text('x',$flank_x+$width,'y',$ay+$size_scale,'stroke','none','fill',$bcolor,
            '-cdata',$down_text,'font-family','Arial', 'font-size',$size_scale,'text-anchor','middle');
}
  print $svg->xmlify;
}

#************************************************************#
#
#                         FUNCTION                           #
#
#************************************************************#

#************#
# Function 1

sub sum{
    my $sum_present = 0;
    foreach (@_)    {
        $sum_present += $_;
    }
    $sum_present;
}

#************#
# Function 2
# this function usde to find the max value in an Array

sub max{
    my $max_present = $_[0];
    foreach (@_)    {
        $max_present = $_ if ($_ > $max_present);
    }
    $max_present;
}

#************#
# Function 3
# This function is used to caclutate the min value and the number of value that show in y axis.

sub axis{
    #useage: axis_spli(a[,b,c]),a is the max value in the axis scale polt
    #b is the ratio of the max value to the length of the Y axis.
    #c is for precision, often use 2,4,8,16
    die"the ratio must between 0.5 to 0.98, if not you should revise your figure.\n"  if ($_[1] > 0.98 || $_[1] < 0.5);
    my ($maxv,$rate,$preci,$sr) = @_;
    $sr && ($maxv /= $sr);
    $rate ||= $scalrate;
    $preci ||= 2;
    $preci = 2**$preci;
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
    my $value_number = int($maxv / $min_value);  # the number of value show in y axis
    ($maxv > $value_number * $min_value) && ($value_number++);
    ($min_value, $value_number);
}

#************#
# Function 4
# This function is used to caclutate the length of the longtest string in an array

sub longest_str{
    my ($group,$cluster,$group_len,$be,$cluhor) = @_;
    my $sig_long = 0;
    my $i = 0;
    my %hash;
    my $clu_len = 0;
    foreach (@{$group}){
        if($be){
            if(/(.+?):/ &&  ! defined $hash{$1}){
                my $pre = $1;
                $cluhor && (length($pre) > $clu_len) && ($clu_len = length($pre));
                $hash{$pre} = 1;
                if(@{$cluster}){
                    push @{$cluster->[-1]},$i-1;
                }
                push @{$cluster},[$pre,$i];
            }
            s/^.+?:\s*//;
        }
        $i++;
        my $long = length($_);
        $be && (push @{$group_len},$long);
        ($long > $sig_long) && ($sig_long = $long);
    }
    if($be && @{$cluster}){
        push @{$cluster->[-1]},$i-1;
        ($sig_long,$clu_len);
    }else{
        $sig_long;
    }
}
#sub5
#===============
sub draw_y_axis{
#===============
    my ($flank_x,$flank_y,$width,$height,$size_scale,$yvalue_min,$yvalue_number,$grid,$micro_scale,$tranxy,$log,$vice) = @_;
    my $yalia_unit = ($tranxy ? $width : $height) / $yvalue_number;
    my $ypolt_unit = $yalia_unit / $yvalue_min;    # The rate of distance in the y ordinate to the real while drawing
    my ($ylab_x, $ylab_y, $ylab, $ss);
    my @group = ('stroke','none','fill','black','font-family','Arial','font-size',$size_scale);
    if ($tranxy){
        $ss = $vice ? -$size_scale/3 : 1.3*$size_scale;
        my $flank_y2 = $vice ? $flank_y : $flank_y + $height;
        push @group,('text-anchor', 'middle');
        for (0 .. $yvalue_number) {
            ($ylab_x, $ylab_y, $ylab) = ($flank_x + $yalia_unit * $_,$flank_y2 + $ss,$_ * $yvalue_min);
#           ($grid && $_ == $yvalue_number) && next;
            $log && ($ylab = $log ** $ylab);
            $svg->text('x',$ylab_x, 'y',$ylab_y,'-cdata',$ylab,@group);
            $svg->line('x1',$ylab_x,'y1',$flank_y2 + $ss / 3,
                   'x2',$ylab_x,'y2',$flank_y2,'stroke','black','stroke-width', $_ ? 1 : 2);
            if($micro_scale && $_ != $yvalue_number){
                foreach(1..4){
                    $ylab_x += $yalia_unit / 5;
                    $svg->line('x1',$ylab_x,'y1',$flank_y2 + $ss / 6,
                   'x2',$ylab_x,'y2',$flank_y2,'stroke','black','stroke-width', $_ ? 1 : 2);
                }
            }
            (!$_ || $_ == $yvalue_number) && next;
            $grid && $svg->line('x1',$ylab_x,'y1',$flank_y,'x2',$ylab_x,'y2',$flank_y + $height,
                'stroke','black','stroke-width', 0.5, 'stroke-dasharray',"3 2");
        }
        $svg->line('x1',$flank_x,'y1',$flank_y2,'x2',$flank_x+$width,'y2',$flank_y2,'stroke','black','stroke-width',1);
    }else{
        $ss = $vice ? -$size_scale : $size_scale;
        my $flank_x2 = $vice ? $flank_x + $width : $flank_x;
        $vice || (push @group,('text-anchor','end'));
        for (0 .. $yvalue_number){
            ($ylab_x, $ylab_y, $ylab) = ($flank_x2 - $ss / 2,$flank_y + $height - $yalia_unit * $_,$_ * $yvalue_min);
            $log && ($ylab = $log ** $ylab);
            $svg->text('x',$ylab_x,'y',$ylab_y + $size_scale / 3,'-cdata',$ylab,@group);
            $svg->line('x1', $flank_x2, 'y1', $ylab_y, 'x2',$flank_x2 - $ss / 3,'y2', $ylab_y, 'stroke', 'black', 'stroke-width', 2);
            if($micro_scale && $_ != $yvalue_number){
                foreach(1..4){
                    $ylab_y -= $yalia_unit/5;
                    $svg->line('x1', $flank_x2, 'y1', $ylab_y, 'x2',$flank_x2 - $ss / 6,'y2', $ylab_y, 'stroke', 'black', 'stroke-width', 2);
                }
            }
            (!$_ || $_ == $yvalue_number) && next;
            $grid && $svg->line('x1', $flank_x, 'y1', $ylab_y, 'x2',$flank_x+$width,'y2', $ylab_y,
                'stroke', 'black', 'stroke-width', 0.5, 'stroke-dasharray',"3 2");
        }
        $svg->line('x1',$flank_x2,'y1',$flank_y,'x2',$flank_x2,'y2',$flank_y+$height,'stroke','black','stroke-width',1);
    }
    $ypolt_unit;
}
#sub5.1
sub caculate_ymax{
    my ($yvalue_number,$yvalue_min,$log) = @_;
    my $y_max = 0;
    for (0..$yvalue_number){
        my $ylab = $_ * $yvalue_min;
        $log && ($ylab = $log ** $ylab);
        ($y_max < length($ylab)) && ($y_max = length($ylab));
    }
    $y_max+1;
}
#sub6 arrow
#==========
sub arrow
#==========
{
    my ($svg,$x,$y,$w,$h,$turn,$wl,$color) = @_;
    $turn ||= 0;#0-down, 1-up, 2-left, 3-right
    $wl ||= 1.5;
    $color ||= 'black';
    my @xy;
    if($turn==0){
        @xy = ($x+$w/2,$y,$x+$w/2,$y+$h);
        push @xy, ($h>$w ? ($x,$y+$h-$w/2,$x+$w,$y+$h-$w/2) : ($x+$w/2-$h/2,$y,$x+$w/+$h/2,$y));
    }elsif($turn==1){
        @xy = ($x+$w/2,$y+$h,$x+$w/2,$y);
        push @xy, ($h>$w ? ($x,$y-$w/2,$x+$w,$y-$w/2) : ($x+$w/2-$h/2,$y+$h,$x+$w+$h/2,$y+$h));
    }elsif($turn==2){
        @xy = ($x+$w,$y+$h/2,$x,$y+$h/2);
        push @xy, ($w > $h ? ($x+$h/2,$y,$x+$h/2,$y+$h) : ($x+$w,$y+$h/2-$w/2,$w+$2,$y+$h/2+$w/2));
    }else{
        @xy = ($x,$y+$h/2,$x+$w,$y+$h/2);
        push @xy, ($w > $h ? ($x+$w-$h/2,$y,$x+$w-$h/2,$y+$h) : ($x,$y+$h/2-$w/2,$x,$y+$h/2+$w/2));
    }
    $svg->line('x1',$xy[0],'y1',$xy[1],'x2',$xy[2],'y2',$xy[3],'stroke',$color,'stroke-width', $wl);
    $svg->line('x1',$xy[4],'y1',$xy[5],'x2',$xy[2],'y2',$xy[3],'stroke',$color,'stroke-width', $wl);
    $svg->line('x1',$xy[6],'y1',$xy[7],'x2',$xy[2],'y2',$xy[3],'stroke',$color,'stroke-width', $wl);
}
