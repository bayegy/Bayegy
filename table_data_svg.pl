#!/use/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/";
use PATHWAY;
(-s "$Bin/Pathway_cfg.txt") || die"error: can't find config at $Bin/, $!\n";
BEGIN {
    my( $svg_lib)= get_pathway("$Bin/../../bin/Pathway_cfg.txt",[qw(SVG_Lib)]);
    unshift @INC, $svg_lib
};
use polygon3D_2;
use SVG;
use Getopt::Long;
my ($flankx,$flanky,$height,$fsize,$range,$symbol,$tsize,$xh) = (100,100,0,15,0,"Identity(%)",10,30);
my ($colors,$simply,$redcut,$fleval,$noline,$D3,$lcolor,$scolor,$least);
GetOptions(
    "flank_x:s" => \$flankx,
    "flank_y:s" => \$flanky,
    "height:i"  => \$height,
    "fsize:f"   => \$fsize,
    "range:s"   => \$range,
    "symbol:s"  => \$symbol,
    "colors:s"  => \$colors,
    "simply"    => \$simply,
    "redcut:f"  => \$redcut,
    "tsize:f"   => \$tsize,
    "fleval"    => \$fleval,
    "noline"    => \$noline,
    "D3"        => \$D3,
    "lcolor:s"  => \$lcolor,
    "scolor:s"  => \$scolor,
    "least:f"   => \$least,
    "xh:f"      => \$xh
);
@ARGV || die"Name: table_data_svg.pl
Usage: perl table_data_svg.pl <in.table> > out.svg
    --flank_x <num>     left,right x border width(e.g. 100,100), default auto
    --flank_y <num>     up,down y border width(e.g. 100,100), default auto
    --height <num>      figure height, default=300
    --fsize <flo>       font size, default=15
    --range <str>       value range(e.g. 0,100), default auto
    --colors <str>      set range colors, default=lightgreen-yello
    --symbol <str>      symbol text, default='Identity(%)'
    --simply            to simply the figure(del the row below)
    --redcut <flo>      cut off for mask the number in red, default not set
    --tsize <flo>       text font-size, default=10
    --fleval            show color fill leval
    --noline            not show the line,
    --lcolor <str>      line color, default=black
    --D3                show 3D view\n\n";
#   --scolor <str>      stroke color at middle, default=none
#   --least <flo>       least fill leval for --fleval, default=0
#   --xh <flo>          height for each sample, default=20
#==========================================================
my (@name1,@name2,@data);
my $num;
my (%nameh1,%nameh2);
for my $f(0..$#ARGV){
    my (@name_1,@name_2,@ranks);
    open IN,$ARGV[$f] || die$!;
    $num = 0;
    while(<IN>){
        if(m/^\s+/ || !@name_1){
            chomp;
            @name_1 = /\t/ ? split/\t/ : split/\s+/;
            shift @name_1;
            if($f){for my $n(@name_1){push @ranks,$nameh1{$n};}}
        }else{
            my @l = split/\s+/;
            push @name_2,(shift @l);
            @ranks && (@l = @l[@ranks]);
            my $k = defined $nameh2{$name_2[-1]} ? $nameh2{$name_2[-1]} : $num;
            for my $i(0..$#l){
                push @{$data[$k]->[$i]}, (split /\//,$l[$i]);
            }
            $num++;
        }
    }
    close IN;
    if(!$f){
        @name1 = @name_1;
        @name2 = @name_2;
        if($#ARGV){
            for my $i(0..$#name1){$nameh1{$name1[$i]} = $i;}
            for my $i(0..$#name2){$nameh2{$name2[$i]} = $i;}
        }
    }
}
my $sn = $simply ? 1 : 0;
if(!$height){
    $height = $xh * @name1;
    ($height < 300) && ($height = 300);
}
#==========================================================
my $width = 2 * $height;
my $unit = $height / ($num - $sn); # del the row bellow
#my $unit = $height/ $height;
my ($flank_x,$flank_x2,$flank_y,$flank_y2);
my @flankxy;
if($flankx && $flankx=~/\S+,\S+/ && $flanky && $flanky=~/\S+,\S+/){
}else{
    @flankxy = mflank_xy($unit,$fsize,\@name1,\@name2,$sn,5);
    push @flankxy,6*$tsize;
}
$flankx && ($flankx=~/(\S+),(\S+)/) && (@flankxy[0,1] = ($1,$2));
$flanky && ($flanky=~/(\S+),(\S+)/) && (@flankxy[2,3] = ($1,$2));
($flank_x,$flank_y) = @flankxy[0,2];
my $pheight = $flankxy[2] + $flankxy[3] + $height + $unit;
my $pwidth = $flankxy[0] + $flankxy[1] + $width;
my $svg = SVG->new(width=>$pwidth, height=>$pheight);
my $sfsize = $unit / 3;
($sfsize > 0.85*$fsize) && ($sfsize = 0.85*$fsize);
$sfsize = sprintf("%.2f",$sfsize);
$colors ||= "$Bin/colors.txt";
my @color = (-s $colors) ? split/\s+/,`less $colors` : get_color("$Bin/rgb_colors.txt",$colors);
my ($x0,$y0) = ($flank_x+$width/2,$flank_y);
my ($minv,$maxv) = $range ? (split/,/,$range) : min_max(@data);
my ($cx,$cy) = ($flank_x + $width/2, 2*$pheight);
for my $i(0..$num-1-$sn){
    my ($x1,$y1) = ($x0 + $i*$unit, $y0 + $i*$unit);
	for my $j(0..$num-1-$i-$sn){
		map_color($svg,$data[$i]->[$num-1-$j],$x1,$y1,$unit,$sfsize,\@color,$minv,$maxv,$redcut,$fleval,$D3,$cx,$cy,$scolor);
		$x1 -= $unit;
		$y1 += $unit;
	}
}
#==========================================================
my ($x1,$y1) = ($flank_x,$flank_y+$height);
my ($x2,$y2) = ($x1+$unit,$y1+$unit);
my ($x3,$y3) = ($flank_x+$width/2,$flank_y);
($fsize > 1.4*$unit) && ($fsize = 1.4*$unit);
$lcolor ||= 'black';
$noline || $svg->line('x1',$x1,'y1',$y1,'x2',$x3,'y2',$y3,'stroke', $lcolor, 'stroke-width',1);
my ($rx,$ry) = ($x3+$unit/2+$fsize/2, $y3+$unit/2-$fsize/2);
my $group = $svg->group("transform","rotate(-45,$rx,$ry)",'stroke','none',
    'fill','black','font-size',$fsize,'font-family','Arial');
$group->text('x',$rx,'y',$ry+$fsize/3,'-cdata',$name1[0]);
$x3 += $unit;
$y3 += $unit;
foreach(0..$num-$sn){
    if($_ == $num-$sn){
        $x2 -= $unit;
        $y2 -= $unit;
    }else{
        $noline || $svg->line('x1',$x3,'y1',$y3,'x2',$x2,'y2',$y2,'stroke', $lcolor, 'stroke-width',1);
        ($rx,$ry) = ($x3+$unit/2+$fsize/2, $y3+$unit/2-$fsize/2);
        $group = $svg->group("transform","rotate(-45,$rx,$ry)",'stroke','none',
            'fill','black','font-size',$fsize,'font-family','Arial');
        ($sn && $_==$num-$sn-1) || $group->text('x',$rx,'y',$ry+$fsize/3,'-cdata',$name1[$_+1]);
        ($rx,$ry) = ($x1+$unit/2-$fsize/2, $y1-$unit/2-$fsize/2);
        $group = $svg->group("transform","rotate(45,$rx,$ry)",'stroke','none',
            'fill','black','font-size',$fsize,'font-family','Arial');
        $group->text('x',$rx,'y',$ry+$fsize/3,'-cdata',$name2[$_+$sn],'text-anchor','end');
    }
    $noline || $svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke', $lcolor, 'stroke-width',1);
    $x1 += $unit;
    $y1 -= $unit;
    $x2 += 2*$unit;
    $x3 += $unit;
    $y3 += $unit;
}
#==========================================================
$x1 = $flank_x + $unit/2;
$y1 = $flank_y + $height + $unit + $tsize;
$svg->text('x',$x1+$width/6,'y',$y1+4.333*$tsize,'stroke','none','fill','black','-cdata',$symbol,
	'font-size',$tsize,'text-anchor','middle', 'font-family', 'Arial');
$svg->line('x1',$x1,'y1',$y1+2*$tsize,'x2',$x1+$width/3,'y2',$y1+2*$tsize,
'stroke', $lcolor, 'stroke-width',0.5);
my $midv = ($minv+$maxv)/2;
for ([0,$minv],[$width/6,$midv],[$width/3,$maxv]){
	$svg->line('x1',$x1+$_->[0],'y1',$y1+2*$tsize,'x2',$x1+$_->[0],'y2',$y1+2.3*$tsize,
	'stroke', 'black', 'stroke-width',0.5);
	$svg->text('x',$x1+$_->[0],'y',$y1+3.2*$tsize,'stroke','none','fill','black','-cdata',$_->[1],
	'font-size',$tsize,'text-anchor','middle', 'font-family', 'Arial');
}

my $gra1 = $svg->gradient(-type => "linear",id => "gradient_1",
gradientUnits => "userSpaceOnUse",x1=>$x1, y1=> 0, x2 =>$x1+$width/3,y2=>0);
$gra1->stop('offset'=>0,style=>"stop-color:$color[0]",'stop-opacity',1);
$gra1->stop('offset'=>"100%",style=>"stop-color:$color[-1]",'stop-opacity',1);
$svg->rect('x',$x1,'y',$y1,'width',$width/3,'height',2*$tsize,'fill','url(#gradient_1)');#,'stroke',$lcolor);



print $svg->xmlify;
#==========================================================
sub mflank_xy{
    my ($unit,$fsize,$data1,$data2,$sn,$def_edge) = @_;
    my ($xedge2,$yedge2) = flank_xy($unit,$fsize,[@{$data1}[0..$#$data1-$sn]],$def_edge);
    my ($xedge1,$yedge1) = flank_xy($unit,$fsize,[reverse @{$data2}[$sn .. $#$data2]],$def_edge);
    ($yedge1 < $yedge2) && ($yedge1 = $yedge2);
    ($xedge1,$xedge2,$yedge1);
}
sub flank_xy{
    my ($unit,$fsize,$data,$def_edge) = @_;
    $def_edge ||= 5;
    my $yend = $unit * 0.5;
    my $xend = $unit * ($#$data + 0.5);
    my ($xedge,$yedge) = ($def_edge, $def_edge);
    foreach(@$data){
        my $sq2 = $fsize*(sqrt(2)*length($_)/4+2);
        my ($x,$y) = ($sq2-$xend,$sq2-$yend);
        ($x > $xedge) && ($xedge = $x);
        ($y > $yedge) && ($yedge = $y);
        $xend -= $unit;
        $yend += $unit;
    }
    ($xedge, $yedge);
}
sub map_color{
	my ($svg,$data,$x,$y,$unit,$fsize,$color,$minv,$maxv,$redcut,$fleval,$D3,$cx,$cy,$stroke_color) = @_;
    $stroke_color ||= 'none';
    for(@$data){
        m/\.\d\d\d\d/ && ($_ = sprintf("%.3f",$_));
    }
	my @id = @$data;
	my $range = $maxv - $minv;
    my $radio;
    my $thick0 = $unit / 2;
    $cx ||= 1000;
    $cy ||= 1000;
    my $id_num = @id;
    for my $i(0..$id_num-1){
        $radio = ($id[$i]-$minv) / $range;
        my @pos;
        my $tex_yp = out_pos($i,$unit,$x,$y,$fsize,$fleval,$radio,$id_num,$least,\@pos);
    	my $col = $color->[int($#$color*$radio)];
        my $thick = $thick0 * $radio;
        $D3 ? poly_3D($svg,'points',$pos[$i],'face-color',$col,'line-color','green','thick',$thick,'cx',$cx,'cy',$cy,'opacity',0.6,'sel_face',[2,0,3]) :
	    $svg->polygon('points',$pos[$i],'fill', $col, 'stroke',$stroke_color);
        ($id_num==2 && $i==1) && ($id[1] = "($id[1])");
        my $text_col = ($redcut && $id[$i]>=$redcut) ? 'red' : 'black';
	    $svg->text('x',$x, 'y',$tex_yp,'stroke','none','fill',$text_col,'-cdata',$id[$i],
    	'font-size',$fsize,'text-anchor','middle', 'font-family', 'Arial');
    }
}
sub out_pos{
    my ($i,$unit,$x,$y,$fsize,$fleval,$radio,$id,$least,$pos) = @_;
    my $unit0 = $unit;
    my $dy = 0;
    if($fleval){
        $least ||= 0;
        my $rr = 1 - $least;
        my $dy = $unit*(1-$rr*$radio-$least)/$id;
        ($id==3) && ($dy *= $i==1 ? 9/8 : 3/4);
        $unit *= ($rr*$radio+$least);
    }
    my ($unit1,$unit2) = (5*$unit0/4,3*$unit0/4);
#    @{$pos} = ($id == 1) ? ([$x,$y+$dy,$x-$unit,$y+$unit0,$x,$y+2*$unit0-$dy,$x+$unit,$y+$unit0]) :
    @{$pos} = ($id == 1) ? ([$x,$y+$unit0-$unit,$x-$unit,$y+$unit0,$x,$y+$unit0+$unit,$x+$unit,$y+$unit0]) :
        ($id == 2) ? ([$x,$y+$dy,$x-$unit,$y+$unit0-$dy,$x+$unit,$y+$unit0-$dy],[$x,$y+2*$unit0-$dy,$x-$unit,$y+$unit0+$dy,$x+$unit,$y+$unit0+$dy]) :
        ([$x,$y+$dy,$x-$unit2+$dy,$y+$unit2-$dy,$x+$unit2-$dy,$y+$unit2-$dy],
         [$x-$unit2+$dy,$y+$unit2+$dy,$x-$unit0+$dy,$y+$unit0,$x-$unit2+$dy,$y+$unit1-$dy,$x+$unit2-$dy,$y+$unit1-$dy,$x+$unit0-$dy,$y+$unit0,$x+$unit2-$dy,$y+$unit2+$dy],
         [$x-$unit2+$dy,$y+$unit1+$dy,$x,$y+2*$unit0-$dy,$x+$unit2-$dy,$y+$unit1+$dy]);
    $noline || ($dy = 0, $unit=$unit0);
    my @tex_yp = ($id == 1) ? ($y+$unit0+$fsize/3) : ($id == 2) ? ($y+$unit-$dy-$fsize/3,$y+$unit+$dy+$fsize) :
        ($y+$unit2-$dy-$fsize/3,$y+$unit0+$fsize/3,$y+$unit1+$dy+$fsize);
    $tex_yp[$i];
}
sub get_color{
    my ($rgb_svg,$color,$num) = @_;
    my @c = split /-/,$color;
    $c[1] || die"erro form at --colors, $!";
    if($c[0]=~/[^\d,]/ || $c[1]=~/[^\d,]/){
        my %rgbh = split/\s+/,`less $rgb_svg`;
        ($c[0]=~/[^\d,]/) && ($c[0] = $rgbh{$c[0]} || '255,255,255');
        ($c[1]=~/[^\d,]/) && ($c[1] = $rgbh{$c[1]} || '255,255,255');
        %rgbh = ();
    }
    my @r1 = split/,/,$c[0];
    my @r2 = split/,/,$c[1];
    my @dis;
    for (0..2){
        push @dis,($r2[$_]-$r1[$_]);
    }
    if(!$num){
        $num = 0;
        for (@dis){(abs($_) > $num) && ($num = abs($_));}
    }
    for (@dis){$_ /= $num;}
    my @out_color;
    foreach (0..$num){
        push @out_color, ( "rgb(" . join(",",int($r1[0]+0.5),int($r1[1]+0.5),int($r1[2]+0.5)) . ")");
        for my $i(0..2){
           $r1[$i] += $dis[$i];
        }
    }
    push @out_color, ( "rgb(" . join(",",int($r1[0]+0.5),int($r1[1]+0.5),int($r1[2]+0.5)) . ")");
   @out_color;
}
sub min_max{
    my ($min,$max);
    foreach my $k(@_){
        for my $j(@{$k}){
            for my $i(@{$j}){
                (!(defined $min) || $i < $min) && ($min = $i);
                (!(defined $max) || $i > $max) && ($max = $i);
            }
        }
    }
    mul_axis($min,$max);
}

#### from COMM.pm /liuwenbin

#sub1.1.1
#============#
sub mul_axis{
    my ($min_x,$max_x,$prec,$rate) = @_;
    $prec ||= 4;
    $rate ||= 0.95;
    my ($u, $n) = axis_split($max_x-$min_x, $rate, $prec); #sub1.1.1.1
    my $min = $u * int($min_x / $u);
    until($min + $u * $n >= $max_x){$n++;}
    ($min,$min+$u*$n);
}
#sub1.1.1.1
#=============#
sub axis_split{
    my ($maxv,$rate,$preci) = @_;
    $rate ||= 0.9;
    die"the ratio must between 0.5 to 0.98, if not you should revise your figure.\n"  if ($rate > 0.98 || $rate < 0.5);
    $preci ||= 2;
    $preci = 2**$preci;
    sprintf("%1.1e", $maxv) =~ /^(.)(.*)e(.*)/;
    my $mbs = $1;    # the MSB of the max value in the plot
    my $mag = $3;    # the order of magnitude of the max value in the plot.
    $mag =~ s/^\+//;
    $mag =~ s/^0*//;
    $mag  ||= 0;
    my $k = $rate / (1 - $rate) / $preci;              # the middle value used to caclutate $min_value-
                       # -you can also change preci into 2 or 1, the y-scal will become more precision
    my $min_value;     # the min value show in y axis
    foreach(2,1,0.5,0.25,0.125,0.1,0.05){
        $min_value = $_;
        ($mbs >= $_ * $k) && last;
    }
    $min_value = $min_value * 10**$mag;
    my $value_number = int($maxv / $min_value);  # the number of value show in y axis
    ($value_number * $min_value == $maxv) || ($value_number++);
    ($min_value, $value_number);
}

