#Describe: pm for drawing 3D polygon
#Author: Wenbin Liu, liuwenbin@genomics.org.cn
package polygon3D_2;
use Math::Trig;
use strict qw(subs refs);
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(poly_3D);
1;
#sub1
############
sub poly_3D
############
{
        #usage: poly_3D(\@array,line_color,face_color,height,line_width,center_x,center_y,opacity,show_line_to_center,not_show_opacity_view)  -> used to be used, now use hash
        my $svg = shift;
        my %hash = @_;
        my ($a,$color1,$color2,$color3,$h,$lw,$llw,$cx,$cy,$opa,$show_line,$nop);
        if($hash{'points'}){
                $a = $hash{'points'};
        }elsif($hash{'x'} && $hash{'y'} && $hash{'height'} && $hash{'width'}){
                my @x=($hash{'x'},$hash{'x'}+$hash{'width'},$hash{'x'}+$hash{'width'},$hash{'x'});
                my @y=($hash{'y'},$hash{'y'},$hash{'y'}+$hash{'height'},$hash{'y'}+$hash{'height'});
                my @xy;
                foreach(0..3){push @xy,($x[$_],$y[$_]);}
                $a=\@xy;
        }
        $color1 = ($hash{'line-color'} || 'black');
        $color2 = ($hash{'face-color'} || 'none');
        $color3 = ($hash{'focul-color'} || 'black');
        $h = ($hash{'thick'} || 10);
        $lw = ($hash{'line-width'} || 2);
        $llw = ($hash{'focul-width'} || $lw/5);
        $cx = ($hash{'cx'} || -1000);
        $cy = ($hash{'cy'} || -1000);
        $opa = ($hash{'opacity'} || 0.3);
        $show_line = ($hash{'show-focul'} || 0);
        $nop = ($hash{'not-opacity'} || 0);
        my $sel = $hash{'sel_face'} || face_sel($a,$cx,$cy);#sub1.1
        my $fface = $sel->[int($#$sel/2)] -1;
        my $b = low_face($a,$h,$cx,$cy,$fface,$show_line,$color3,$llw);#sub1.2
        draw_face($svg,$a,$b,$color1,$color2,$opa,$lw,$sel,$nop);#sub1.3
}
#sub1.1
############
sub face_sel
############
{
        my ($a,$cx,$cy) = @_;
        my @xy = @{$a};
        my $n = @xy/2;
        my @inter = (0,0,0,0);
        foreach(0..$n-1){
                ($xy[2*$_] >= $cx) && ($inter[0]=1);
                ($xy[2*$_] <= $cx) && ($inter[1]=1);
                ($xy[2*$_+1] >= $cy) && ($inter[2]=1);
                ($xy[2*$_+1] <= $cy) && ($inter[3]=1);
        }
        (sum(@inter)==4) && return([0],0);#sub1.1.1
        my @k;
        foreach(0..$n-1){
                my ($y,$x) = ($xy[2*$_+1]-$cy,$xy[2*$_]-$cx);
                my $r = $y/sqrt($y**2+$x**2);
                $k[$_]= asin($r);
                ($x<0) && ($k[$_] = &pi-$k[$_]);
        }
        my ($maxk,$mink,$max,$min)= edge_face(@k);
        if($maxk - $mink > &pi){
                foreach(@k){($_ < 0) && ($_ += 2*&pi);}
                ($maxk,$mink,$max,$min)= edge_face(@k);#sub1.1.1
        }
        my @sel;
        if($max<$min){
            @sel = ($max..$min-1,0);
        }elsif($min>1){
            @sel = ($max..$n,1..$min-1,0);
        }else{
             @sel = ($max..$n,0);
        }
        \@sel;
}
#sub1.1.1
########
sub sum
########
{
        my $sum=0;
        foreach(@_){$sum+=$_;}
        $sum;
}
#sub1.1.2
#############
sub edge_face
#############
{
        my @k=@_;
        my $n=$#k;
        my ($maxk,$mink,$max,$min)=($k[0],$k[0],1,1);
        foreach(1..$n){
                ($k[$_] > $maxk) && ($maxk=$k[$_],$max=$_+1);
                ($k[$_] > $mink) || ($mink=$k[$_],$min=$_+1);
        }
        ($maxk,$mink,$max,$min);
}
#sub1.2
#############
sub low_face
#############
{
        my ($a,$h,$cx,$cy,$ff,$sl,$cl,$lw) = @_;
        my @xy = @{$a};
        @xy[0,1,2*$ff,2*$ff+1] = @xy[2*$ff,2*$ff+1,0,1];
        my @xy2;
        my ($k,$i) = (0,0);
        foreach(0..@xy/2-1){
                push @xy2,line_pos(@xy[2*$_,2*$_+1],$cx,$cy,$h,$k,$i,$sl,$cl,$lw);#sub1.2.1
        }
        @xy2[2*$ff,2*$ff+1,0,1] = @xy2[0,1,2*$ff,2*$ff+1];
        \@xy2;
}
#sub1.2.1
#############
sub line_pos
#############
{
        my ($x1,$y1,$x2,$y2,$h,$k,$i,$sl,$cl,$lw) = @_;
        $sl && $svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke',$cl,'stroke-width',$lw);
        #($y1==$y2) && die"Error while useing sub poly_3D, input points y can't equal to cy\n";
        if($y1==$y2 && $x1==$x2){
                return($x1,$y1);
        }elsif($y1==$y2){
                $i || ($k=$h/($x2-$x1),(($k<0) && ($k=-$k)),$_[5]=$k,$_[6]=1);
                my $t = ($x2>$x1) ? $x1+$h :  $x1-$h;
                return($t,$y1);
        }else{
                $i || ($k=$h/($y2-$y1),(($k<0) && ($k=-$k)),$_[5]=$k,$_[6]=1);
                return(($x2-$x1)*$k+$x1,($y2-$y1)*$k+$y1);
        }
}

#sub1.3
##############
sub draw_face
##############
{
        my @xy;
        my ($svg,$a,$b,$color1,$color2,$opa,$lw,$sl,$nop) = @_;
        my @xy1 = @{$a};
        my @xy2 = @{$b};
        my $llw = $lw/4;
        my $n = @xy1/2-1;
        foreach(0..$n){
                my @sel = (2*$_..2*$_+3);
                ($_ == $n) && (@sel = (2*$n,2*$n+1,0,1));
                $xy[$_] = [@xy1[@sel],@xy2[@sel[2,3,0,1]]];
                $nop || $svg->polygon('points',[@{$xy[$_]}],'fill','none','stroke',$color1,'stroke-width',$llw);
        }
        unshift @xy,$a;
        foreach(@xy[@{$sl}]){#sub1.2.1
                $svg->polygon('points',[@{$_}],'fill-opacity',$opa,'stroke-opacity',1,'fill',$color2,'stroke',$color1,'stroke-width',$lw);
        }
}

__END__
