#/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Cwd 'abs_path';
use File::Basename;

my $string1= <<'END_MESSAGE1';
DATASET_MULTIBAR
#In multi-value bar charts	 each ID is associated to multiple numeric values	 which are displayed as a stacked or aligned bar chart
#lines starting with a hash are comments and ignored during parsing
#select the separator which is used to delimit the data below (TAB	SPACE or COMMA).This separator must be used throughout this file (except in the SEPARATOR line	 which uses space).

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
SEPARATOR TAB
#SEPARATOR SPACE
#SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL	example multi bar chart

#dataset color (can be changed later)
COLOR	#ff0000

#define colors for each individual field column (if

#field labels

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#dataset scale: you can simply set the values where the scale will be drawn
#DATASET_SCALE	2000	10000	20000
#or you can specify value	 label and color for each scale line (dash separated	 format: VALUE-LABEL-COLOR) 
#DATASET_SCALE	2000-2k line-#0000ff	10000-line at 10k-#ff0000	20000-3rd line-#00ff00

#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#maximum width
WIDTH	1000

#left margin	 used to increase/decrease the spacing to the next dataset. Can be negative	 causing datasets to overlap.
MARGIN	0

#always show internal values; if set	 values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
SHOW_INTERNAL	0

#bar height factor; Default bar height will be slightly less than the available space between leaves	 but you can set a multiplication factor here to increase/decrease it (values from 0 to 1 will decrease it	 values above 1 will increase it)
HEIGHT_FACTOR	1

#Bars are aligned to the node lines by default. Using BAR_SHIFT	 you can move them all up/down by
BAR_SHIFT	0

#align individual fields; if set to 1	 individual bar charts will not be stacked
ALIGN_FIELDS	0
END_MESSAGE1

#FIELD_LABELS	rl0	rl1	rl2	rl3	rl4	rl5	rl6	rl7	rl8
#FIELD_COLORS	#2a9087	#5c2936	#913e40	#2366a1	#658238	#a4f4e0	#a4fc16	#34ec51	#b778a5
my $string4 = <<END_MESSAGE4;
DATASET_SCALE	100	300	500	700

#Internal tree nodes can be specified using IDs directly	 or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA
END_MESSAGE4

my $string2 = <<'END_MESSAGE2';
LABELS
#use this template to change the leaf labels, or define/change the internal node names (displayed in mouseover popups)

#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file (except in the SEPARATOR line, which uses space).

SEPARATOR TAB
#SEPARATOR SPACE
#SEPARATOR COMMA

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA
#NODE_ID,LABEL

#Examples
END_MESSAGE2

my $string3 = <<'END_MESSAGE3';
TREE_COLORS
SEPARATOR TAB
DATA
END_MESSAGE3

my $filename = shift;
my $abs_path_filename = abs_path($filename);
my $dir_name = dirname ($abs_path_filename);
print $dir_name, "\n";

my $n = 0;
my $multibar_n = 0;
my $multibar_color = '';
my $outputfile1 = $dir_name. "/tol_multibar.txt";
my $outputfile2 = $dir_name. "/tol_labels.txt";
my $outputfile3 = $dir_name. "/tol_ranges.txt";
my @samplenames = ();
my @colors = ('#ba1c7a', '#b4c060', '#667332', '#8700c4', '#6e744c', '#1b7679', '#36c0d5', '#bfc2da', '#1d98b7', '#415a40', '#88f99e', '#63c100', '#ebbb20', '#1cdff2', '#d61e3d', '#490cd2', '#bbdcd9', '#087907', '#a690ec', '#f41b61', '#87d2e1', '#9a4c1d', '#e2c55f', '#7b148a', '#a4f0ec', '#be90cf', '#8fa6fe', '#35d17b', '#e1da74', '#eb30a5', '#3f16ab', '#8684e1', '#8cb5d1', '#3eeabb', '#c2d78a', '#06929d', '#4246b7', '#015cdf', '#ae9948', '#a2c83c', '#9b895b', '#298914', '#837fb0', '#78ba5c', '#a6cb3d', '#59f0c8', '#f982a0', '#de64b7', '#b3e605', '#da4180');
my @multibar_colors = ('#ba1c7a', '#b4c060', '#667332', '#8700c4', '#6e744c', '#1b7679', '#36c0d5', '#bfc2da', '#1d98b7', '#415a40', '#88f99e', '#63c100', '#ebbb20', '#1cdff2', '#d61e3d', '#490cd2', '#bbdcd9', '#087907', '#a690ec', '#f41b61', '#87d2e1', '#9a4c1d', '#e2c55f', '#7b148a', '#a4f0ec', '#be90cf', '#8fa6fe', '#35d17b', '#e1da74', '#eb30a5', '#3f16ab', '#8684e1', '#8cb5d1', '#3eeabb', '#c2d78a', '#06929d', '#4246b7', '#015cdf', '#ae9948', '#a2c83c', '#9b895b', '#298914', '#837fb0', '#78ba5c', '#a6cb3d', '#59f0c8', '#f982a0', '#de64b7', '#b3e605', '#da4180');

#my $key = '';
my (%hash, %hashcolor) = ();


sub MAIN(){
	open FH, "<$filename" or die $!;
	while (<FH>){
		chomp;
		if (/^#OTU/){
			@samplenames = split(/\t/);
		}
		else{
			next if /^#/;
			my @array = split(/\t/);
			#print $array[0], "\t", $array[-1], "\n";
			my $temp = $_;
			#put \Q for the regular expression match for special characters.
			$temp =~ s/$array[0]|\Q$array[-1]//g;
			my $tax_short = get_short_tax($array[-1]);
			#print $array[-1], "\n";
			#print $temp, "\n";
			$hash{$array[0]}{'count'} = $temp;
			$hash{$array[0]}{'tax'} = $tax_short;
			if (defined $hashcolor{$tax_short}){
				next;
			}
			else{
				$hashcolor{$tax_short} = shift @colors;
			}
			#print $hashcolor{$tax_short}, "\n";
			$multibar_n = $#array;
		}	
	}

	close(FH);
}


MAIN();

sub get_short_tax() {
	my ($taxonomy) = @_;
	#print $taxonomy, "\n";
	my @taxarray = split (/;/, $taxonomy);
	my $i = 0;
	for ($i = $#taxarray; $i>0; $i--){
		#print $taxarray[$i], "\n";
		if ($taxarray[$i] =~ /__\w+/){
			last;
		}
		else{
			next;
		}
	}
	my $taxshort = $taxarray[$i];
	#print $taxshort, "\n";
	return($taxshort);
}

#FIELD_LABELS	rl0	rl1	rl2	rl3	rl4	rl5	rl6	rl7	rl8
#FIELD_COLORS	#2a9087	#5c2936	#913e40	#2366a1	#658238	#a4f4e0	#a4fc16	#34ec51	#b778a5

open OFH1, ">$outputfile1" or die $!;
print OFH1 $string1;
print OFH1 'FIELD_LABELS';
for (my $i = 1; $i<$multibar_n; $i++){
	#print OFH1 "\t", 'rl', $i;
	print OFH1 "\t", $samplenames[$i];
}
print OFH1 "\nFIELD_COLORS";
for (my $i = 0; $i<$multibar_n - 1; $i++){
	$multibar_color = pop @multibar_colors;
	print OFH1 "\t", $multibar_color;
}
print OFH1 "\n";
print OFH1 $string4;
foreach my $key (sort keys %hash){
	#print $key, $hash{$key}{'count'}, "\n";
	print OFH1 $key, $hash{$key}{'count'}, "\n";
}
close(OFH1);

open OFH2, ">$outputfile2" or die $!;
print OFH2 $string2;
foreach my $key (sort keys %hash){
	print OFH2 $key, "\t", $hash{$key}{'tax'}, "\n";
}
close(OFH2);

#print Dumper @colors;
open OFH3, ">$outputfile3" or die $!;
print OFH3 $string3;
foreach my $key (sort keys %hash){
	print OFH3 $key, "\t", 'range', "\t", $hashcolor{$hash{$key}{'tax'}}, "\t", 'Range', $n, "\n";
	$n ++;
}
close(OFH3);

__DATA__
# Constructed from biom file
#OTU ID	A1	A2	A3	CE1	CE2	CE3	CG1	CG2	CG3	DB1	DB2	DB3	FS1	FS2	FS3	LB1	LB2	LB3	taxonomy
9c10fe267094f807cc5163257e332048	146.0	150.0	7068.0	17.0	837.0	266.0	205.0	81.0	135.0	100.0	23.0	334.0	305.0	32.0	178.0	55.0	248.0	104.0	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Photobacterium; s__damselae
6fc1e9c5930f11516a4fe2e986779b9f	107.0	13.0	0.0	10.0	14.0	36.0	19.0	44.0	22.0	32.0	20.0	21.0	29.0	12.0	19.0	237.0	350.0	52.0	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Alteromonadales; f__Shewanellaceae; g__Shewanella
d43f35dadba69469b5a5648d2a61c1c8	88.0	451.0	573.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Alteromonadales
4f16fba9dad1680e78d8a8ce769da6d7	0.0	0.0	47.0	0.0	37.0	20.0	254.0	0.0	33.0	606.0	260.0	28.0	41.0	15.0	0.0	81.0	66.0	23.0	k__Bacteria; p__Verrucomicrobia; c__Verrucomicrobiae; o__Verrucomicrobiales; f__Verrucomicrobiaceae; g__Luteolibacter; s__
6231baca9ef44eb92fa56d72e35e3d35	12.0	47.0	0.0	0.0	0.0	0.0	42.0	18.0	18.0	167.0	313.0	102.0	85.0	29.0	27.0	306.0	41.0	48.0	k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Tenacibaculum; s__
e7020a746d083ea9b54177e10ec73776	64.0	2453.0	1357.0	0.0	44.0	180.0	748.0	37.0	57.0	67.0	62.0	87.0	336.0	17.0	58.0	128.0	182.0	137.0	k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__; g__; s__
