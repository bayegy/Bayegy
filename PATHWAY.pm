package  PATHWAY;
use strict qw(subs refs);
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw( get_pathway );

#Usage:
#use lib "/System/Pipline/DNA/DNA_Micro/16S_pipeline/16S_pipeline_V1.10/lib/00.Commbin/";
#use PATHWAY;
#
#Example1:
#my ($qimme,$mother,$usearch) = get_path("config.txt",[qimme mother usearch],$Bin,$lib);
#Example2:
#my %path = get_path("config.txt",0,$Bin);
#my ($qimme,$mother,$usearch) = ($path{qimme},$path{mother},$path{usearch});
sub get_pathway{
    my ($config,$outpath,$Bin,$Lib) = @_;
    ($config && -s $config) || return(0);
    my %DIR;
    $Bin && ($DIR{BIN} = $Bin);
    $Lib && ($DIR{LIB} = $Lib);
    my %path;
    open CONFIGPATH,$config || die$!;
    while(<CONFIGPATH>){
        /^#/ && next;
        /\S/ || next;
        chomp;
        if(/(\S+)\s+==\s+(.+)/){
            $DIR{$1} = $2;
        }elsif(/(\S+)\s+=\s+(.+)/){
            my ($key,$value) = ($1,$2);
            if($value !~ /^\//){
                my @val = split/\s+/,$value;
                for(@val){
                    (m#/# && m#^([^/]+)/# && $DIR{$1}) &&
                        (s#^([^/]+)/#$DIR{$1}/#);
                }
                $value = "@val";
            }
            $path{$key} = $value;
        }
    }
    close CONFIGPATH;
    if($outpath){
        my @out;
        for (@$outpath){
            if($path{$_}){
                push @out,$path{$_};
            }else{
                die"error: can't find key work: $_, at $config\n";
            }
        }
        return(@out);
    }else{
        return(%path);
    }
}

1;

__END__
