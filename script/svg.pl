#!/usr/bin/perl
use strict;
use List::Util qw[min max];
if($#ARGV<1){
    print "$0 infile outfile\n";
    exit(0);
}
my $colWidth=4;
my $row=0;
my $col=0;
my $scale=200;
open OUT,">" ,$ARGV[1];
print OUT "<svg>\n";
open IN,"<",$ARGV[0];
my $pathid=1;
my $nLines = <IN>;
for (my $ii=0;$ii<$nLines;$ii++){
    my $nseg=<IN>;
    for(my $jj=0;$jj<$nseg;$jj++){
        my $npt = <IN>;
        my $line = <IN>;
        my @toks = split /\s+/, $line;
        print OUT '<path stroke-width="1px" fill="none" stroke="#000000"';
        my @mn=(10,10);
        my @mx=(-10,-10);
        for(my $kk=0;$kk<=$#toks;$kk++){
            my @xy=split /,/,$toks[$kk];
            for (my $axis=0;$axis<2;$axis++){
                $mn[$axis] = min($xy[$axis],$mn[$axis]);
                $mx[$axis] = max($xy[$axis],$mx[$axis]);
            }
        }
        my @translate = (0,0);
        for (my $axis=0;$axis<2;$axis++){
            $translate[$axis]=-$mn[$axis];#-($mx[$axis]+$mn[$axis])/2;
        }
        my @spread = (0,0);
        $spread[1]=$row*$scale;
        $spread[0]=$col*$scale;
        $col++;
        if($col==$colWidth){
            $col=0;
            $row++;
        }

        print OUT "\nd=\"M ";
        for(my $kk=0;$kk<=$#toks;$kk++){
            my @xy=split /,/,$toks[$kk];
            for (my $axis=0;$axis<2;$axis++){
                $xy[$axis]+=$translate[$axis];
                $xy[$axis]*=$scale;
                $xy[$axis]+=$spread[$axis];
                
              
                print OUT $xy[$axis];
                if($axis==0){
                    print OUT ",";
                }else{
                    print OUT " ";
                }
            }        
        }
        print OUT "z\"\n";
        print OUT "id=\"path$pathid\"\n/>\n";
        $pathid++;
    }
    <IN>;
}
close IN;
print OUT "</svg>\n";
close OUT;
