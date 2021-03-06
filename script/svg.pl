#!/usr/bin/perl
use strict;
use List::Util qw[min max];
if($#ARGV<1){
    print "$0 infile outfile\n";
    exit(0);
}
my $colWidth=5;
my $row=0;
my $col=0;
my $scale=1000;
my $spreadscale=400;
my $width=1800;
open OUT,">" ,$ARGV[1];
print OUT "<svg>\n";
open IN,"<",$ARGV[0];
my $pathid=1;
my $nLines = <IN>;
my $left=0;
my $top=0;
 my @spread = (0,0);
$spread[1]=$top;
       
for (my $ii=0;$ii<$nLines;$ii++){
    my $nseg=<IN>;
    if($nseg>0){
        print OUT "<g>\n";      
        <IN>;#normal
    }#else{
      #  next;
   # }
    for(my $jj=0;$jj<$nseg;$jj++){
        my $npt = <IN>;
        my $line = <IN>;
        <IN>;#vertex indices
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
        my @center = (0,0);
        print OUT "\nd=\"M ";
        $spread[0]=$left;#$col*$spreadscale;
        for(my $kk=0;$kk<=$#toks;$kk++){
            my @xy=split /,/,$toks[$kk];
            for (my $axis=0;$axis<2;$axis++){
                $xy[$axis]+=$translate[$axis];
                $xy[$axis]*=$scale;
                $xy[$axis]+=$spread[$axis];
               
                $center[$axis]+=$xy[$axis];
              
                print OUT $xy[$axis];
                if($axis==0){
                    print OUT ",";
                }else{
                    print OUT " ";
                }
            }
            if($xy[0]>$left){
                $left=$xy[0];
            }
            if($xy[1]>$top){
                $top=$xy[1];
            }
        }
        $left+=5;
        $col++;

        if($left>=$width){
            $col=0;
            $left=0;
            $row++;
            $spread[1]=$top+5;
        }

        print OUT "z\"\n";
        print OUT "id=\"path$pathid\"\n/>\n";
        for (my $axis=0;$axis<2;$axis++){
            $center[$axis]/=($#toks+1);
        }
        $center[0]-=15;
        print OUT "<text x=\"$center[0]\" y=\"$center[1]\"\n"; 
        print OUT "font-family=\"Arial\" font-size=\"30\" fill=\"blue\" >\n";
        print OUT "$ii\n</text>\n";
        print OUT "</g>\n";
        $pathid++;

    }
    <IN>;
}
close IN;
print OUT "</svg>\n";
close OUT;
