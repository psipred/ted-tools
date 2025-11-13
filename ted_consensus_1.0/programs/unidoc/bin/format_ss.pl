#!/usr/bin/perl
use Math::Trig;

my $ss2=$ARGV[0];
my $horiz=$ARGV[1];

my @rst=`cat $ss2`; chomp(@rst);

my $Lch=0;
my %pred=();
foreach my $r(@rst)
{
    next if($r =~ /^#/);
    next if($r !~ /\d+/);

    if($r =~ /(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
    {
    $Lch++;
    #print "$1\t$2\t$3\t$4\t$5\t$6\n";
    $pred{$1, 1}=$2;
    $pred{$1, 2}=$3;
    $pred{$1, 3}=$4;
    $pred{$1, 4}=$5;
    $pred{$1, 5}=$6;
    }
}

my $i=0;
my %conf=();
my @rst=`cat $horiz`;chomp(@rst);
foreach my $r(@rst)
{
    if ($r=~ /Conf:\s+(\S+)/)
    {
    my @a=split(//, $1);
    foreach my $b(@a)
    {
        $i++;
        $conf{$i}=$b;
    }
    }
}

my %secondary=(
            '1'=>'C',
            '2'=>'H',
            '4'=>'E',
            'C'=>'1',
            'H'=>'2',
            'E'=>'4',
    );

my %ts=(
     'GLY'=>'G',
     'ALA'=>'A',
     'VAL'=>'V',
     'LEU'=>'L',
     'ILE'=>'I',
     'SER'=>'S',
     'THR'=>'T',
     'CYS'=>'C',
     'MET'=>'M',
     'PRO'=>'P',
     'ASP'=>'D',
     'ASN'=>'N',
     'GLU'=>'E',
     'GLN'=>'Q',
     'LYS'=>'K',
     'ARG'=>'R',
     'HIS'=>'H',
     'PHE'=>'F',
     'TYR'=>'Y',
     'TRP'=>'W',

     'ASX'=>'B',
     'GLX'=>'Z',
     'UNK'=>'X',

     'G'=>'GLY',
     'A'=>'ALA',
     'V'=>'VAL',
     'L'=>'LEU',
     'I'=>'ILE',
     'S'=>'SER',
     'T'=>'THR',
     'C'=>'CYS',
     'M'=>'MET',
     'P'=>'PRO',
     'D'=>'ASP',
     'N'=>'ASN',
     'E'=>'GLU',
     'Q'=>'GLN',
     'K'=>'LYS',
     'R'=>'ARG',
     'H'=>'HIS',
     'F'=>'PHE',
     'Y'=>'TYR',
     'W'=>'TRP',

     'a'=>'CYS',
     'b'=>'CYS',
     'c'=>'CYS',
     'd'=>'CYS',
     'e'=>'CYS',
     'f'=>'CYS',
     'g'=>'CYS',
     'h'=>'CYS',
     'i'=>'CYS',
     'j'=>'CYS',
     'k'=>'CYS',
     'l'=>'CYS',
     'm'=>'CYS',
     'n'=>'CYS',
     'o'=>'CYS',
     'p'=>'CYS',
     'q'=>'CYS',
     'r'=>'CYS',
     's'=>'CYS',
     't'=>'CYS',
     'u'=>'CYS',
     'v'=>'CYS',
     'w'=>'CYS',
     'x'=>'CYS',
     'y'=>'CYS',
     'z'=>'CYS',

     'B'=>'ASX',
     'Z'=>'GLX',
     'X'=>'CYS',
    );




open(OUT1, ">seq.dat");
for(my $i=1;$i<=$Lch;$i++)
{
    printf OUT1 "%5d   %3s%5d%5d\n", $i, $ts{$pred{$i,1}}, $secondary{$pred{$i, 2}},$conf{$i};
    # printf OUT2 "%4d %1s %1s %6.3f %6.3f %6.3f \n", $i, $pred{$i, 1}, $pred{$i, 2}, $pred{$i, 3},$pred{$i, 4},$pred{$i, 5};
}


close(OUT1);
# close(OUT2);
