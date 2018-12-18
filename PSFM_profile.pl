# This file is used to generate sequence profile from multiple sequence alignments (MSA)

@AA=("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");
@MSA=();
undef %ResidueMatrix;
undef %timeMatrix; 
undef %finalFreq;
open(IN,"MSA.txt")||print "can not open MSA.txt";
	while($line=<IN>){
		chomp($line);
		@wds=split(/\s+/,$line);
		push(@MSA,$wds[1]);
	}
close IN;

$sequence=$MSA[0];
$sequence=~s/\-//g;

$query_sequence=$MSA[0];
$thisL = length($query_sequence);
$template_position=0;
for($fnum=1;$fnum<@MSA;$fnum++){
		$template_position++;
		$template_sequence=$MSA[$fnum];
		$pos=0;
		$thisL=length($template_sequence);
		for($i=1;$i<=$thisL;$i++)
		{
        		$R1=substr($query_sequence,$i-1,1);
        		$R2=substr($template_sequence,$i-1,1);
        		$pos++ if($R1 ne '-');
        		if($R1 ne '-' && $R2 ne '-')
						{
	   					$ResidueMatrix{$template_position,$pos}=$R2;
        		}
     		}
}
$template_position++;
for($i=1;$i<=length($sequence);$i++)
{
		$R=substr($sequence,$i-1,1);
    $ResidueMatrix{$template_position,$i}=$R;
}
for($i=1;$i<=length($sequence);$i++){
    for($j=1;$j<=$template_position;$j++){
			$R=$ResidueMatrix{$j,$i};
			$timeMatrix{$R,$i}++;
    }
}

@w=();
for($i=1;$i<=$template_position;$i++){
		$w[$i]=0;
    for($j=1;$j<=length($sequence);$j++){
				$r=0;
				for($k=0;$k<@AA;$k++){
						$A=$AA[$k];
						$value=$timeMatrix{$A,$j};
				    if($value>0){
				    	$r++; 
				    }
				}
				$A=$ResidueMatrix{$i,$j};
				$s=$timeMatrix{$A,$j};
				$w1=1.0/($r*$s);  #1/(r*s) where r = the number of different residues in the position and s = the number of times the particular residue appears in the position
				$w[$i]+=$w1;
    }
    $wSum+=$w[$i];
}
for($i=1;$i<=$template_position;$i++){
    $w[$i]/=$wSum;
}
for($i=1;$i<=$template_position;$i++){
    for($j=1;$j<=length($sequence);$j++){
			$A=$ResidueMatrix{$i,$j};
			$finalFreq{$j,$A}+=$w[$i];
    }
}
open(OUT,">protein.PSFMV2.txt");
#printf OUT "       A          R         N        D          C       Q       E         G         H        I        L        K         M        F        P         S       T       W         Y       V\n";
printf OUT "    %9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s\n","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V";
for($i=1;$i<=length($sequence);$i++){
		$R=substr($sequence,$i-1,1);
    printf OUT "%3d %3s",$i,$R;
    $factory=0;
    for($k=0;$k<@AA;$k++){
    		$R=$AA[$k];
        $factory+=$finalFreq{$i,$R};
    }
    for($k=0;$k<@AA;$k++){
    		$R=$AA[$k];
        printf OUT "%9.6f",$finalFreq{$i,$R}/$factory;
    }
    print OUT "\n";
}
close(OUT);
