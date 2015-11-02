#!/usr/bin/perl


open (DOV, "photo.dat");

$cep = "HW Car";

while ($radek = <DOV>) {
	$name = substr($radek, 0, 8);
	$name =~ s/\s+$//;
	if (lc($cep) eq lc($name)) {
		$jd = substr($radek, 18, 10);
		$v  = substr($radek, 29, 6);
		$ub = substr($radek, 36, 6);
		$bv = substr($radek, 43,6);
		$vr = substr($radek, 50,6);
		$ri = substr($radek, 57,6);
		$vi = substr($radek, 64,6);
		
		if ($v =~ /\d+/)  {
			print "3 150 $jd $v -10.000 577 HW_Car\n";
		} else {next;}
	
		if ($bv =~ /\d+/) {
			$b = $bv + $v;
			printf "2 150 $jd %6.3f -10.000 577 HW_Car\n", $b;
		}
		if ($ub =~ /\d+/) {
			$u = $ub + $b;
			printf "1 150 $jd %6.3f -10.000 577 HW_Car\n", $u;
		}

		if ($vr =~ /\d+/) {
			$r = $v - $vr;
			printf "4 150 $jd %6.3f -10.000 577 HW_Car\n", $r;
		}

		if ($ri =~ /\d+/) {
			$i = $r - $ri;
			printf "5 150 $jd %6.3f -10.000 577 HW_Car\n", $i;
		} else {
			if ($vi =~ /\d+/) {
				$i = $v - $vi;
				printf "5 150 $jd %6.3f -10.000 577 HW_Car\n", $i;		
			}
		}


	}


}
