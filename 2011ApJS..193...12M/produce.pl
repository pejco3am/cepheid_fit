#!/usr/bin/perl

open (DOV, "catalog.dat");

while ($radek = <DOV>) {
	@casti = split(/\s+/, $radek);
	$ret = $casti[1];
	$ret =~ s/\_//;
	$name[$casti[0]] = $ret;
}



close(DOV);


open (DOV, "table3.dat");

while ($radek = <DOV>) {
	$name = substr($radek, 0, 9);
	$name =~ s/\s//g;
	$name =~ s/V0/V/;


	$matched = 0;
	for ($i=0;$i<@name;$i++) {
		if (lc($name) eq lc($name[$i])) {
			#print $name, "  ", $i, "\n";
			$match_id = $i;
			$matched = 1;
		}
	}

#	if ($matched == 0) {print $name, "\n";}
	
	$jd = substr($radek, 12, 9);

	$j = substr($radek, 27,6);
	$j_err = substr($radek, 34,5);
	$h = substr($radek, 40,6);
	$h_err = substr($radek, 47,5);
	$k = substr($radek, 53,6);
	$k_err = substr($radek, 60,5);


	print "6 1001 $jd  $j $j_err $match_id $name[$match_id]\n";
	print "7 1001 $jd  $h $h_err $match_id $name[$match_id]\n";
	print "8 1001 $jd  $k $k_err $match_id $name[$match_id]\n";




}

