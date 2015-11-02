#!/usr/bin/perl


open(DOV, "ven.dat");

$i = 1;
$first = 1;

while ($radek = <DOV>) {
	$radek =~ s/^\s+//;
	@casti = split(/\s+/, $radek);
	if ($first == 1) {
		$name[$i] = $casti[6];
		print $i, "  ", $name[$i], "\n";
		$i++;
		$first = 0;
	} else {
	if ($casti[6] ne $name[$i-1]) {
		$name[$i] = $casti[6];
		print $i, "  ", $name[$i], "\n";
		$i++;	
	}
	}


}
