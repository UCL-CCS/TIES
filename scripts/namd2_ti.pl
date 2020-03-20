#!/usr/bin/perl -w
use strict;
use Getopt::Long; # for command line options
$| = 1; # regular print buffer flush

my $integralsizelimit = 30; #arbitrary size limit on single trapezoid area
my %warnings; #put stuff in here if consistency checks fail

#command line options
my $tailsample = ''; my $blocksample = ''; my $spline = 0;
my ($unsafe, $keepequilibration, $fulloutput, $dgonly, $help, $oldtrapezoidal);

if (! GetOptions ('blocksample=s' => \$blocksample, 'tailsample=s' => \$tailsample, 'unsafe' => \$unsafe, 'keepequilibration' => \$keepequilibration, 'fulloutput' => \$fulloutput, 'dgonly' => \$dgonly, 'help|?' => \$help, 'oldtrapezoidal' => \$oldtrapezoidal, 'spline=s' => \$spline))
  {print ("Run with -h or --help for instructions\n"); exit();}
if ($help) {helpmsg(); exit(); }

if ($tailsample && (! ($tailsample > 0 && $tailsample <= 100)) || ($tailsample eq '0')) {
  die("tailsample: specify percentage (up to 100)\n");
}

my @outfiles;
if (scalar(@ARGV)) {
  if (! -e $ARGV[0]) {die "Specified location $ARGV[0] not found - run with --help for help message\n";}
  if (scalar(@ARGV) == 1 && -d $ARGV[0]) {
    print "Looking for files *ti.out in directory $ARGV[0]...\n";
    @outfiles = glob("$ARGV[0]*ti.out");
  }
  else {    
    @outfiles = @ARGV;
  }
}
else {
  @outfiles = glob("*ti.out");
}

if (scalar(@outfiles == 0)) {  print "No TI output files found. Run with -h for help message.\n";exit()}
if (!$dgonly) {print ("Found ",scalar(@outfiles)," output files\n");}

my (%elec1, %elec2, %vdw1, %vdw2, %pme); #key = respective lambda value; value = array of dU/dlambda values
my %timesteps; #key = 'global' lambda value; value = number of timesteps (not the number of data points)
my ($lambda, $vdwl1, $vdwl2, $elecl1, $elecl2); # current lambda for each component
my $outputfreq = 0; #frequency of output in timesteps
my $prevstepnr; #previous MD step, to check if anything is missing or inserted
my $datalinecount;
my @outputfreqs; #add in all outputfreqs found, warn if they aren't the same across the board
my $linenr;
my %lambda_hash; #key = global lambda as set in conf file; value = hash of individual vdw and elec lambdas associated with each global lamda.


my $elecl1_zero_highestlambda = 0;
my $elecl2_zero_lowestlambda = 1;
my $vdwl1_one_lowestlambda = 1;
my $vdwl2_one_highestlambda = 0;


if (!$dgonly) {print ("Reading files: ");}

foreach my $file (sort{ -M $b <=> -M $a } @outfiles) {   # cycle through files... oldest ones first to give block averaging / tail sample a chance of making sense over multiple files
  $outputfreq = 0; $prevstepnr = 0; $linenr = 0; #$datalinecount = 0;
  if ($fulloutput) {print "$file...   ";}
  elsif (!$dgonly) {print ".";}
  open (FILE, "$file");
  my $first = <FILE>;
  if (! ($first =~ /\#       STEP      Elec_d[EU]\/dl/) ) {
    print ("Output file $file not recognised as a NAMD TI output file\n"); next;
  }
  while (my $line = <FILE>) {  # cycle through every line...
    $linenr ++;
    if ($line =~ /^\#/) { # starting with #, this is some kind of control information
      if ($line =~ /WINDOW: LAMBDA (\S+)/) { $lambda = $1; $outputfreq = 0; $prevstepnr = 0; $datalinecount = 0;}
      #plugging in some more consistency checks. If in any case, a given global lambda value is associated with more than one elec/vdw lambda value, something has gone wrong and the data is no longer self-consistent...
      if ($line =~ /1 VDW LAMBDA (\S+)/) { 
        $vdwl1 = $1; 
        if (($vdwl1 == 1) && ($lambda < $vdwl1_one_lowestlambda)) {$vdwl1_one_lowestlambda = $lambda; delete($vdw1{$vdwl1}); }
        if (defined(${lambda_hash{$lambda}}{vdw1}) && (${$lambda_hash{$lambda}}{vdw1} != $vdwl1)) {$warnings{lambdaconsistency} = 1;}
        else {${$lambda_hash{$lambda}}{vdw1} = $vdwl1; }
      }
      if ($line =~ /2 VDW LAMBDA (\S+)/) { 
        $vdwl2 = $1; 
        if (($vdwl2 == 1) && ($lambda > $vdwl2_one_highestlambda)) {$vdwl2_one_highestlambda = $lambda; delete($vdw2{$vdwl2}); }
        if (defined(${lambda_hash{$lambda}}{vdw2}) && (${$lambda_hash{$lambda}}{vdw2} != $vdwl2)) {$warnings{lambdaconsistency} = 1;}
        else {${$lambda_hash{$lambda}}{vdw2} = $vdwl2; }
      }
      if ($line =~ /1 ELEC LAMBDA (\S+)/) { 
        $elecl1 = $1; 
        if (($elecl1 == 0) && ($lambda > $elecl1_zero_highestlambda)) {$elecl1_zero_highestlambda = $lambda; delete($elec1{$elecl1}); }
        if (defined(${lambda_hash{$lambda}}{elec1}) && (${$lambda_hash{$lambda}}{elec1} != $elecl1)) {$warnings{lambdaconsistency} = 1;}
        else {${$lambda_hash{$lambda}}{elec1} = $elecl1; }
      }
      if ($line =~ /2 ELEC LAMBDA (\S+)/) { 
        $elecl2 = $1; 
        if (($elecl2 == 0) && ($lambda < $elecl2_zero_lowestlambda)) {$elecl2_zero_lowestlambda = $lambda; delete($elec2{$elecl2}); }
        if (defined(${lambda_hash{$lambda}}{elec2}) && (${$lambda_hash{$lambda}}{elec2} != $elecl2)) {$warnings{lambdaconsistency} = 1;}
        else {${$lambda_hash{$lambda}}{elec2} = $elecl2; }
      }
      if ($line =~ /\#(\d+) STEPS OF EQUILIBRATION.*COMPLETED/) { #go back and throw out equilibration data, unless requested to keep
        if (! $keepequilibration) {
          my $numbertodelete;
          if (!$outputfreq) {$numbertodelete = $datalinecount;} # covers us if there's little or no output between start of run and end of equilibration
          else {$numbertodelete = $1 / $outputfreq; }
          splice(@{$elec1{$elecl1}}, scalar(@{$elec1{$elecl1}}) - $numbertodelete,$numbertodelete); # delete last $numbertodelete values from end of array (x5)
          splice(@{$elec2{$elecl2}}, scalar(@{$elec2{$elecl2}}) - $numbertodelete,$numbertodelete);
          splice(@{$vdw1{$vdwl1}}, scalar(@{$vdw1{$vdwl1}}) - $numbertodelete,$numbertodelete);
          splice(@{$vdw2{$vdwl2}}, scalar(@{$vdw2{$vdwl2}}) - $numbertodelete,$numbertodelete);
          splice(@{$pme{$lambda}}, scalar(@{$pme{$lambda}}) - $numbertodelete,$numbertodelete);
          $timesteps{$lambda} -= $numbertodelete * $outputfreq;
          $warnings{equilibrationdeleted} ++;
        }
        $outputfreq = 0; $prevstepnr = 0;
      }
    }
    if ($line =~ m/TI\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ ) { # match the actual data lines
      $datalinecount ++;
      if (!$outputfreq) { # get the outputfreq as stepnr - prevstepnr, if we don't have it already
        if (!$prevstepnr) {$prevstepnr = $1;} 
        else {$outputfreq = $1 - $prevstepnr; push @outputfreqs, $outputfreq; $prevstepnr = $1;} 
      }
      else { # check if MD step follows the last one
        if (($1 != ($prevstepnr + $outputfreq)) && (!$unsafe)) {die("\nDiscontinuity in data, file $file, line $linenr - probably indicates missing or corrupted data. Run with --unsafe to ignore.\n$line");}
        $prevstepnr = $1;
      }
      # for each of these dE/dls: add the instantaneous dE/dl value onto an array, which sits inside hash with key=lambda, value=dE/dl array
      if ($lambda >= $elecl1_zero_highestlambda) {  push (@{$elec1{$elecl1}}, $2); }
      if ($lambda <= $elecl2_zero_lowestlambda) {  push (@{$elec2{$elecl2}}, $6); }
      if ($lambda <= $vdwl1_one_lowestlambda) {  push (@{$vdw1{$vdwl1}}, $4); }
      if ($lambda >= $vdwl2_one_highestlambda) {  push (@{$vdw2{$vdwl2}}, $8); }
      push (@{$pme{$lambda}}, 0); 
      $timesteps{$lambda} += $outputfreq;
    }
  }
  if ($tailsample) { # if requested, delete a portion of the data, using only the end portion of the sample
    $timesteps{$lambda} -= int( scalar(@{$pme{$lambda}}) * (1 - ($tailsample / 100))) * $outputfreq;
    splice(@{$elec1{$elecl1}},0, int(scalar(@{$elec1{$elecl1}}) * (1 - ($tailsample / 100))));
    splice(@{$elec2{$elecl2}},0, int(scalar(@{$elec2{$elecl2}}) * (1 - ($tailsample / 100))));
    splice(@{$vdw1{$vdwl1}},0, int(scalar(@{$vdw1{$vdwl1}}) * (1 - ($tailsample / 100))));
    splice(@{$vdw2{$vdwl2}},0, int(scalar(@{$vdw2{$vdwl2}}) * (1 - ($tailsample / 100))));
    splice(@{$pme{$lambda}},0, int( scalar(@{$pme{$lambda}}) * (1 - ($tailsample / 100))));
  }
  close(FILE);
}

if (!scalar(keys(%lambda_hash))) {
  print "No usable TI output files found - check presence of files in current directory make sure the correct extension is specified (--extension)\n";
  exit();
}

if ($warnings{lambdaconsistency}) {
  if (!$unsafe) {
    die ("\nInconsistency in lambda values found\nIndividual component lambdas (vdw and elec) should each match one-to-one with a single 'global' lambda value - this was found not to be the case. \nThis could be caused by using different values of tiElecLambdaStart or tiVdwLambdaEnd within one data set.\nRun with --unsafe to ignore\n.");
  }
}

#count timesteps for each value, print report if --fulloutput
my ($ts_min, $ts_max, $ts_avg, $total, $count);
if ($fulloutput) {if (!$dgonly) {print "\n\nNumber of MD steps per lambda value\n\n  lambda         MD steps\n";}}
foreach my $lval (sort {$a <=> $b} (keys(%timesteps))) {
  if ( !defined($ts_min) && !defined($ts_max)) {$ts_min = $ts_max = $timesteps{$lval};}
  if ($timesteps{$lval} < $ts_min) {$ts_min = $timesteps{$lval};}
  if ($timesteps{$lval} > $ts_max) {$ts_max = $timesteps{$lval};}
  $total += $timesteps{$lval}; $count ++;
  if ($fulloutput && !$dgonly) {
    printf("%9.7f\t% 8i\n",$lval, $timesteps{$lval});
  }
}
if (!$total) {print "\nApparently, no data has been collected from any of the output files. This could \nmean that the files are valid in structure but contain no data - or this is a \nbug in the script.\n"; exit();}

$ts_avg = $total / $count;
#discard any obviously undersampled data points (<25% of average nr), unless --unsafe
foreach my $lval (sort {$a <=> $b} (keys(%timesteps))) {
  if ($timesteps{$lval} < ($ts_avg * 0.25)) {
    if (!$unsafe) {
      delete($elec1{${$lambda_hash{$lval}}{elec1}});
      delete($elec2{${$lambda_hash{$lval}}{elec2}});
      delete($vdw1{${$lambda_hash{$lval}}{vdw1}});
      delete($vdw2{${$lambda_hash{$lval}}{vdw2}});
      delete($pme{$lval});
    }
    my @a; push @a, $lval, $timesteps{$lval};
    push @{$warnings{undersampled}}, \@a;
  }
}

# consistency check: are tiOutFreq frequencies consistent through all data?
my $prev;
foreach my $fval (@outputfreqs) {
  if (defined($prev) && ($fval != $prev)) {$warnings{outputfreq} = 1;}
  $prev = $fval;
}
if ($warnings{outputfreq}) {print "\nWarning: different tiOutFreq (output frequency) values were found. \nThis means some output files contain data sampled at a higher frequency than others.\nThis data will skew the average as the extra data points may contribute disproportionately to the average.\nBlock sampling (if requested) is disabled.\n";}

my    (%elec1avg, %elec2avg, %vdw1avg, %vdw2avg, %pmeavg); #key = lambda value; value = average dE/dl

if ($blocksample && !$warnings{$outputfreq}) {
  print "\nCalculating integrals for every $blocksample samples plus cumulative total\n";
  print "\n    Timestep     Block dG    Cumulative\n";
  # calculate running TI totals
  # value of blocksample = number of data points (each outputfreq steps) to sample at once
  my $ranOutOfData = 0; #keep going only until data at any lambda value runs out
  my $block = 0;
  my $cumulative=0;
  while (!$ranOutOfData) {
    foreach my $lval (keys(%elec1)) {$elec1avg{$lval} = Array_avg(\@{$elec1{$lval}},$blocksample,$block,\$ranOutOfData); }
    foreach my $lval (keys(%elec2)) {$elec2avg{$lval} = Array_avg(\@{$elec2{$lval}},$blocksample,$block,\$ranOutOfData); }
    foreach my $lval (keys(%vdw1)) {$vdw1avg{$lval} = Array_avg(\@{$vdw1{$lval}},$blocksample,$block,\$ranOutOfData); }
    foreach my $lval (keys(%vdw2)) {$vdw2avg{$lval} = Array_avg(\@{$vdw2{$lval}},$blocksample,$block,\$ranOutOfData); }
    foreach my $lval (keys(%pme)) {$pmeavg{$lval} = Array_avg(\@{$pme{$lval}},$blocksample,$block,\$ranOutOfData); }
    if (! $ranOutOfData) {
      my ($elec1_val, $elec1_report) = integrate_hash(\%elec1avg);
      my ($elec2_val, $elec2_report) = integrate_hash(\%elec2avg);
      my ($vdw1_val, $vdw1_report) = integrate_hash(\%vdw1avg);
      my ($vdw2_val, $vdw2_report) = integrate_hash(\%vdw2avg);
      #my ($pme_val, $pme_report) = integrate_hash(\%pmeavg);
      $cumulative += $elec1_val - $elec2_val + $vdw1_val - $vdw2_val; #+ $pme_val;
      printf "% 12i   % 10.4f   %10.4f\n",$blocksample*$block*$outputfreq,$elec1_val - $elec2_val + $vdw1_val - $vdw2_val,$cumulative/($block+1); 
    }
    $block ++;
  }
}


# average all the arrays for final output
foreach my $lval (sort {$a <=> $b} (keys(%elec1))) {$elec1avg{$lval} = Array_avg(\@{$elec1{$lval}}); } #print "elec1 lam $lval values ",scalar(@{$elec1{$lval}}),"\n";}
foreach my $lval (sort {$a <=> $b} (keys(%elec2))) {$elec2avg{$lval} = Array_avg(\@{$elec2{$lval}}); } #print "elec2 lam $lval values ",scalar(@{$elec2{$lval}}),"\n";}
foreach my $lval (sort {$a <=> $b} (keys(%vdw1))) {$vdw1avg{$lval} = Array_avg(\@{$vdw1{$lval}}); } #print "vdw1 lam $lval values ",scalar(@{$vdw1{$lval}}),"\n";}
foreach my $lval (sort {$a <=> $b} (keys(%vdw2))) {$vdw2avg{$lval} = Array_avg(\@{$vdw2{$lval}}); } #print "vdw2 lam $lval values ",scalar(@{$vdw2{$lval}}),"\n";}
#foreach my $lval (sort {$a <=> $b} (keys(%pme))) {$pmeavg{$lval} = Array_avg(\@{$pme{$lval}}); } #print "pme   lam $lval values ",scalar(@{$pme{$lval}}),"\n";}

#calculate integrals, store report data
if ($spline) {print "Spline interpolated electrostatics, partition 1\n";}
my ($elec1_val, $elec1_report) = integrate_hash(\%elec1avg);
if ($spline) {print "Spline interpolated electrostatics, partition 2\n";}
my ($elec2_val, $elec2_report) = integrate_hash(\%elec2avg);
if ($spline) {print "Spline interpolated vdW, partition 1\n";}
my ($vdw1_val, $vdw1_report) = integrate_hash(\%vdw1avg);
if ($spline) {print "Spline interpolated vdW, partition 2\n";}
my ($vdw2_val, $vdw2_report) = integrate_hash(\%vdw2avg);
#if ($spline) {print "Spline interpolated PME\n";}
#my ($pme_val, $pme_report) = integrate_hash(\%pmeavg);
my $timesteps_report;


if ($fulloutput && !$dgonly) { # print more detailed information about each energy component
  #print "\n\nPME Electrostatics\n$pme_report\n";
  print "Partition 1 electrostatics\n$elec1_report\n";
  print "Partition 1 vdW\n$vdw1_report\n";
  print "Partition 2 electrostatics\n$elec2_report\n";
  print "Partition 2 vdW\n$vdw2_report\n";
}

if (!$dgonly) { # main output section

  if ($warnings{undersampled}) {  #warning if undersampled data points were deleted
    my @us = @{$warnings{undersampled}};
    print "\n\n  Some data points had small number of samples (less than 25% of average):\n    lambda   timesteps\n";
    foreach my $aref (@us) {
      my @a = @{$aref};
      printf("  % 8.6f % 6i\n",$a[0],$a[1]);
    }
    if (!$unsafe) {
      print "  These data points were discarded - run with --unsafe to keep\n";
    }
    else {
      print "  The data points were kept due to the use of the --unsafe flag.\n";
    }
  }

  print  "\n|-----------------------------------------------|\n";
  print  "|         |    elec   |    vdW    |   Subtotal  |\n";
  print  "|-----------------------------------------------|\n";
  printf ("| Part. 1 |% 10.4f |% 10.4f | %10.4f  | \n",$elec1_val,$vdw1_val, $elec1_val + $vdw1_val);
  printf ("| Part. 2 |% 10.4f |% 10.4f | %10.4f  | \n", $elec2_val,$vdw2_val, $elec2_val + $vdw2_val);
  print  "|-----------------------------------------------|\n";
  printf ("| Subtotal|%10.4f |%10.4f |%10.4f   | \n", $elec1_val - $elec2_val,$vdw1_val - $vdw2_val, $elec1_val - $elec2_val + $vdw1_val - $vdw2_val);
  print  "|-----------------------------------------------|\n";
#  printf ("|   PME   |           |           |% 10.4f   | \n", $pme_val);
#  print  "|-----------------------------------------------|\n";
  printf "Total deltaG for transition lambda 0 -> 1: % 10.4f\n\n",$elec1_val - $elec2_val + $vdw1_val - $vdw2_val; #+ $pme_val; 

  print "Number of MD steps used for average:  min $ts_min max $ts_max average ",int($ts_avg);
  if (scalar(@outfiles) == 1) {print ", from 1 file";}
  else {print ", from ",scalar(@outfiles)," files";}
  if (!$tailsample) {print "\nFull dataset used, ";}
  if ($tailsample) {print "\n$tailsample\% of dataset used, ";}
  if (! $warnings{equilibrationdeleted}) {print "no initial equilibration data was discarded"; if (!$keepequilibration) {print " (none found)";}print "\n";}
  else {print "initial equilibration data was discarded\n";} 

  #some final consistency checks / warnings
  if ($warnings{sizelimit}) {print "\nWarning: one or more integrands were worryingly large (>$integralsizelimit); run with --fulloutput and check consistency.\n";}
  if ($warnings{lambdaconsistency}) { print "\nWarning: inconsistent lambda values - data is not self-consistent - disable --unsafe for more detail\n";} 

  # does the data cover the range lambda 0-1?
  my $max = 0; my $min = 1;
  foreach my $lval (keys(%lambda_hash)) {
    if ($lval > $max) {$max = $lval;}
    if ($lval < $min) {$min = $lval;}
  }
  if ($min == $max) {print "\nData only found for a single lambda value ($min) - no integral to calculate. Add data for more lambda values.\n";}
  else {
    if ($min != 0) {print "\nWarning: the first lambda value is not 0 - transition is incomplete. Starts from $min.\n";}
    if ($max != 1) {print "\nWarning: the last lambda value is not 1 - transition is incomplete. Ends at $max.\n\n";}
  }
} else {printf ("% 10.4f\n",$elec1_val - $elec2_val + $vdw1_val - $vdw2_val);}  #'dgonly' output, just the dG value
    

sub Array_avg {
  # return average of a given portion of an array (or all of it)
  # numblocks is the 
  my $count = 0; my $tot = 0;
  my @a = @{$_[0]};
  if (scalar(@_) == 1) { #just average the array
    foreach my $val (@a) {
      $tot += $val;
      $count ++;
    }
  return $tot / $count;
  }
  elsif (scalar(@_) == 4) { #means we want average only for a given block of the aray
    my $numblocks = $_[1];
    my $blocknumber = $_[2];
    if ($numblocks * $blocknumber + $numblocks > scalar(@a)) { #if we requested data beyond the end of the array...
      ${$_[3]} = 1; #arg [3] is reference to a flag indicating data has run out for at least one point
      return "outofdata";  
    }
    for (my $i = $numblocks * $blocknumber; $i < $numblocks * $blocknumber + $numblocks; $i++) { #otherwise average it
      $tot += $a[$i];
      $count ++;
    }
  return $tot / $count;
  }
}

sub integrate_hash {  # integration from hash with key=lambda, value = dE/dl ensemble avg
  my %hash = %{$_[0]};
  my @lambdas = sort {$a <=> $b} keys(%hash);
  my $report;
  if (scalar(@lambdas) <= 1) {
    $report = "\n  Lambda\t  dE/dl\t          Area\t         Total\n";
    if (scalar(@lambdas)) {
      $report .= sprintf("%9.7f\t%8.6f\t%8s\t%8.6f\n",$lambdas[0], $hash{$lambdas[0]}, "N/A", $total); 
    }
    else {
      $report .= sprintf("%9.7f\t%8.6f\t%8s\t%8.6f\n",0, 0, "N/A", $total);     
    }
    return 0, $report;
  }
  my @dedl;
  foreach my $lamval (@lambdas) {push @dedl, $hash{$lamval};}
  my $total = 0;
  $report = "\n  Lambda\t  dE/dl\t          Area\t         Total\n";

  if ($oldtrapezoidal) { # old trapezoidal integration, as far as i can see this is worse than default (spline) but still available
    for (my $i = 0; $i < scalar(@lambdas) - 1; $i ++) {
      my $area;
      my $lambda = $lambdas[$i];
      my $nextlambda = $lambdas[$i+1];
      my $currentavg = $hash{$lambda};
      my $nextavg = $hash{$nextlambda};
      my $current = 0;
      $current += ($currentavg * ($nextlambda - $lambda));
      $current -= ($currentavg - $nextavg) * ($nextlambda - $lambda) / 2;
      $area += $current;
      $report .= sprintf("%9.7f\t%8.6f\t%8.6f\t%8.6f\t\n",$lambda, $currentavg, $current, $total + $current);
      if ($area >= $integralsizelimit || $area <= $integralsizelimit * (-1)) {
        $warnings{sizelimit} = 1;
        $total += $area;
      }
      else {
        $total += $area;
      }
    }
    $report .= sprintf("%9.7f\t%8.6f\t%8s\t%8.6f\t\n",$lambdas[scalar(@lambdas)-1], $hash{$lambdas[scalar(@lambdas)-1]}, "N/A", $total); 
    return $total, $report;
  }

  else {  #cubic spline interpolation
  
    # first derivatives
    # these are estimated crudely... slope of the line from previous dU/dl to the next one
    # can do better...
    my @firstderiv;
    my $num = scalar(@lambdas);
    $firstderiv[0] = ($dedl[1] - $dedl[0])/($lambdas[1]-$lambdas[0]);
    $firstderiv[$num-1] = ($dedl[$num-1] - $dedl[$num-2])/($lambdas[$num-1]-$lambdas[$num-2]);
    for (my $i = 1; $i < $num - 1; $i++) {
      $firstderiv[$i] = ($dedl[$i+1]-$dedl[$i-1]) / ($lambdas[$i+1] - $lambdas[$i-1]);
    }
    # cubic spline code stolen straight from NAMD's nonbondedUtil... 
    # long term project is to implement a smooting spline, weighted according to number of points per sample
    my (@aa,@ab,@ac,@ad);
    for (my $i = 0; $i < $num-1; $i++) {
      my $x = $lambdas[$i+1] - $lambdas[$i];
      my $a = $dedl[$i];
      my $b = $firstderiv[$i];
      my $anext = $dedl[$i+1];
      my $bnext = $firstderiv[$i+1];
      my $c = ( 3.0 * ($anext - $a) - $x * (2.0 * $b + $bnext) ) / ( $x * $x );
      my $d = ( -2.0 * ($anext - $a) + $x * ($b + $bnext) ) / ( $x * $x * $x );
      push @aa, $a;
      push @ab, $b;
      push @ac, $c;
      push @ad, $d;
    }
    push @aa, $dedl[$num-1];
    push @ab, $firstderiv[$num-1];
    push @ac, 0;
    push @ad, 0;
  
    # integrate (analytical integration of each of the cubic spline portions)
    for (my $i = 0; $i < scalar(@lambdas)-1; $i++) {
      my $lambda=$lambdas[$i];
      my $nextlambda=$lambdas[$i+1];
      my $a = $aa[$i];  my $b = $ab[$i];  my $c = $ac[$i];  my $d = $ad[$i];
      my $integral = (($nextlambda - $lambda)*(12*$a + ($nextlambda - $lambda)*(6*$b + ($nextlambda - $lambda)*(4*$c + 3*$d*$nextlambda - 3*$d*$lambda))))/12;
      $report .= sprintf("%9.7f\t%8.6f\t%8.6f\t%8.6f\n",$lambda, $dedl[$i], $integral, $total+$integral);
      if ($integral >= $integralsizelimit || $integral <= $integralsizelimit * (-1)) {
        $warnings{sizelimit} = 1;
        $total += $integral;
      }
      else {
        $total += $integral;
      }
    }
    $report .= sprintf("%9.7f\t%8.6f\t%8s\t%8.6f\n",$lambdas[scalar(@lambdas)-1], $dedl[scalar(@lambdas)-1], "N/A", $total); 
    if (!$spline) {return $total, $report;}
    
    # print interpolated values if requested
    my @int_x;
    for (my $i = 0; $i < 1; $i += (1/$spline)) {
      push @int_x, $i;
    }
    push @int_x,1;
    foreach my $i (@int_x) {
      my $xbelow;
      my $indexbelow;
      for (my $j = 0; $j < scalar(@lambdas); $j++) {
        if (!defined($xbelow)) {$xbelow = $lambdas[$j]; $indexbelow = 0;}
        elsif ($i >= $lambdas[$j]) {
          $xbelow = $lambdas[$j]; $indexbelow = $j;
        }
      }
      my $diffa = $i - $xbelow;
  
      my $a = $aa[$indexbelow];
      my $b = $ab[$indexbelow];
      my $c = $ac[$indexbelow];
      my $d = $ad[$indexbelow];
      my $val = $d*$diffa**3 + $c*$diffa**2 + $b*$diffa + $a;
      printf ("INTERPOLATED: %8.7f %8.4f\n",$i,$val);    
    }
    return $total, $report;
  }
}

  


sub helpmsg {
  print "\nScript for calculation of free energy change for an alchemical
transformation using the Thermoydynamic Integration functionality of NAMD

estimate the integral of dU/dlambda from lambda 0 -> 1
  
Input is pretty flexible - dU/dl data can all be in a single file describing the
whole transition, or in multiple files with a single lamda value each, or in 
multiple files, each containing data for single or multiple lambda values.

This flexibility means it's easy to add / remove data points from the calculation
- for example, you might want to do a quick initial run for a small set of lambda
values and realise you need more sampling in a certain range - you can simply add
more output files and the data will be included.

However, care should be taken that data is self-consistent, i.e. produced with 
the same global parameters for TI calculations (tiOutFreq, tiElecLambdaStart, 
tiVdwLambdaEnd etc.). Some consistency checks are performed but don't rely on 
these to catch everything. If the script doesn't like something in the data 
that you know to be acceptable, you can generally overrule the consistency 
check with the option -u (--unsafe).

Undersampled data points - if a file is found that contains far less data than 
the rest of the dataset (arbitrarily set at 25% of the average for all files), 
the data is discarded, as poorly sampled data is judged to reduce the quality 
of the final dG estimate. Overrule with -u (--unsafe).

See the NAMD user's guide for more detail about the meaning of the data in the 
TI output files and the meaning of this script's output.

File input: default is to include all files ending 'ti.out'. Alternatively
files can be included from the command line, e.g. 
  \"NAMD_ti.pl alch*ti.out\"
Specifying a directory will include all files ending 'ti.out' in that directory.

  Command line options:

  -k --keepequilibration
  NAMD offers the option to discard the first x number of MD steps after 
  starting or changing lambda value, as specified by 'tiEquilSteps' in the 
  NAMD configuration file. If you want to include these values in the average 
  after all, specify -k or --keepequilibration.

  -t --tailsample
  Offers the option to use only the last x percent of the data for each lambda
  value for calculating the average, giving more flexibility in which data 
  should be regarded as equilibrated and counted towards the average. For 
  example, specifying -t=60 would throw out the first 40% of the data and use 
  the remaining 60% to calculate the average.
  --tailsample and --keepequilibration may be used together. Default (-k not 
  specified) is first to throw out the equilibration data specified by 
  'tiequilsteps', and then to use the last x percent of the *remaining* data. 
  Using -k -t=60 would use the last 60% of *all* the data.
  Care should be taken if data for a single lambda value is spread over 
  multiple files: if data isn't read in the same order as it was produced, the 
  'tail' data used might not be the latest data. File input is conducted 
  according to file modification date with the newest files read last, so block
  sampling can be used successfully across multiple files as long as newer data
  has a more recent timestamp than older data. 

  -b --blocksample
  Calculates rolling TI dG totals and cumulative averages over the course of a 
  run. For example, running with -b=10, a TI total plus cumulative average will
  be calcuated for every 10 data points (corresponding to 10*tiOutFreq MD 
  steps). This output can be graphed to assess variance of the data and 
  convergence of the average.
  This type of block sampling over a time course is designed for data sets 
  where each lambda value is simulated in a separate run from a single starting
  point - the starting structures for each lambda value should in principle be 
  identical. You might or might not get something meaningful out with other 
  sampling schemes but this is not recommended and results should be 
  interpreted with care!
  Block sampling will stop once the end of the shortest dataset is reached
  - e.g. if a single lambda value has only half as much data as the rest of the
  dataset, block analysis will stop halfway.
  Care should be taken if data for a single lambda value is spread over 
  multiple files (see 'tailsample' above) - oldest files are read first. 
  
  -u --unsafe
  Ignore failed data consistency checks, and don't discard any data points that 
  were judged to be undersampled. If this option is necessary, it's probably 
  wise to look at the underlying causes.

  -f --fulloutput
  Display a detailed summary of the data, with dU/dl ensemble averages for each 
  value of lambda for each of the components. Unless you're very confident about
  what you're doing, it's always worthwhile to graph dU/dl against lambda and 
  see how the integral looks - this will help judge if finer or coarser sampling
  would be better at any stage of the transformation.
  
  -d --dgonly
  Output nothing but the final dG value - if you're feeling lucky

  -o --oldtrapezoidal
  Integration of dU/dl values is by default done by means of a cubic spline 
  integration. Running with -o uses a simpler trapezoidal integration instead. 
  The cubic spline method seems to cope better with less than perfect data but 
  this hasn't been tested extensively. There are lots of ways to estimate 
  integrals for arbitrary data, some without doubt better than what's 
  implemented here, so users are encouraged to experiment.

  -s --spline
  Outputs a list of interpolated lambda vs. dE/dl values for each energy 
  component, from a cubic spline interpolation. The default dG output (when 
  trapezoidal integration is not requested) is the integral of this curve 
  (calculated analytically). Specify how many points to write out with -s=x.

  -h --help
  Display this message.\n\n";
}
