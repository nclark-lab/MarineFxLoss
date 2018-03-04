sub in_fasta {
	#input array of lines
	#output hash reference of sequences
	my %seq;
	my $name='';
	my $s='';	
	foreach my $line (@_) {
		chomp $line;
		if ($line =~ /^>(\S+)/) {
			$seq{$name} = $s if $name;
			$name = $1;
			$s = '';
		} else {
			$line =~ s/\s+//g;
			$s .= $line;
		}
	}
	$seq{$name} = $s;
	return \%seq;
}

sub out_fasta {
	#input: reference to hash containing multiple sequences (not required to be aligned or even same length)
	#output: a fasta formatted string, which is ready for printing to screen or writing to a file.
	my %seq = %{$_[0]};
	my $out = '';
	foreach my $k (keys %seq) {
		$out .= ">$k\n";
		$out .= "$seq{$k}\n";
	}
	return $out;
}

sub pruneTree {
    my  $tree=$_[0];
    my $l=$_[1];
    my $numclass='1234567890eE\.\-+';
    my ($d1,$d2, $posE, $posS, $pcount);
    my ($nstr, $nstrs, $dsum);

# $l is to the right
    if ($tree=~/:([$numclass]*),$l:[$numclass]*\):([$numclass]*)/){
	$d1=$1;
	$d2=$2;
	$posE=$-[0];
	$posS=$posE;
	$pcount=0;
	while($pcount>=0){
	    $posS--;  
	    if(substr($tree, $posS,1)eq")"){
		$pcount++;
	    }
	    if(substr($tree, $posS,1)eq"("){
	$pcount--;
	    }
	}
#nstr is the string for the neighbor tree
	$nstr=substr($tree, $posS+1, $posE-$posS-1);
	$nstrs=quotemeta $nstr;
	$dsum=$d1+$d2;
	$tree=~s/\($nstrs:[$numclass]*,$l:[$numclass]*\):[$numclass]*/$nstr:$dsum/;

    }
#$l is on left
    elsif ($tree =~/\($l:[$numclass]*,/ ){
	$posS=$+[0];
	$posE=$posS;
	$pcount=0;
	while($pcount>=0){
	    $posE++;  
	    if(substr($tree, $posE,1)eq"("){
		$pcount++;
	    }
	    if(substr($tree, $posE,1)eq")"){
	$pcount--;
	    }
	}
	
	
	$nstr=substr($tree, $posS, $posE-$posS+1);

#remove distance if neighbor is a single leaf
	if(substr($nstr,0,1) ne "("){
	    $nstr=~s/:[$numclass]*\)//;	
	}
	$nstrs=quotemeta $nstr;
	$tree=~/\($l:[$numclass]*,$nstrs:([$numclass]*)\):([$numclass]*)/;
	$d1=$1;
	$d2=$2;
	$dsum=$d1+$d2;
	$tree=~s/\($l:[$numclass]*,$nstrs:[$numclass]*\):[$numclass]*/$nstr:$dsum/;

    
    }
    
    else{
	print STDERR "Warning: no leaf $l, nothing to do\n";
    }
     return($tree)
}
 


sub out_paml {
	#PAML version of phylip format, I allow names up to 30 characters
	my $hash = $_[0];
	my @keys = sort keys %{$hash};
	my $n_seqs = scalar @keys;
	#Count non-white space length
	my $counter = $hash->{$keys[0]};
	$counter =~ s/\s+//g;
	my $length = length $counter;
	my $out = " $n_seqs  $length\n";
	foreach (@keys) {
		my $name = $_;
		$name = substr $_, 0, 30 if length $_ > 30;
		my $space = " " x (32 - length $name);
		$out .= "$name$space$hash->{$_}\n";
	}
	return $out;
}

sub range {
	#Given array, rank and return max and min values.
	return 'na' unless @_;
	my @e = sort {$a <=> $b} @_;	#sort from highest to lowest
	return ($e[0] , $e[-1]);
}


sub max_min {
	#Given array, rank and return max and min values.
	return 'na' unless @_;
	my @e = sort {$b <=> $a} @_;	#sort from highest to lowest
	return ($e[0] , $e[-1]);
}

######################
sub composite_LD {
	#based on Weir, BS (1996). Genetic Data Analysis II (Sunderland, MA: Sinauer Associates, Inc.).
	#I took it from Rohlfs, Swanson, Weir. 2010. AJHG 86, 674-685.
######################
	#Arguments: Array of 9 positive integers. The genotype counts in contingency table:
	#	AABB, AABb, AAbb		0 , 1 , 2
	#	AaBB, AaBb, Aabb		3 , 4 , 5
	#	aaBB, aaBb, aabb		6 , 7 , 8
	
	my @g = @_; #Array of 9 genotypes
	
	#Count total counted genotypes. (Chromosomes = 2N)
	my $N = 0;
	foreach (@g) { $N += $_; }
	return 'na' unless $N;

	#Determine allele freq locus A and B
	my $p1 = (2 * ($g[0] + $g[1] + $g[2]) + $g[3] + $g[4] + $g[5]) / (2*$N) ;
	my $q1 = (2 * ($g[0] + $g[3] + $g[6]) + $g[1] + $g[4] + $g[7]) / (2*$N) ;

	my $deltaAB = (2*$g[0] + $g[1] + $g[3] + $g[4]/2) / $N - 2*$p1*$q1;
	my $p11 = ($g[0] + $g[1] + $g[2]) / $N;
	my $q11 = ($g[0] + $g[3] + $g[6]) / $N;
	return 'na' if (($p1*(1-$p1) + $p11 - $p1**2) == 0 || ($q1*(1-$q1) + $q11 - $q1**2) == 0);
	my $chi2 = $N * $deltaAB**2 / ($p1*(1-$p1) + $p11 - $p1**2) / ($q1*(1-$q1) + $q11 - $q1**2);
	
	return ($chi2); #One degree of freedom
}

######################
sub genotype_association {
	#I took this from Rohlfs, Swanson, Weir. 2010. AJHG 86, 674-685.
	#Chi-square with 4 degrees of freedom.
	#Returns 'na' in the case of division by zero from expected values of zero.
######################
	#Arguments: Array of 9 positive integers. The genotype counts in contingency table:
	#	AABB, AABb, AAbb		0 , 1 , 2
	#	AaBB, AaBb, Aabb		3 , 4 , 5
	#	aaBB, aaBb, aabb		6 , 7 , 8
	
	#Count total counted genotypes.
	my $N = 0;
	foreach (@_) { $N += $_; }
	return 'na' unless $N;

	my @g = map({$_ / $N} @_); #Array of 9 genotype frequencies

	#Determine marginal genotype freqs: locus A
	my $p11 = ($g[0] + $g[1] + $g[2]);
	my $p12 = ($g[3] + $g[4] + $g[5]);
	my $p22 = 1 - $p11 - $p12;

	#Determine marginal genotype freqs: locus B
	my $q11 = ($g[0] + $g[3] + $g[6]);
	my $q12 = ($g[1] + $g[4] + $g[7]);
	my $q22 = 1 - $q11 - $q12;
	
	#Compute expected genotype frequencies.
	my @exp = ( $p11 * $q11 , $p11 * $q12 , $p11 * $q22 ,
				$p12 * $q11 , $p12 * $q12 , $p12 * $q22 ,
				$p22 * $q11 , $p22 * $q12 , $p22 * $q22 );
	
	foreach my $i (0..8) {
		return 'na' if $exp[$i] == 0;
	}
	my $chi2 = 0;
	foreach my $i (0..8) {
		$chi2 += ($g[$i] - $exp[$i])**2 / $exp[$i];
	}

	return ($chi2 * $N); #Four degrees of freedom	
}


sub rank_percentile {
	#Determine the rank of a value in a distribution
	#Args: value, @distribution
	#Returns percentile (more positive = closer to 0) "Top 1%"

	return ("NA") if scalar @_ < 3;
	my $value = shift @_;
	my @distri = sort {$b <=> $a} @_;
	my $n_distri = scalar @distri;

	my $tally = 0;
	foreach (@distri) {
		last if $value >= $_;
		$tally++;
	}
	return ($tally / $n_distri);
}

sub rank_percentile_presorted {
	#Determine the rank of a value in a distribution
	#Args: value, @distribution already sorted to start with highest value (i.e. {$b <=> $a}
	#Returns percentile (more positive = closer to 0) "Top 1%"

	return ("NA") if scalar @_ < 3;
	my $value = shift @_;
	my @distri = @_;
	my $n_distri = scalar @distri;

	my $tally = 0;
	foreach (@distri) {
		last if $value >= $_;
		$tally++;
	}
	return ($tally / $n_distri);
}

sub input_to_hash {
	# Reads in a tab-delimited file and stores 2 columns in a hash.
	# Specify the file path, and the two columns in arguments.
	# Returns a hash, not a reference.
	# Arguments: (1) file path, (2) element# for key, (3) element# for value
	# Usage example: my %hash = input_to_hash("$ARGV[0]", 0, 1);
	my %hash;
	foreach my $line (read_input($_[0])) {
		chomp $line;
		next if $line =~ /^\s*$/;
		next if $line =~ /^\s*#/;
		my @e = split /\t/, $line;
		$hash{ $e[$_[1]] } = $e[$_[2]] ;
	}
	return %hash;
}
sub input_to_hash_ref {
	# Reads in a tab-delimited file and stores 2 columns in a hash.
	# Specify the file path, and the two columns in arguments.
	# Returns a hash, not a reference.
	# Arguments: (1) file path, (2) element# for key, (3) element# for value, (4) delimiter, optional, default "\t"
	# Usage example: my $hashref = input_to_hash("$ARGV[0]", 0, 1);
	my %hash;
	foreach my $line (read_input($_[0])) {
		my $delim = "\t";
		if ($_[3]) { $delim = $_[3]; }
		chomp $line;
		next if $line =~ /^\s*$/;
		next if $line =~ /^\s*#/;
		my @e = split /$delim/, $line;
		$hash{ $e[$_[1]] } = $e[$_[2]] ;
	}
	return \%hash;
}

sub table_2Dhash {
	# Reads in a tab-delimited file and stores table in a 2D hash.
	# Specify the file path
	# Returns a reference to a 2D hash.
	# Arguments: (1) file path
	my %hash;
	my @col;
	foreach my $line (read_input_clean($_[0])) {
		my @e = split /\t/, $line;
		if (@col) {
			my $row = shift @e;
			for (my $i=0; $i<scalar @e; $i++) {
				$hash{$row}->{$col[$i]} = $e[$i];
			}
		} else {
			shift @e;
			@col = @e;
		}
	}
	return \%hash;
}

sub phylip_sequential2fasta {
	my @fasta;
	my ($number, $length, $seq) = in_phylip_sequential(read_input($_[0]));
	foreach my $k (sort keys %{$seq}) {
		push @fasta, ">$k", $seq->{$k};
	}
	return (\@fasta);
}

sub limit_stdev {
	#@_ = 0each value, 1pop mean, 2pop stdev, 3imposed limit.
	#normalize mean
	$_[0] -= $_[1];

	#limit
	if (abs($_[0] / $_[2]) > $_[3]) {
		return (sign_indicator($_[0]) * $_[2] * $_[3] );
	} else {
		return $_[0];
	}
}

sub sign_indicator {
	if ($_[0] < 0) {
		return (-1);
	} else {
		return (1);
	}
}

sub regression_standalone {
	# This subroutine returns the y-intercept ($q), the slope ($m) and correlation coefficient ($r) for a given set of X and Y values.
	#
	# Argument 1 is a reference to an array of values for variable X
	# Argument 2 is a reference to an array of values for variable Y
	
	my @y = @{$_[0]};
	my @x = @{$_[1]};
	
	#sums and sums of squares, ...etc
	my $N = scalar @x;
	my $sum_x = 0;
	$sum_x += $_ foreach @x;
	
	my $sum_y = 0;
	$sum_y += $_ foreach @y;

	if ($sum_x == 0) {		#Adds pseudo count to any sum of x that is zero. This is b/c some future calculations divide by the sum of x.
		$sum_x = 0.0000001;
	}
	
	my $sum_xy;
	for (0..$N-1) {
		$sum_xy += $x[$_] * $y[$_];
	}
	
	my $sum_x2 = 0;
	$sum_x2 += $_**2 foreach @x;
	
	my $sum_y2 = 0;
	$sum_y2 += $_**2 foreach @y;
	
	#calculate regression line with least squares fit
	my $q = ( $sum_y/$sum_x - $sum_xy/$sum_x2) / ( $N/$sum_x - $sum_x/$sum_x2 ); #intercept
	my $m = ( $sum_y - $N * $q) / $sum_x;#slope

	my @yc = map($q + $m * $_ , @x); #predicted Y values from correlation line
	my $sum_d2 = 0;	#sum of squared difference between y and yc (rms?)
	for (0..$N-1) {
		$sum_d2 += ($y[$_] - $yc[$_])**2;
	}
	
	my $Sy = ($sum_d2 / $N)**0.5; #Standard error of the estimate
	my $sigma_y = ( $sum_y2 / $N - ($sum_y/$N)**2 )**0.5; #standard deviation of y
	
	my $pre_r = (1 - $Sy**2 / $sigma_y**2 );
	if ($pre_r < 0) {
		$pre_r = 0;
	}
	my $r = $pre_r**0.5;
	$r *= -1 if $m < 0; #give proper sign to correlation coefficient
	
	return ($q, $m, $r);
}

sub bin2dec {
	# Convert a binary number to a decimal number
	my @e = split //, $_[0];
	my $sum = 0;
	my $factor = 1;
	while (@e) {
		$sum += $factor * pop @e;		
		$factor *= 2;
	}	
	return $sum;
}

sub outputText {
	my ($path, @output) = @_;
	open OUT, "$path" or die "Can't open $path for output.\n";
	print OUT @output;
	close OUT;
}


sub draw_without_replacement {
	#Arguments: number of draws, reference to array from which to draw.
	my @draws;
	my @background = @{$_[1]};
	for(1..$_[0]) {
		push @draws, (splice @background, (rand @background), 1);
	}
	return @draws;
}

sub median {
	return ('na') unless scalar @_;
	#Give an array of values. Returns the median value.
	my @s = sort { $a <=> $b } @_;
	return (($s[scalar @s / 2] + $s[(scalar @s - 1) / 2]) / 2);
}

sub mean_variance {
	return ('na','na') unless scalar @_;
	my $sum = 0;
	$sum += $_ foreach @_;
	my $mean = $sum / scalar @_;
	
	#MSE - mean squared error
	my $var = 0;
	$var += ($_ - $mean)**2 foreach @_;
	$var /= (scalar @_ - 1);
	return ($mean, $var);
}

sub contingency_table_permutations {
	# Table:  ________
	#		  |a | b |
	#		  --------
	#		  |c | d |
	#		  --------
	my ($ta, $tb, $tc, $td, $permutations) = @_;
	$permutations = 100000 unless $_[4];  #Default of 100-thousand permutations if none specified.
	
	#Set counts for tests >= 'a' and 'c'
	my $count_a = 0;
	my $count_c = 0;
	
	my $ac = $ta + $tc;  #marginal
	foreach (1..$permutations) {
		#Set/re-set marginals that are subtracted from.
		my $ab = $ta + $tb;
		my $cd = $tc + $td;

		#Set simulated 'a' and 'c' cells.
		my ($ia, $ic) = (0,0);
		
		for (1..$ac) {
			my $freq_ab = $ab / ($ab + $cd);
			if (rand() < $freq_ab) {
				$ia++;
				$ab--;
			} else {
				$ic++;
				$cd--;
			}
		}
		$count_a++ if $ia >= $ta;
		$count_c++ if $ic >= $tc;
	}

	my $p_left = $count_c / $permutations;
	my $p_right= $count_a / $permutations;
	my @p =   sort {$a <=> $b}  ($p_left, $p_right);
	my $p_two = $p[0] + (1 - $p[1]);   #I'm not sure what this p-value means.
	return ($p_left , $p_right, $p_two);
}


sub arithmetic_mean {
	return ('na') unless scalar @_;
	my $sum = 0;
	foreach my $value (@_) {
		$sum += $value;
	}
	return ($sum/(scalar @_));
}

sub summation {
	return ('na') unless scalar @_;
	my $sum = 0;
	foreach my $value (@_) {
		$sum += $value;
	}
	return $sum;
}

sub shuffle {  # randomly reorders the array thru draws without replacement
	my @s;
	while(@_) {
		push @s, (splice @_, (rand @_), 1);
	}
	return @s;
}

sub bootstrap_array {  #Bootstraps elements in array to equal the number in the original array.
						# Similar to shuffle, except it draws *with* replacement.
	my @b;
	do { 
		push @b, $_[rand @_];
	} until (scalar @b == scalar @_);
	return @b;
}

######################
sub read_first_fasta {
######################
	my $seqFile = '';
	foreach (read_input($_[0])) { #Remove all comment lines and blank lines.
		$seqFile .= $_ unless $_ =~ m/^\s+$|^\s*#/;
	}
	my @seqFile = split /^>/m, $seqFile;
	my @first = split /\n/, $seqFile[1];
	my $title = ">" . shift @first;
	my $seq = join '', @first;
	$seq =~ s/\s//g;
	return ($title, $seq);
}

############
sub prompt {		#Returns STDIN without newline.
############
	my $prompt = shift @_;
	print STDERR $prompt;
	my $response = <STDIN>;
	chomp $response;
	return $response;
}

sub prompt_choice {	#Returns STDIN only if an appropriate choice is made.
	#Arguments are ["prompt"] ["choices"]...
	my $prompt = shift @_;
	my @choices = @_;
	my ($response, $match, $before);
	until ($match) {
		print STDERR "\nSorry, your response is not an appropriate choice.\n" if $before;
		print STDERR $prompt;
		$response = <STDIN>;
		chomp $response;
		foreach (@choices) {
			$match = 1 if $_ eq $response;
		}
		$before = 1;
	}
	return $response;
}



sub getHomeDir {
	my @env = `env`;
	my $env = join '', @env;
	my $homeDir = 0;
	$homeDir = $1 if $env =~ /^HOME=(\S+)$/m;
	return $homeDir;
}

sub output_string { #same as output_text
	my ($path, @output) = @_;
	open OUT, "$path" or die "Can't open $path for output.\n";
	print OUT @output;
	close OUT;
}
sub output_text { #same as output_string
	my ($path, @output) = @_;
	open OUT, "$path" or die "Can't open $path for output.\n";
	print OUT @output;
	close OUT;
}
sub output_array {
	my ($path, @output) = @_;
	open OUT, "$path" or die "Can't open $path for output.\n";
	print OUT $_,"\n" foreach @output;
	close OUT;
}


################
sub readINPUT {
################
	my ($inputName) = @_;
	open INPUT, $inputName or die "Could not open $inputName";
	my @input = <INPUT>;
	close INPUT;
	return @input;
}
################
sub read_input {
################
	my ($inputName) = @_;
	open INPUT, $inputName or die "Could not open $inputName";
	my @input = <INPUT>;
	close INPUT;
	return @input;
}

################
sub read_input_unix {
################
	my ($inputName) = @_;
	open INPUT, $inputName or die "Could not open $inputName";
	my $string = join '' , <INPUT>;
	close INPUT;
	$string =~ s/\r\n|\n|\r|\cM/\n/g;	#convert DOS, old Mac, and MacOSX newlines to Unix newlines
	my @input = split /\n/ , $string;
	return @input;
}


################
sub read_input_clean {
################
	my ($inputName) = @_;
	open INPUT, $inputName or die "Could not open $inputName";
	my @input;
	while (<INPUT>) {
		next if /^#/;
		next if /^\s*$/;
		chomp;
		push @input, $_;
	}
	close INPUT;
	return @input;
}

################
sub getARGUMENTS {
################
	my $option = '';
	my @A;
	$optWithArg = shift @_;
	
	foreach (@_) {
		if ($option) {
			$option = '';
			next;
		}
		if ($_=~/^-\S+/) {
			$option = $_;
			unless ($optWithArg =~ /$option-/) {
				$option = '';
			}
		}else {
			push @A, $_;
		}
	}
	return (@A);
}

################
sub getOPTIONS {
################
	my $option = '';
	my %options;
	$optWithArg = shift @_;

	foreach (@_) {
		if ($option) {
			$options{$option} = $_;
			$option = '';
			next;
		}
		if ($_=~/^-\S+/) {
			$option = $_;
			unless ($optWithArg =~ /$option-/) {
				$options{$option} = 1;
				$option = '';
			}
		}
	}
	return (%options);
}

###############
sub in_phylip_interleaved {
###############
	my %hash;
	my @phylip = @_;
	my (@names, @seq);
	my $numberLine = shift @phylip;
	my ($number, $length) = ($1,$2) if $numberLine =~ /^\s*(\d+)\s+(\d+)/;
	#title lines
	foreach (0..$number-1) {
		if ($phylip[$_] =~ /^(\S+)\s+(\S.+)$/) {
			push @names, $1;
			push @seq, $2;
		}
	}
	#rest of the sequence
	my $i = 0;
	foreach ($number..(scalar@phylip)-1) {
		next unless $phylip[$_] =~ /\S/;
		$i = 0 if $i == $number;
		$seq[$i] .= $phylip[$_];
		$i++;
	}
	#remove whitespace from seq
	foreach (@seq) {
		$_ =~ s/\s//g;
	}
	for (0..$number-1) {
		$hash{$names[$_]} = $seq[$_];
	}
	return ($number, $length, \%hash);
}

########################
sub in_phylip_sequential {
########################
	my %hash;
	my @alignment;
	foreach (@_) {
		push @alignment, $_ unless $_ =~ /^\s*$/;
	}
	
	my $numberLine = shift @alignment;
	my ($number, $length) = ($1,$2) if $numberLine =~ /\s*(\d+)\s+(\d+)/;
	#sequence lines
	my $file = join '', @alignment;
	my @char = split m//, $file;

	for (1..$number) {
		my $name = '';
		my $this_char = '';
		do {
			$this_char = shift @char;
			$name = $this_char unless $this_char =~ /\s/;
		} until ($name);
		
		do {
			$this_char = shift @char;
			$name .= $this_char unless $this_char =~ /\s/;
		} until ($this_char =~ /\s/);
		
		
		my $seq = '';
		do {
			$this_char = shift @char;
			$seq .= $this_char unless $this_char =~ /\s/;
		} until (length $seq == $length);

		$hash{$name} = $seq;
	}

	return ($number, $length, \%hash);
}

sub out_phylip_sequential {
	#Pass a hash reference of the sequences.
	#Returns seq phylip as a STRING.
	my $hash = $_[0];
	my @keys = sort keys %{$hash};
	my $n_seqs = scalar @keys;
	#Count non-white space length
	my $counter = $hash->{$keys[0]};
	$counter =~ s/\s+//g;
	my $length = length $counter;
	my $out = " $n_seqs  $length\n";
	foreach (@keys) {
		my $name = $_;
		$name = substr $_, 0, 8 if length $_ > 8;
		my $space = " " x (10 - length $name);
		$out .= "$name$space$hash->{$_}\n";
	}
	return $out;
}


############
sub revCOMP
############
{	my ($sequence) = @_;
	$sequence =~ tr /ACGTacgt/TGCAtgca/;
	$sequence = reverse $sequence;
	return($sequence);
}


#########################################
sub calcPERCENT
#Calculates the percentage of base(s) in sequence.
#########################################
{	my ($sequence, $base) = @_;
	my $count = 0;
	@sequence = split '', uc $sequence;
	@bases = split '', uc $base;

	foreach $b (@bases)
	{	foreach $s (@sequence)
		{	if ($s eq $b)
			{	$count++
			}
		}
	}
	
	$percent = ($count / (length $sequence))*100;
	return($percent);
}

###############
sub allPERCENT {
###############
	my ($sequence) = @_;
	@sequence = split '', uc $sequence;
	@bases = (A,C,G,T,U,M,R,W,S,Y,K,V,H,D,B,N);
	
	
	foreach $b (@bases) {
		$count{$b} = 0;
		foreach $s (@sequence) {
			if ($s eq $b) {
				$count{$b}++
			}
		}
	}
	
	return (%count);
}


###############
sub allPERCENT_AA {
###############
	my ($sequence) = @_;
	@sequence = split '', uc $sequence;
	@bases = (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y);	
	
	foreach $b (@bases) {
		$count{$b} = 0;
		foreach $s (@sequence) {
			if ($s eq $b) {
				$count{$b}++
			}
		}
	}
	
	return (%count);
}
	
##################
sub whichFRAME
		#determines the best reading frame and returns it with any stops within it.
##################
{	my ($thisSeq ) = @_;
	my (@one, @two, @three, @stops) = ();
	
	##print "\nSUB:TEST: length: " , length $thisSeq , "\n";
	#Count stop codons.
	for ($j = 0; $j < (length $thisSeq)-2 ; $j += 3)
	{
		$thisCodon = substr $thisSeq, $j, 3;
		##print $thisCodon;
		if ($thisCodon =~ /TAA|TAG|TGA/i)
		{	push @one, ($j+1);
		}
		$thisCodon = substr $thisSeq, ($j+1), 3;
		if ($thisCodon =~ /TAA|TAG|TGA/i)
		{	push @two, ($j+2);
		}
		$thisCodon = substr $thisSeq, ($j+2), 3;
		if ($thisCodon =~ /TAA|TAG|TGA/i)
		{	push @three, ($j+3);
		}
	}

	$frameDisplay = " one : @one\n two : @two\nthree: @three";####$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#

	#determine which readingframe by least stop codons.
	$n1 = scalar @one;	$n2 = scalar @two;	$n3 = scalar @three;
	if ($n1 < $n2 && $n1 < $n3) {	$frame = 1;	push @stops, @one;}
	if ($n2 < $n1 && $n2 < $n3) {	$frame = 2;	push @stops, @two;}
	if ($n3 < $n1 && $n3 < $n2) {	$frame = 3;	push @stops, @three;}
	if ($n1 == $n2 && $n1 < $n3) {	$frame = 12; push @stops, @one, @two;	}
	if ($n1 == $n3 && $n1 < $n2) {	$frame = 13; push @stops, @one, @three;	}
	if ($n2 == $n3 && $n2 < $n1) {	$frame = 23; push @stops, @two, @three;	}
	if ($n1 == $n2 && $n1 == $n3) {	$frame = 123; push @stops, @one, @two, @three;	}
	
	#return best frame and positions of any stop codons it may contain.
	return  ($frame, $frameDisplay, @stops);
}

##################
sub whichFRAMEMT
		#determines the best reading frame and returns it with any stops within it.
##################
{	my ($thisSeq ) = @_;
	my (@one, @two, @three, @stops) = ();
	
	##print "\nSUB:TEST: length: " , length $thisSeq , "\n";
	#Count stop codons.
	for ($j = 0; $j < (length $thisSeq)-2 ; $j += 3)
	{
		$thisCodon = substr $thisSeq, $j, 3;
		##print $thisCodon;
		if ($thisCodon =~ /TAA|TAG/i)
		{	push @one, ($j+1);
		}
		$thisCodon = substr $thisSeq, ($j+1), 3;
		if ($thisCodon =~ /TAA|TAG/i)
		{	push @two, ($j+2);
		}
		$thisCodon = substr $thisSeq, ($j+2), 3;
		if ($thisCodon =~ /TAA|TAG/i)
		{	push @three, ($j+3);
		}
	}

	$frameDisplay = " one : @one\n two : @two\nthree: @three";####$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#

	#determine which readingframe by least stop codons.
	$n1 = scalar @one;	$n2 = scalar @two;	$n3 = scalar @three;
	if ($n1 < $n2 && $n1 < $n3) {	$frame = 1;	push @stops, @one;}
	if ($n2 < $n1 && $n2 < $n3) {	$frame = 2;	push @stops, @two;}
	if ($n3 < $n1 && $n3 < $n2) {	$frame = 3;	push @stops, @three;}
	if ($n1 == $n2 && $n1 < $n3) {	$frame = 12; push @stops, @one, @two;	}
	if ($n1 == $n3 && $n1 < $n2) {	$frame = 13; push @stops, @one, @three;	}
	if ($n2 == $n3 && $n2 < $n1) {	$frame = 23; push @stops, @two, @three;	}
	if ($n1 == $n2 && $n1 == $n3) {	$frame = 123; push @stops, @one, @two, @three;	}
	
	#return best frame and positions of any stop codons it may contain.
	return  ($frame, $frameDisplay, @stops);
}



##################
sub getNAMES {
	#determines the geneExonName and geneName from the informationLine
##################
	my ($infoLine) = $_[0];
	if ($infoLine =~ /Gene_((NM_[0-9]+)_?[0-9]*)\s/ ) {
		$geneExonName = $1;
		$geneName = $2;
	} elsif ($infoLine =~ /Gene_(([A-Z]+[0-9]+)_?[0-9]*)\s/ ) {
		$geneExonName = $1;
		$geneName = $2;
	} else {
		$geneExonName = 'NoName';
		$geneName = 'NoName';
	}	
	return ($geneExonName, $geneName);
}


####################
sub returnSubstring
####################
{	my ($seq, $first, $last) = @_;
	if ($last) {
		return substr $seq, $first-1, ($last-$first+1);
	}
	else {
		return substr $seq, $first, 20;
	}
}

#############################
sub populate_geno_table {
#############################
	my ($n_chr, $g1, $g2) = @_;
	
	#Initiate Genotype table (3x3)
	my $table = [[0,0,0],[0,0,0],[0,0,0]];
	
	for (my $i=0; $i < $n_chr; $i+=2) {
		my $ind1 = substr $g1, $i, 2;
		my $ind2 = substr $g2, $i, 2;
		next if $ind1 =~ /-/;
		next if $ind2 =~ /-/;
		
		my $s='';
		#Set $s, first array position
		$s = 0 if $ind1 eq '00';
		$s = 1 if $ind1 eq '01';
		$s = 1 if $ind1 eq '10';
		$s = 2 if $ind1 eq '11';
		
		my $t='';
		#Set $t, second array position
		$t = 0 if $ind2 eq '00';
		$t = 1 if $ind2 eq '01';
		$t = 1 if $ind2 eq '10';
		$t = 2 if $ind2 eq '11';
		
		$table->[$s][$t]++;
	}
	
	return $table;
}

#################
sub calculate_D { #from Walsh textbook.
#################
	my ($geno) = @_;
	
	my $N = 0;
	foreach (@{$geno}) {foreach (@{$_}) {$N += $_;} }
	
	my $nA1 = 0;
	foreach (@{$geno}) {
		$nA1 += $_->[0];
		$nA1 += $_->[0];
		$nA1 += $_->[1];
	}
	my $pA1 = $nA1 / (2*$N);

	my $nB1 = 0;
	foreach (@{$geno->[0]}) {
		$nB1 += $_;
		$nB1 += $_;
	}
	foreach (@{$geno->[1]}) {
		$nB1 += $_;
	}
	my $pB1 = $nB1 / (2*$N);

	my $D = $N/($N-1)*(((4*$geno->[0][0] + 2*($geno->[1][0]+$geno->[0][1]) + $geno->[1][1]) / (2*$N)) - 2*$pA1*$pB1);
	my $var_D = $pA1*(1-$pA1)*$pB1*(1-$pB1)/($N-1) + ((2*$pA1-1) * (2*$pB1-1) * $D) / (2*$N) + $D**2 / ($N*($N-1));
	my $se_D = $var_D**0.5;
	return ($D, $var_D, $se_D);
}



########################
sub genotype_LD_chi2 { #from my own derivation of expected genotype frequences under linkage equilibrium and Hardy-Weinberg equilibrium.
#THIS IS WRONG!
######################## It uses a chi2 goodness of fit test.
	my $geno = $_[0];
	my $N = 0;
	foreach (@{$geno}) {foreach (@{$_}) {$N += $_;} } #Count total counted genotypes. (Chromosomes = 2N)
	#Determine allele freq locus 1
	my $nP1 = 0;
	foreach (@{$geno}) {
		$nP1 += 2 * $_->[0];
		$nP1 += $_->[1];
	}
	my $p1 = $nP1 / (2*$N);
	my $p2 = 1- $p1;
	
	my $nQ1 = 0;
	foreach (@{$geno->[0]}) {
		$nQ1 += 2 * $_;
	}
	foreach (@{$geno->[1]}) {
		$nQ1 += $_;
	}
	my $q1 = $nQ1 / (2*$N);
	my $q2 = 1- $q1;
	
	#Compute expected genotype frequencies.
	my $exp_freq = [[$p1**2*$q1**2 , 2*$p1*$p2*$q1**2 , $p2**2*$q1**2]
				   ,[2*$p1**2*$q1*$q2 , 4*$p1*$p2*$q1*$q2 , 2*$p2**2*$q1*$q2]
				   ,[$p1**2*$q2**2 , 2*$p1*$p2*$q2**2 , $p2**2*$q2**2]];
		
	my @exp = ();
	foreach (@{$exp_freq}) {
		my @row =  map {$_ * $N} @{$_}  ;
		push @exp , \@row;
	}
	
	my $chi2 = 0;
	for my $i (0..2) {
		for my $j (0..2) {
			$chi2 += ($geno->[$i][$j] - $exp[$i]->[$j])**2 / $exp[$i]->[$j];
		}
	}
	my $df = 6; #For this test, degrees of freedom are 6.  (9 observed classes - 2 est params[p1,q1] - 1)
	return ($df , $chi2);
}


sub least_squares_regression {
	my @y = @{$_[0]};
	my @x = @{$_[1]};
	
	#sums and sums of squares, ...etc
	my $N = scalar @x;
	my $sum_x = 0;
	$sum_x += $_ foreach @x;
	
	my $sum_y = 0;
	$sum_y += $_ foreach @y;
	
	my $sum_xy;
	for (0..$N-1) {
		$sum_xy += $x[$_] * $y[$_];
	}
	
	my $sum_x2 = 0;
	$sum_x2 += $_**2 foreach @x;
	
	my $sum_y2 = 0;
	$sum_y2 += $_**2 foreach @y;
	
	#calculate linear regression with least squares fit
	my $q = ( $sum_y/$sum_x - $sum_xy/$sum_x2) / ( $N/$sum_x - $sum_x/$sum_x2 ); #intercept
	my $m = ( $sum_y - $N * $q) / $sum_x; #slope

	my @yc = map($q + $m * $_ , @x); #predicted Y values from correlation line
	my $sum_d2 = 0;	#sum of squared difference between y and yc
	for (0..$N-1) {
		$sum_d2 += ($y[$_] - $yc[$_])**2;
	}
	
	my $Sy = ($sum_d2 / $N)**0.5; #Standard error of the estimate
	my $sigma_y = ( $sum_y2 / $N - ($sum_y/$N)**2 )**0.5; #standard deviation of y
	my $r = (1 - $Sy**2 / $sigma_y**2 )**0.5; #correlation coefficient 'r'
	$r *= -1 if $m < 0; #give proper sign to correlation coefficient
	my $rms = 'unk'; #unknown
	my $r2 = $r**2;
	
	return ($q, $m, $r, $rms, $r2);
}


# From example11-7.pl

# call_stride
#
# -given a PDB filename, return the output from the "stride"
#     secondary structure prediction program

sub call_stride {

    use strict;
    use warnings;

    my($filename) = @_;

    # The stride program options
    my($stride) = '/usr/local/bin/stride';
    my($options) = '';
    my(@results) = (  );

    # Check for presence of PDB file
    unless ( -e $filename ) {
        print "File \"$filename\" doesn\'t seem to exist!\n";
        exit;
    }

    # Start up the program, capture and return the output
    @results = `$stride $options $filename`;

    return @results;
}



# From Chapter 8

#
# codon2aa
#
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }elsif ($codon =~ /N|-/) {
    	print STDERR "Bad codon \"$codon\"! Contains an \"N\"!\n";
    	return "?";
    }else{
        print STDERR "Bad codon \"$codon\"!!\n";
        exit;
    }
}

sub codon2aa_silent {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }elsif ($codon =~ /N|-/) {
    	return "?";
    }else{
        exit;
    }
}

#
# codon2aa
#
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup

sub codon2aaMT {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => 'W',    # tryptophan  Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'M',    # methionine  Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }elsif ($codon =~ /N/) {
#    	print STDERR "Bad codon \"$codon\"! Contains an \"N\"!\n";
    	return "?";
    }else{

            print STDERR "Bad codon \"$codon\"!!\n";
            return "?";
    }
}



# From example6-3.pl

sub countG {
    # return a count of the number of G's in the argument $dna

    # initialize arguments and variables
    my($dna) = @_;

    my($count) = 0;

    # Use the fourth method of counting nucleotides in DNA, as shown in
    # Chapter Four, "Motifs and Loops"
    $count = ( $dna =~ tr/Gg//);

    return $count;
}



# From Chapter 8

# dna2peptide 
#
# A subroutine to translate DNA sequence into a peptide

sub dna2peptide {

    my($dna) = @_;

    use strict;
    use warnings;
    use PerlSubs;     # see Chapter 6 about this module

    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($dna,$i,3) );
    }

    return $protein;
}


sub dna2peptide_silent {

    my($dna) = @_;

    use strict;
    use warnings;
    use PerlSubs;     # see Chapter 6 about this module

    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa_silent( substr($dna,$i,3) );
    }

    return $protein;
}


# dna2peptideMT 
#
# A subroutine to translate DNA sequence into a peptide in the metazoan mitochondrion.

sub dna2peptideMT {

    my($dna) = @_;

    use strict;
    use warnings;
    use PerlSubs;     # see Chapter 6 about this module

    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aaMT( substr($dna,$i,3) );
    }

    return $protein;
}


# From example11-5.pl

# extractSEQRES
#
#-given an scalar containing SEQRES lines,
#    return an array containing the chains of the sequence

sub extractSEQRES {

    use strict;
    use warnings;

    my($seqres) = @_;

    my $lastchain = '';
    my $sequence = '';
    my @results = (  );
    # make array of lines

    my @record = split ( /\n/, $seqres);
    
    foreach my $line (@record) {
        # Chain is in column 12, residues start in column 20
        my ($thischain) = substr($line, 11, 1);
        my($residues)  = substr($line, 19, 52); # add space at end
    
        # Check if a new chain, or continuation of previous chain
        if("$lastchain" eq "") {
            $sequence = $residues;
        }elsif("$thischain" eq "$lastchain") {
            $sequence .= $residues;
    
        # Finish gathering previous chain (unless first record)
        }elsif ( $sequence ) {
            push(@results, $sequence);
            $sequence = $residues;
        }
        $lastchain = $thischain;
    }

    # save last chain
    push(@results, $sequence);
    
    return @results;
}



# From example12-2.pl

# extract_HSP_information
#
# -parse a HSP from a BLAST output alignment section
#        - return array with elements:
#    Expect value
#    Query string
#    Query range 
#    Subject string
#    Subject range

sub extract_HSP_information {

    my($HSP) = @_;
    
    # declare and initialize variables
    my($expect) = '';
    my($query) = '';
    my($query_range) = '';
    my($subject) = '';
    my($subject_range) = '';

    ($expect) = ($HSP =~ /Expect = (\S+)/);

    $query = join ( '' , ($HSP =~ /^Query.*\n/gm) );

    $subject = join ( '' , ($HSP =~ /^Sbjct.*\n/gm) );

    $query_range = join('..', ($query =~ /(\d+).*\D(\d+)/s));

    $subject_range = join('..', ($subject =~ /(\d+).*\D(\d+)/s));

    $query =~ s/[^acgt]//g;

    $subject =~ s/[^acgt]//g;

    return ($expect, $query, $query_range, $subject, $subject_range);
}



# From Chapter 8

# extract_sequence_from_fasta_data
#
# A subroutine to extract FASTA sequence data from an array

sub extract_sequence_from_fasta_data {

    my(@fasta_file_data) = @_;

    use strict;
    use warnings;

    # Declare and initialize variables
    my $sequence = '';

    foreach my $line (@fasta_file_data) {

        # discard blank line
        if ($line =~ /^\s*$/) {
            next;

        # discard comment line
        } elsif($line =~ /^\s*#/) {
            next;

        # discard fasta header line
        } elsif($line =~ /^>/) {
            next;

        # keep line, add to sequence string
        } else {
            $sequence .= $line;
        }
    }

    # remove non-sequence data (in this case, whitespace) from $sequence string
    $sequence =~ s/\s//g;

    return $sequence;
}


sub extract_all_sequences_from_fasta_data {
    my(@fasta_file_data) = @_;
    use warnings;

    # Declare and initialize variables
	my %seq;
    my $sequence = '';
	my $name = '';

    foreach my $line (@fasta_file_data) {

        # discard blank line
        if ($line =~ /^\s*$/) {
            next;

        # discard comment line
        } elsif($line =~ /^\s*#/) {
            next;

        # discard fasta header line
        } elsif($line =~ /^>/) {
        	$sequence =~ s/\s//g;
        	$seq{$name} = $sequence unless $name eq '';
        	$name = '';
        	$sequence = '';
        	
        	$name = $1 if $line =~ /^>(\S+)/;

        # keep line, add to sequence string
        } else {
            $sequence .= $line;
        }
    }
    $sequence =~ s/\s//g;
    $seq{$name} = $sequence;
    return (\%seq);
}


sub extract_many_sequences_from_fasta_data {
    my(@fasta_file_data) = @_;
    use warnings;

    # Declare and initialize variables
    my @sequences = ();
    my $sequence = '';
	my $names;
	
    foreach my $line (@fasta_file_data) {

        # discard blank line
        if ($line =~ /^\s*$/) {
            next;

        # discard comment line
        } elsif($line =~ /^\s*#/) {
            next;

        # discard fasta header line
        } elsif($line =~ /^>/) {
        	$sequence =~ s/\s//g;
        	push @sequences, $sequence if $sequence;
        	$names .= $1 if $line =~ /^>(\S+)/;
        	$names .= ',';
        	$sequence = '';

        # keep line, add to sequence string
        } else {
            $sequence .= $line;
        }
    }
    $sequence =~ s/\s//g;
    push @sequences, $sequence if $sequence;
    return ($names, @sequences);
}


# From example10-5.pl

# get_annotation_and_dna
#
#   - given filehandle to open GenBank library file, get next record

sub get_annotation_and_dna {

    my($record) = @_;

    my($annotation) = '';
    my($dna) = '';

    # Now separate the annotation from the sequence data
    ($annotation, $dna) = ($record =~ /^(LOCUS.*ORIGIN\s*\n)(.*)\/\/\n/s);

    # clean the sequence of any whitespace or / characters 
    #  (the / has to be written \/ in the character class, because
    #   / is a metacharacter, so it must be "escaped" with \)
    $dna =~ s/[\s\/\d]//g;

    return($annotation, $dna)
}



# From Chapter 8

# A Subroutine to Read FASTA Files

# get_file_data
#
# A subroutine to get data from a file given its filename

sub get_file_data {

    my($filename) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}



# From example10-5.pl

# get_next_record
#
#   - given GenBank record, get annotation and DNA

sub get_next_record {

    my($fh) = @_;

    my($offset);
    my($record) = '';
    my($save_input_separator) = $/;

    $/ = "//\n";

    $record = <$fh>;

    $/ = $save_input_separator;

    return $record;
}



# From example11-5.pl

# iub3to1
#
#-change string of 3-character IUB amino acid codes (whitespace separated)
#    into a string of 1-character amino acid codes

sub iub3to1 {

    my($input) = @_;
    
    my %three2one = (
      'ALA' => 'A',
      'VAL' => 'V',
      'LEU' => 'L',
      'ILE' => 'I',
      'PRO' => 'P',
      'TRP' => 'W',
      'PHE' => 'F',
      'MET' => 'M',
      'GLY' => 'G',
      'SER' => 'S',
      'THR' => 'T',
      'TYR' => 'Y',
      'CYS' => 'C',
      'ASN' => 'N',
      'GLN' => 'Q',
      'LYS' => 'K',
      'ARG' => 'R',
      'HIS' => 'H',
      'ASP' => 'D',
      'GLU' => 'E',
    );

    # clean up the input
    $input =~ s/\n/ /g;

    my $seq = '';
    
    # This use of split separates on any contiguous whitespace
    my @code3 = split(' ', $input);

    foreach my $code (@code3) {
        # A little error checking
        if(not defined $three2one{$code}) {
            print "Code $code not defined\n";
            next;
        }
        $seq .= $three2one{$code};
    }
    return $seq;
}



# From example9-1.pl

# Example 9-1 Translate IUB ambiguity codes to regular expressions 
# IUB_to_regexp
#
# A subroutine that, given a sequence with IUB ambiguity codes,
# outputs a translation with IUB codes changed to regular expressions
#
# These are the IUB ambiguity codes
# (Eur. J. Biochem. 150: 1-5, 1985):
# R = G or A
# Y = C or T
# M = A or C
# K = G or T
# S = G or C
# W = A or T
# B = not A (C or G or T)
# D = not C (A or G or T)
# H = not G (A or C or T)
# V = not T (A or C or G)
# N = A or C or G or T 

sub IUB_to_regexp {

    my($iub) = @_;

    my $regular_expression = '';

    my %iub2character_class = (
    
        A => 'A',
        C => 'C',
        G => 'G',
        T => 'T',
        R => '[GA]',
        Y => '[CT]',
        M => '[AC]',
        K => '[GT]',
        S => '[GC]',
        W => '[AT]',
        B => '[CGT]',
        D => '[AGT]',
        H => '[ACT]',
        V => '[ACG]',
        N => '[ACGT]',
    );

    # Remove the ^ signs from the recognition sites
    $iub =~ s/\^//g;

    # Translate each character in the iub sequence
    for ( my $i = 0 ; $i < length($iub) ; ++$i ) {
        $regular_expression
          .= $iub2character_class{substr($iub, $i, 1)};
    }

    return $regular_expression;
}



# From example11-3.pl

# list_recursively
#
#   list the contents of a directory,
#              recursively listing the contents of any subdirectories

sub list_recursively {

    my($directory) = @_;

    my @files = (  );
    
    # Open the directory
    unless(opendir(DIRECTORY, $directory)) {
        print "Cannot open directory $directory!\n";
        exit;
    }
    
    # Read the directory, ignoring special entries "." and ".."
    #
    @files = grep (!/^\.\.?$/, readdir(DIRECTORY));
    
    closedir(DIRECTORY);
    
    # If file, print its name
    # If directory, recursively print its contents

    # Notice that we need to prepend the directory name!
    foreach my $file (@files) {
    
        # If the directory entry is a regular file
        if (-f "$directory/$file") {
    
            print "$directory/$file\n";
        
        # If the directory entry is a subdirectory
        }elsif( -d "$directory/$file") {

            # Here is the recursive call to this subroutine
            list_recursively("$directory/$file");
        }
    }
}



# From example7-3.pl
#   and
# From example7-4.pl

# make_random_DNA
#
# Make a string of random DNA of specified length.
#
# WARNING: make sure you call srand to seed the
#  random number generator before you call this function.

sub make_random_DNA {

    # Collect arguments, declare variables
    my($length) = @_;

    my $dna;

    for (my $i=0 ; $i < $length ; ++$i) {
        $dna .= randomnucleotide(  );
    }

    return $dna;
}



# From example7-3.pl
#   and
# From example7-4.pl

# make_random_DNA_set
#
# Make a set of random DNA
#
#   Accept parameters setting the maximum and minimum length of
#     each string of DNA, and the number of DNA strings to make
#
# WARNING: make sure you call srand to seed the
#  random number generator before you call this function.

sub make_random_DNA_set {

    # Collect arguments, declare variables
    my($minimum_length, $maximum_length, $size_of_set) = @_;

    # length of each DNA fragment
    my $length;
    
    # DNA fragment
    my $dna;

    # set of DNA fragments
    my @set;

    # Create set of random DNA
    for (my $i = 0; $i < $size_of_set ; ++$i) {

        # find a random length between min and max
        $length = randomlength ($minimum_length, $maximum_length);

        # make a random DNA fragment
        $dna = make_random_DNA ( $length );

        # add $dna fragment to @set
        push( @set, $dna );
    }

    return @set;
}



# From example9-3.pl

#
# Find locations of a match of a regular expression in a string
#
# 
# return an array of positions where the regular expression
#  appears in the string
#

sub match_positions {

    my($regexp, $sequence) = @_;

    use strict;

    use PerlSubs;     # see Chapter 6 about this module

    #
    # Declare variables
    #

    my @positions = (  );

    #
    # Determine positions of regular expression matches
    #
    
    while ( $sequence =~ /$regexp/ig ) {

        push ( @positions, pos($sequence) - length($&) + 1);
    }

    return @positions;
}



# From example7-4.pl

# matching_percentage
#
# Subroutine to calculate the percentage of identical bases in two
# equal length DNA sequences

sub matching_percentage {

    my($string1, $string2) = @_;

    # we assume that the strings have the same length
    my($length) = length($string1);
    my($position);
    my($count) = 0;

    for ($position=0; $position < $length ; ++$position) {
        if(substr($string1,$position,1) eq substr($string2,$position,1)) {
            ++$count;
        }
    }

    return $count / $length;
}



# From example7-2.pl

# A subroutine to perform a mutation in a string of DNA
#
# WARNING: make sure you call srand to seed the
#  random number generator before you call this function.

sub mutate {

    my($dna) = @_;

    my(@nucleotides) = ('A', 'C', 'G', 'T');

    # Pick a random position in the DNA
    my($position) = randomposition($dna);

    # Pick a random nucleotide
    my($newbase) = randomnucleotide(@nucleotides);

    # Insert the random nucleotide into the random position in the DNA
    # The substr arguments mean the following:
    #  In the string $dna at position $position change 1 character to
    #  the string in $newbase
    substr($dna,$position,1,$newbase);

    return $dna;
}


# From Chapter 7

# mutate_better
#
# Subroutine to perform a mutation in a string of DNA-version 2, in which
#  it is guaranteed that one base will change on each call
#
# WARNING: make sure you call srand to seed the
#  random number generator before you call this function.

sub mutate_better {

    my($dna) = @_;
    my(@nucleotides) = ('A', 'C', 'G', 'T');

    # Pick a random position in the DNA
    my($position) = randomposition($dna);

    # Pick a random nucleotide
    my($newbase);

    do {
        $newbase = randomnucleotide(@nucleotides);

    # Make sure it's different than the nucleotide we're mutating
    }until ( $newbase ne substr($dna, $position,1) );

    # Insert the random nucleotide into the random position in the DNA
    # The substr arguments mean the following:
    #  In the string $dna at position $position change 1 character to
    #  the string in $newbase
    substr($dna,$position,1,$newbase);

    return $dna;
}


# From example11-4.pl

# Example 11-4   Demonstrate File::Find

sub my_sub {
    -f and (print $File::Find::name, "\n");
}



# From example10-5.pl

# open_file
#
#   - given filename, set filehandle

sub open_file {

    my($filename) = @_;
    my $fh;

    unless(open($fh, $filename)) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh;
}



# From example10-1.pl

# parse1
#
# -parse annotation and sequence from GenBank record

sub parse1 {

    my($annotation, $dna, $filename) = @_;

    # $annotation-reference to array
    # $dna       -reference to scalar
    # $filename  -scalar
    
    # declare and initialize variables
    my $in_sequence = 0; 
    my @GenBankFile = (  );
    
    # Get the GenBank data into an array from a file
    @GenBankFile = get_file_data($filename);
    
    # Extract all the sequence lines
    foreach my $line (@GenBankFile) {

        if( $line =~ /^\/\/\n/ ) { # If $line is end-of-record line //\n,
            last; #break out of the foreach loop.
        } elsif( $in_sequence) { # If we know we're in a sequence,
            $$dna .= $line; # add the current line to $$dna.
        } elsif ( $line =~ /^ORIGIN/ ) { # If $line begins a sequence,
            $in_sequence = 1; # set the $in_sequence flag.
        } else{ # Otherwise
            push( @$annotation, $line); # add the current line to @annotation.
        }
    }
    
    # remove whitespace and line numbers from DNA sequence
    $$dna =~ s/[\s0-9]//g;
}



# From example11-6.pl

# parseATOM
#
# -extract x, y, and z coordinates, serial number and element symbol
#     from PDB ATOM record type
#      Return a hash with key=serial number, value=coordinates in a string

sub parseATOM {

    my($atomrecord) = @_;

    use strict;
    use warnings;
    my %results = (  );

    # Turn the scalar into an array of ATOM lines
    my(@atomrecord) = split(/\n/, $atomrecord);

    foreach my $record (@atomrecord) {
       my $number  = substr($record,  6, 5);  # columns 7-11
       my $x       = substr($record, 30, 8);  # columns 31-38
       my $y       = substr($record, 38, 8);  # columns 39-46
       my $z       = substr($record, 46, 8);  # columns 47-54
       my $element = substr($record, 76, 2);  # columns 77-78

        # $number and $element may have leading spaces: strip them
        $number =~ s/^\s*//;
        $element =~ s/^\s*//;

        # Store information in hash
        $results{$number} = "$x $y $z $element";
    }

    # Return the hash
    return %results;
}



# From example11-5.pl

# parsePDBrecordtypes
#
#-given an array of a PDB file, return a hash with
#    keys   = record type names
#    values = scalar containing lines for that record type 

sub parsePDBrecordtypes {

    my @file = @_;

    use strict;
    use warnings;
    
    my %recordtypes = (  );
    
    foreach my $line (@file) {
    
        # Get the record type name which begins at the
        # start of the line and ends at the first space

        # The pattern (\S+) is returned and saved in $recordtype
        my($recordtype) = ($line =~ /^(\S+)/);
    
        # .= fails if a key is undefined, so we have to
        # test for definition and use either .= or = depending
        if(defined $recordtypes{$recordtype} ) {
            $recordtypes{$recordtype} .= $line;
        }else{
            $recordtypes{$recordtype} = $line;
        }
    }
    
    return %recordtypes;
}



# From example9-2.pl

# Example 9-2 Subroutine to parse a REBASE datafile 
# parseREBASE-Parse REBASE bionet file
#
# A subroutine to return a hash where
#    key   = restriction enzyme name
#    value = whitespace-separated recognition site and regular expression

sub parseREBASE {

    my($rebasefile) = @_;

    use strict;
    use warnings;
    use PerlSubs;     # see Chapter 6 about this module

    # Declare variables
    my @rebasefile = (  );
    my %rebase_hash = (  );
    my $name;
    my $site;
    my $regexp;

    # Read in the REBASE file
    my $rebase_filehandle = open_file($rebasefile);

    while(<$rebase_filehandle>) {

        # Discard header lines
        ( 1 .. /Rich Roberts/ ) and next;

        # Discard blank lines
        /^\s*$/ and next;
    
        # Split the two (or three if includes parenthesized name) fields
        my @fields = split( " ", $_);

        # Get and store the name and the recognition site

        # Remove parenthesized names, for simplicity's sake,
        # by not saving the middle field, if any,
        # just the first and last
        $name = shift @fields;

        $site = pop @fields;

        # Translate the recognition sites to regular expressions
        $regexp = IUB_to_regexp($site);

        # Store the data into the hash
        $rebase_hash{$name} = "$site $regexp";
    }

    # Return the hash containing the reformatted REBASE data
    return %rebase_hash;
}



# From example10-6.pl

# parse_annotation
#
#  given a GenBank annotation, returns a hash  with
#   keys: the field names
#   values: the fields

sub parse_annotation {

    my($annotation) = @_; 
    my(%results) = (  );

    while( $annotation =~ /^[A-Z].*\n(^\s.*\n)*/gm ) {
        my $value = $&;
        (my $key = $value) =~ s/^([A-Z]+).*/$1/s;
        $results{$key} = $value;
    }

    return %results;
}



# From example12-1.pl

# parse_blast
#
# -parse beginning and ending annotation, and alignments,
#     from BLAST output file

sub parse_blast {

    my($beginning_annotation, $ending_annotation, $alignments, $filename) = @_;

    # $beginning_annotation-reference to scalar
    # $ending_annotation   -reference to scalar
    # $alignments          -reference to hash
    # $filename            -scalar
    
    # declare and initialize variables
    my $blast_output_file = '';
    my $alignment_section = '';
    
    # Get the BLAST program output into an array from a file
    $blast_output_file = join( '', get_file_data($filename));

    # Extract the beginning annotation, alignments, and ending annotation
    ($$beginning_annotation, $alignment_section, $$ending_annotation)
    = ($blast_output_file =~ /(.*^ALIGNMENTS\n)(.*)(^  Database:.*)/ms);
    
    # Populate %alignments hash
    # key = ID of hit
    # value = alignment section
    %$alignments = parse_blast_alignment($alignment_section);
}



# From example12-1.pl

# parse_blast_alignment
#
# -parse the alignments from a BLAST output file,
#       return hash with
#       key = ID
#       value = text of alignment

sub parse_blast_alignment {

    my($alignment_section) = @_;
    
    # declare and initialize variables
    my(%alignment_hash) = (  );

    # loop through the scalar containing the BLAST alignments,
    # extracting the ID and the alignment and storing in a hash
    #
    # The regular expression matches a line beginning with >,
    # and containing the ID between the first pair of | characters;
    # followed by any number of lines that don't begin with >

    while($alignment_section =~ /^>.*\n(^(?!>).*\n)+/gm) {
        my($value) = $&;
        my($key) = (split(/\|/, $value)) [1];
        $alignment_hash{$key} = $value;
    }

    return %alignment_hash;
}



# From example12-2.pl

# parse_blast_alignment_HSP
#
# -parse beginning annotation, and HSPs,
#     from BLAST alignment
#     Return an array with first element set to the beginning annotation,
#    and each successive element set to an HSP

sub parse_blast_alignment_HSP {

    my($alignment ) = @_;

    # declare and initialize variables
    my $beginning_annotation  = '';
    my $HSP_section  = '';
    my @HSPs = (  );
    
    # Extract the beginning annotation and HSPs
    ($beginning_annotation, $HSP_section )
        = ($alignment =~ /(.*?)(^ Score =.*)/ms);

    # Store the $beginning_annotation as the first entry in @HSPs
    push(@HSPs, $beginning_annotation);

    # Parse the HSPs, store each HSP as an element in @HSPs
    while($HSP_section =~ /(^ Score =.*\n)(^(?! Score =).*\n)+/gm) {
        push(@HSPs, $&);
    }

    # Return an array with first element = the beginning annotation,
    # and each successive element = an HSP
    return(@HSPs);
}



# From Example 10-8

# parse_features
#
#  extract the features from the FEATURES field of a GenBank record

sub parse_features {

    my($features) = @_;   # entire FEATURES field in a scalar variable

    # Declare and initialize variables
    my(@features) = ();   # used to store the individual features

    # Extract the features
    while( $features =~ /^ {5}\S.*\n(^ {21}\S.*\n)*/gm ) {

        my $feature = $&;
        push(@features, $feature);

    }

    return @features;
}



# From example11-7.pl

# parse_stride
#
#-given stride output, extract the primary sequence and the
#    secondary structure prediction, returning them in a
#    two-element array.

sub parse_stride {

    use strict;
    use warnings;

    my(@stridereport) = @_;
    my($seq) = '';
    my($str) = '';
    my $length;

    # Extract the lines of interest
    my(@seq) = grep(/^SEQ /, @stridereport);

    my(@str) = grep(/^STR /, @stridereport);

    # Process those lines to discard all but the sequence
    #  or structure information
    for (@seq) { $_ = substr($_, 10, 50) }
    for (@str) { $_ = substr($_, 10, 50) }

    # Return the information as an array of two strings
    $seq = join('', @seq);
    $str = join('', @str);

    # Delete unwanted spaces from the ends of the strings.
    # ($seq has no spaces that are wanted, but $str may)
    $seq =~ s/(\s+)$//;

    $length = length($1);

    $str =~ s/\s{$length}$//;

    return( ($seq, $str) );
}



# From Chapter 8

# print_sequence
#
# A subroutine to format and print sequence data 

sub print_sequence {

    my($sequence, $length) = @_;

    use strict;
    use warnings;

    # Print sequence in lines of $length
    for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $length ) {
        print substr($sequence, $pos, $length), "\n";
    }
}



# From example7-2.pl
#   and
# From example7-3.pl
#   and
# From example7-4.pl

# randomelement
#
# randomly select an element from an array
#
# WARNING: make sure you call srand to seed the
#  random number generator before you call this function.

sub randomelement {

    my(@array) = @_;

    return $array[rand @array];
}



# From example7-4.pl
#   and
# From example7-3.pl

# randomlength
#
# A subroutine that will pick a random number from
# $minlength to $maxlength, inclusive.
#
# WARNING: make sure you call srand to seed the
#  random number generator before you call this function.

sub randomlength {

    # Collect arguments, declare variables
    my($minlength, $maxlength) = @_;

    # Calculate and return a random number within the
    #  desired interval.
    # Notice how we need to add one to make the endpoints inclusive,
    #  and how we first subtract, then add back, $minlength to
    #  get the random number in the correct interval.
    return ( int(rand($maxlength - $minlength + 1)) + $minlength );
}



# From example7-2.pl
#   and
# From example7-3.pl
#   and
# From example7-4.pl

# randomnucleotide
#
# Select at random one of the four nucleotides
#
# WARNING: make sure you call srand to seed the
#  random number generator before you call this function.

sub randomnucleotide {

    my(@nucleotides) = ('A', 'C', 'G', 'T');

    # scalar returns the size of an array. 
    # The elements of the array are numbered 0 to size-1
    return randomelement(@nucleotides);
}



# From example7-2.pl

# randomposition
#
# A subroutine to randomly select a position in a string.
#
# WARNING: make sure you call srand to seed the
#  random number generator before you call this function.

sub randomposition {

    my($string) = @_;

    # Notice the "nested" arguments:
    #
    # $string is the argument to length
    # length($string) is the argument to rand
    # rand(length($string))) is the argument to int
    # int(rand(length($string))) is the argument to return
    # But we write it without parentheses, as permitted.
    #
    # rand returns a decimal number between 0 and its argument.
    # int returns the integer portion of a decimal number.
    #
    # The whole expression returns a random number between 0 and length-1,
    #  which is how the positions in a string are numbered in Perl.
    #

    return int rand length $string;
}



# From Chapter 8

# revcom 
#
# A subroutine to compute the reverse complement of DNA sequence

sub revcom {

    my($dna) = @_;

    # First reverse the sequence
    my $revcom = reverse $dna;

    # Next, complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $revcom =~ tr/ACGTacgt/TGCAtgca/;

    return $revcom;
}



# From example10-5.pl

# search_annotation
#
#   - search annotation with regular expression

sub search_annotation {

    my($annotation, $regularexpression) = @_;

    my(@locations) = (  );

    # note the /s modifier-. matches any character including newline
    while( $annotation =~ /$regularexpression/isg ) {
        push( @locations, pos );
    }

    return (@locations);
}



# From example10-5.pl

# search_sequence
#
#   - search sequence with regular expression

sub search_sequence {

    my($sequence, $regularexpression) = @_;

    my(@locations) = (  );

    while( $sequence =~ /$regularexpression/ig ) {
        push( @locations, pos );
    }

    return (@locations);
}

#
# ".pm" files need to end with a statement that evaluates as "true"
#



# From Chapter 8

# translate_frame
#
# A subroutine to translate a frame of DNA

sub translate_frame {

    my($seq, $start, $end) = @_;

    my $protein;

    # To make the subroutine easier to use, you won't need to specify
    #  the end point-it will just go to the end of the sequence
    #  by default.
    unless($end) {
        $end = length($seq);
    }

    # Finally, calculate and return the translation
        return dna2peptide ( substr ( $seq, $start - 1, $end -$start + 1) );
}



1;
