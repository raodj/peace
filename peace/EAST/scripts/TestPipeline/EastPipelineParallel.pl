#!/usr/bin/perl -w
use strict;
use redhawk;
use Cwd;

# Parameters relationg to redhawk
my $job_limit = 200;
my $pause     = 10;

# Parameters relating to checkpoint information
my %job_hash;   # Maps job description to job id

# Paramters relating to script
my $executable = "python ./MultiGenes.py"; 

my $geneSource = "all_zf_cdnas.reduced.fa";
    
# testing parameters
my %testing_params = (
    gene_source => ["all_zf_cdnas.reduced.fa"],
    num_genes => [10],
    min_length => [2000],
    max_length => [20000],
    error_rate => [0,5],  #[(0..10)],
    use_quality => [0],
    coverage => [10,30],    #[map {5*$_} (3..10)]
    data_type => ["sanger"] #sanger, 454, illumina
    );
my $num_trials = 3;

#***********************
# Set user arguments
my $reset_batch = 0;
while (substr($ARGV[0], 0, 1) eq "-") {
    my $switch = shift @ARGV;
    if (grep {$switch eq $_} qw(-r -reset)) {
	$reset_batch = 1;
    }
    else {
	die("Bad switch: $switch");
    }
}

my $batch_file = shift @ARGV;

# Open batch file
my $bfp;
if (!$reset_batch && -e $batch_file) {
    load_batch_file($batch_file);
    open($bfp, ">>$batch_file");
}
else {
    open($bfp, ">$batch_file");
}

recursiveTestCall([sort keys %testing_params], {});

sub load_batch_file {
    my $fp;
    open($fp, $batch_file);
    while (my $line = <$fp>) {   # Load job_hash
	chomp($fp);
	my ($unique_id, $job_file, $job_id, $dir, $output_file) = split(/\s/, $line);
	#if (job_running($job_file, $job_id) || job_done($output_file)) {
    if (job_done($output_file)) {
	    $job_hash{$unique_id} = [$job_file, $job_id, $dir, $output_file];
	}
    }
    close($fp);
}

sub job_done {
    my $output_file = $_[0];
    if (not -e $output_file) {
      return 0;  #File does not exist; return false
    }

    my $line = `tail -1 $output_file`;
    chomp($line);
    return $line eq 'DONE';
}

sub recursiveTestCall {
    my @keys = @{$_[0]};
    my %params = %{$_[1]};

    if (@keys == 0) {
	for my $t ((1..$num_trials)) {
	    $params{trial_num} = $t;
	    callSingleTest(%params);
	}
	delete $params{trial_num};
    }
    else {
	my $key = shift @keys;
	my @values = @{$testing_params{$key}};
	for my $v (@values) {
	    $params{$key} = $v;
	    recursiveTestCall(\@keys, \%params);
	}
	delete $params{$key}
    }
}


sub callSingleTest {
    my %params = @_;
    
    my $gene_source = $params{gene_source};
    my $num_genes = $params{num_genes};
    my $min_len = $params{min_length};
    my $max_len = $params{max_length};
    my $error_rate = $params{error_rate};
    my $coverage = $params{coverage};
    my $use_quality = $params{use_quality};
    my $trial_num = $params{trial_num};
    my $data_type = $params{data_type};

    my $unique_id = make_id(%params);
    my $output_dir = get_directory($unique_id);
    my $output_file = "$output_dir/output.$unique_id";

    unless (exists $job_hash{$unique_id}) {

	my $execute = "$executable $gene_source $num_genes $min_len $max_len $error_rate $coverage $use_quality $output_file $unique_id\_$$ $data_type";
	print $execute, "\n";
	my $job_file = create_batch_file(
	    file => "MultiGenes.$unique_id.$$",
	    executable => "module load python-2.6; module load biopython; module load bioperl; $execute",
	    output_location => $output_dir
	    );
	
	wait_on_job_limit($job_limit, $pause);
	$job_hash{$unique_id} = [$job_file, execute_batch($job_file), $output_dir, $output_file];
	print $bfp "$unique_id ".join(" ", @{$job_hash{$unique_id}})."\n";
    }
}

sub make_id {
    my %params = @_;
    return join("_", ("ng".$params{num_genes}, "e".$params{error_rate}, "c".$params{coverage}, "t".$params{trial_num}));
}

sub get_directory {
    my $dir = "results.d";
    mkdir($dir) unless -e $dir;
    for my $i (split "_", $_[0]) {
	$dir = "$dir/$i.d";
	mkdir($dir) unless -e $dir;
    }

    return $dir;
}
