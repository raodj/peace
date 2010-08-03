#!/usr/bin/perl
use strict;

my $base_dir = $ARGV[5];
my $data_type = $ARGV[6]; #sanger, 454, illumina,hybrid

my $oriGeneFile = $ARGV[0];
my $numOfEstSanger = $ARGV[1];
my $numOfEst454 = $ARGV[7];
my $numOfEstIllumina = $ARGV[8];
my $errorRate = $ARGV[2];
my $id = $ARGV[3];
my $useQualityFile = $ARGV[4]; #1-use quality file for Cap3, 0-not use.
my $output; #store all the analysis results

my $geneFile = $oriGeneFile;
my $estFile = "estFile_$id.fa";
my $mstFile = "mstFile_$id.fa";
my $estOutFile = "estFile_$id.out";
my $eastConsensusFile = "$estFile.east.contigs";
my $eastSingletonFile = "$estFile.east.singlets";
my $eastNum = "$estFile.east.num"; #number of used ests
my $eastDebugFile = "$estFile.east.debug";

my $cap3ConsensusFile = "$estFile.cap.contigs";
my $cap3SingletonFile = "$estFile.cap.singlets";
my $cap3DebugFile = "$estFile.cap.debug";

my $tgiclConsensusFile = "$estFile.tgicl.contigs";
my $tgiclSingletonFile = "$estFile.tgicl.singlets";

my $velvetConsensusFile = "$estFile.velvet.contigs";

my $miraConsensusFile = "$estFile.mira.contigs";

	chdir "$base_dir/peace";

	my @contigs;
	my $inContigFile;

    #PEACE+EAST
	chdir "$base_dir/peace";
    my $peace_time;
    my $peace_command;
    my $estFileSanger = "$estFile.sanger";
    my $estFile454 = "$estFile.454";
    my $estFileIllumina = "$estFile.illumina";
    
    $peace_command = `python wrapper.py $geneFile $numOfEstSanger $errorRate $mstFile $estFileSanger`;
    $peace_command = `python metasimWrapper.py $geneFile $numOfEst454 454 $mstFile $estFile454`; 
    $peace_command = `python metasimWrapper.py $geneFile $numOfEstIllumina Illumina $mstFile $estFileIllumina`; 
    $peace_command = `python ModifyComment.py $estFile454`;
	open OUTPUTFILE, ">$estFile" or die ("Fail to open the output file $estFile!"); 
    my $line;
	open INPUTFILE, "<$estFileSanger" or die ("Input file $estFileSanger can not be opened!");
    while (defined($line = <INPUTFILE>)) {
        print OUTPUTFILE $line;
    }
	close INPUTFILE;
	open INPUTFILE, "<$estFile454" or die ("Input file $estFile454 can not be opened!");
    while (defined($line = <INPUTFILE>)) {
        print OUTPUTFILE $line;
    }
	close INPUTFILE;
	open INPUTFILE, "<$estFileIllumina" or die ("Input file $estFileIllumina can not be opened!");
    while (defined($line = <INPUTFILE>)) {
        print OUTPUTFILE $line;
    }
	close INPUTFILE;
	close OUTPUTFILE;

	my $time1 = time();
    $peace_command = `mpiexec ./peace --estFile $estFile --output $estOutFile --output-mst-file $mstFile`;
	my $time2 = time();
    $peace_time = $time2 - $time1; #Peace time
    #$peace_command = `mpiexec ./peace --estFile $estFile --output $estOutFile --output-mst-file $mstFile`;
    #if ($peace_command =~ m/real(\s)+((\d)+)m((\d)+\.(\d)+)s/) {
	#	my $minute = &trim($2);
	#	my $second = &trim($4);
	#	$peace_time = int($minute*60 + $second);
	#}
	chdir "$base_dir/ESTAssemblyC++";
	system("mv $base_dir/peace/$geneFile .");
	system("cp -f $base_dir/peace/$estFile .");
	system("mv $base_dir/peace/$mstFile .");
    system("rm $base_dir/peace/$estOutFile");
	$time1 = time();
	#must ouput to a file, or else it will be printed to the screen and be written into the final output file.
	system("./Main $estFile $mstFile $eastConsensusFile $eastSingletonFile $eastNum > $eastDebugFile");
	$time2 = time();
	my @time;
	$time[0] = $time2 - $time1; #EAST time
	
    #CAP3
	chdir "$base_dir/cap3";
	system("cp -f $base_dir/peace/$estFile .");
	if ($useQualityFile == 1) {
		system("perl writeQualityScore.pl $estFile $errorRate");
	}
	$time1 = time();
	system("cap3 $estFile > $cap3DebugFile"); #must ouput to a file, or else it will be printed to the screen and be written into the final output file.
	$time2 = time();
	$time[1] = $time2 - $time1; #CAP3 time 
    
    #TGICL
    chdir "$base_dir/tgicl";
    system("mkdir $id");
	chdir "$base_dir/tgicl/$id";
	system("mv $base_dir/peace/$estFile .");
	$time1 = time();
	system("tgicl $base_dir/tgicl/$id/$estFile > /dev/null 2> /dev/null"); #must ouput to a file, or else it will be printed to the screen and be written into the final output file.
	$time2 = time();
	$time[2] = $time2 - $time1; #tgicl time 
    system("mv asm_1/contigs $base_dir/tgicl/$tgiclConsensusFile");
    system("mv $estFile.singletons $base_dir/tgicl/$tgiclSingletonFile");
	chdir "$base_dir/tgicl";
    system("rm -rf $id");

    #Velvet
    chdir "$base_dir/velvet";
    system("mkdir $id");
	chdir "$base_dir/velvet/$id";
    system("cp $base_dir/peace/$estFileSanger .");
    system("cp $base_dir/peace/$estFile454 .");
    system("cp $base_dir/peace/$estFileIllumina .");

	$time1 = time();
    system("velveth . 21 -fasta -short $base_dir/velvet/$id/$estFileIllumina -long $base_dir/velvet/$id/$estFile454 $base_dir/velvet/$id/$estFileSanger > /dev/null 2>/dev/null"); #must ouput to a file, or else it will be printed to the screen and be written into the final output file.
    #system("velveth . 21 -fasta -short $base_dir/velvet/$id/$estFileIllumina -long $base_dir/velvet/$id/$estFile454 > /dev/null 2>/dev/null"); #must ouput to a file, or else it will be printed to the screen and be written into the final output file.
    system("velvetg . -cov_cutoff auto > /dev/null 2>/dev/null");
	$time2 = time();
	$time[3] = $time2 - $time1; #velvet time 
    system("mv contigs.fa $base_dir/velvet/$velvetConsensusFile");
	chdir "$base_dir/velvet";
    system("rm -rf $id");

    #Mira
    chdir "$base_dir/mira";
    system("mkdir $id");
	chdir "$base_dir/mira/$id";
	system("mv $base_dir/peace/$estFile .");
    system("mv $base_dir/peace/$estFileSanger msd_in.sanger.fasta");
    system("mv $base_dir/peace/$estFile454 msd_in.454.fasta");
    system("mv $base_dir/peace/$estFileIllumina msd_in.solexa.fasta");

	$time1 = time();
    system("mira -fasta -project=msd -DI:lrt=/tmp -job=denovo,est,sanger,454,solexa SANGER_SETTINGS -LR:wqf=no:mxti=no -AS:epoq=no 454_SETTINGS -LR:wqf=no:mxti=no -AS:epoq=no SOLEXA_SETTINGS -LR:wqf=no:ft=fasta -AS:epoq=no > /dev/null 2>/dev/null"); #must ouput to a file, or else it will be printed to the screen and be written into the final output file.
    #system("mira -fasta -project=msd -DI:lrt=/tmp -job=denovo,est,454,solexa 454_SETTINGS -LR:wqf=no:mxti=no -AS:epoq=no SOLEXA_SETTINGS -LR:wqf=no:ft=fasta -AS:epoq=no > /dev/null 2>/dev/null"); #must ouput to a file, or else it will be printed to the screen and be written into the final output file.
	$time2 = time();
	$time[4] = $time2 - $time1; #velvet time 
    system("mv msd_assembly/msd_d_results/msd_out.unpadded.fasta $base_dir/mira/$miraConsensusFile");
	chdir "$base_dir/mira";
    system("rm -rf $id");

	print "$time[1]\t$peace_time\t$time[0]\t$time[2]\t$time[3]\t$time[4]";

	#copy the intermediate files for debugging later
	#system("cp -f $geneFile $base_dir/analysis/.");
	#system("cp -f $eastDebugFile $base_dir/analysis/.");
	#system("cp -f $mstFile $base_dir/analysis/.");
	#system("cp -f $estFile $base_dir/analysis/.");
	#system("cp -f $eastConsensusFile $base_dir/analysis/.");
	#system("cp -f $base_dir/cap3/$cap3ConsensusFile $base_dir/analysis/.");
	
	#remove the intermediate files. Notice: CANNOT remove the original gene file $oriGeneFile, it will be used by other jobs.
	chdir "$base_dir/ESTAssemblyC++";
    system("mv $base_dir/cap3/$cap3ConsensusFile $base_dir/.");
    system("mv $base_dir/cap3/$cap3SingletonFile $base_dir/.");
    system("mv $base_dir/ESTAssemblyC++/$eastConsensusFile $base_dir/");
    system("mv $base_dir/ESTAssemblyC++/$eastSingletonFile $base_dir/");
    system("mv $base_dir/ESTAssemblyC++/$estFile $base_dir/");
    system("mv $base_dir/tgicl/$tgiclConsensusFile $base_dir/");
    system("mv $base_dir/tgicl/$tgiclSingletonFile $base_dir/");
    system("mv $base_dir/velvet/$velvetConsensusFile $base_dir/");
    system("mv $base_dir/mira/$miraConsensusFile $base_dir/");
	system("rm $mstFile");
    system("rm $geneFile");
	system("rm $eastNum");
	system("rm $eastDebugFile");
	system("rm $base_dir/cap3/$estFile*");

# read gene from the input gene file, put the gene in one line.
sub getGene {
	my $fileName = $_[0];
	open INPUTFILE, "<$fileName" or die ("Input file $fileName can not be opened!");
	my $gene;
	my $line;
	
	while (defined($line = <INPUTFILE>)) {
		chomp $line;
		if (!($line =~ />(.)*/i)) {
			#print ("$line\n");
			$gene = $gene . trim($line);
		}
	};

	close INPUTFILE;
	return $gene;
}

# read consensus from the input file, one consensus in one line.
# one parameter: the name of the input file. The file is in fasta format.
# output: an array which includes all the consensus read from the input file in fasta format.
sub getConsensus {
	my $contigFile = $_[0];
	open INPUTFILE, "<$contigFile" or die ("Input file $contigFile can not be opened!");
	my @genes;
	my $contig = "";
	my $line;
	my $index = 0;
	
	while (defined($line = <INPUTFILE>)) {
		chomp $line;
		if ($line =~ />(.)*/i) { #a line of comment
			if (!$contig eq "") {
				$genes[$index++] = $contig;
				$contig = "";
			}
		} else {
			$contig = $contig . trim($line); #put one contig in one line
		}
	};
	if (!$contig eq "") {
		$genes[$index] = $contig;
	}
	
	close INPUTFILE;
	return @genes;
}

# write all the contigs into an output file, one contig in  one line.
# two parameter: output file name, an array which includes all the contigs.
sub writeContigToFile {
	my $fileName = shift @_;
	my @contigs = @_;
	my $size = @contigs;
	open OUTPUTFILE, ">$fileName" or die ("Fail to open the output file $fileName!"); 
	for (my $i=0; $i<$size; $i++) {
		print OUTPUTFILE "$contigs[$i]\n";
	}
	
	close OUTPUTFILE;
}

#delete preceding and trailing white spaces
sub trim
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
