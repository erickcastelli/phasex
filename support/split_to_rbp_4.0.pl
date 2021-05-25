# Coded by Erick C. Castelli
#Version 4.0
use strict;
use threads;
use Thread::Semaphore;
use Getopt::Std;

our ($opt_v, $opt_o, $opt_t, $opt_g, $opt_m, $opt_r, $opt_b, $opt_p, $opt_a, $opt_x);
getopts('v:o:t:g:m:r:b:p:a:x');

if ($opt_x) {goto output;}

my $error = 0;

my $bcftools = `which bcftools`;
if ($bcftools eq "") {print "bcftools must be installed and available!\n\n";$error = 1;}

if ($opt_a eq "") 
{
	$opt_a = "/home/lab/Dropbox/gembio-apps/threads/mthreadscmd.pl";
	if (! -e $opt_a) 
	{
		print "mthreadscmd.pl must be available! Indicate it with -a\n\n";$error = 1;
	}
}


if ($opt_g eq "") {$opt_g = "./GenomeAnalysisTK-3.8.1.jar";}

if (! -e $opt_v) {print "VCF file not indicated!\n";$error = 1;}
if ($opt_o eq "") {print "Output folder not indicated!\n";$error = 1;}
if ($opt_t eq "") {$opt_t = 6;print "Threads not indicated. Default: 6\n";}
if (! -e $opt_g) {print "ReadBackedPhasing jar file not indicated or not detected!\n";$error = 1;}
if ($opt_m eq "") {$opt_m = 4;print "Memory per thread not indicated. Default: 2\n";}
if (! -e $opt_r) {print "Reference file not indicated!\n";$error = 1;}
if (! -e $opt_b) {print "BAM folder not indicated!\n";$error = 1;}
if ($opt_p eq "") {print "ReadBackedPhasing parameters not indicated. Default: --phaseQualityThresh 2000\n";}
if ($opt_p eq "") {$opt_p = "--phaseQualityThresh 2000";}


if ($error eq 1) {print "\n";help();}

my $outfolder = $opt_o;
my $threads = $opt_t;
my $rbp = $opt_g;
my $memory = $opt_m;
my $reference = $opt_r;
my $bams = $opt_b . "/";

my @bams = <$bams/*.bam>;
if (scalar(@bams) eq 0) {@bams = <$bams/*/*.bam>;}

my $cmd = "mkdir $outfolder";
system ($cmd);
my $cmd = "mkdir $outfolder/split-vcf";
system ($cmd);
my $cmd = "mkdir $outfolder/rbp";
system ($cmd);
my $cmd = "mkdir $outfolder/cmd";
system ($cmd);
my $cmd = "mkdir $outfolder/output";
system ($cmd);


my $list = `bcftools query -l $opt_v`;
my @list = split("\n",$list);

open (OUT, ">$outfolder/cmd/split.txt");
foreach (@list)
{
	my $cmd = "bcftools view -Ov -s $_ -o $outfolder/split-vcf/$_.vcf $opt_v";
	print OUT $cmd . "\n";
#	system($cmd);
	#last;
}
close OUT;

open (OUT, ">$outfolder/cmd/rbp.txt");
foreach (@list)
{
	my $cmd = "java -Xmx$memory\g -jar $rbp -T ReadBackedPhasing -U -R $reference";
	my $sample = $_;	
	foreach (@bams)
	{
		if ($_ !~ /$sample/){next;}		
		$cmd .= " -I $_";
	}
	$cmd .= " --variant $outfolder/split-vcf/$_.vcf -o $outfolder/rbp/$_.rbp.vcf $opt_p" ;
	print OUT $cmd . "\n";
	#last;
}
close (OUT);

my $cmd = "perl $opt_a -c $outfolder/cmd/split.txt -t $threads";
system($cmd);
my $cmd = "perl $opt_a -c $outfolder/cmd/rbp.txt -t $threads";
system($cmd);


output:
if ($opt_o eq "") {print "Output folder not indicated!\n";$error = 1;}
my $outfolder = $opt_o;

my @head = ();
my %db_snp = ();
my %db = ();
my %samples = ();

my $format_target = "GT:AD:DP:GQ:HP:PQ:PL";
my $input_folder = "$outfolder/rbp/";

my @files = <$input_folder/*.vcf>;
my @formats = split(":",$format_target);

my @snp_list = ();

open (IN, $files[0]);
while (<IN>)
{
	chomp $_;
	my $line = $_;
	if (substr($_,0,2) eq "##")
	{
		my $line = $_;
		if ($_ =~ /fileformat/) {push(@head,$_);next;}
		if ($_ =~ /##contig/) {push(@head,$_);next;}
		if ($_ =~ /##reference/) {push(@head,$_);next;}
		if ($_ =~ /##FORMAT/) 
		{
			foreach (@formats)
			{
				if ($line =~ /ID=$_/){push(@head,$line);last;}
			}
		}
		next;
	}
	if (substr($_,0,2) eq "#C") {next;}
	if ($_ eq "") {next;}

	my @data = split("\t",$_);
	my $snp = join("\t",@data[0..4]);
	$snp .= "\t.\t.\t.\t$format_target";
	my $key = $data[1] . $data[3] . $data[4];
	$db_snp{$key} = $snp;
	push(@snp_list,$key);	
}
close IN;



foreach (@files)
{
	open (IN, $_);
	my $sample = "";
	while (<IN>)
	{
		chomp $_;
		my $line = $_;
		if (substr($_,0,2) eq "##") {next;}
		if (substr($_,0,2) eq "#C") 
		{
			my @data = split("\t",$_);
			$sample = $data[9];
			$samples{$sample} = 1;
			next;
		}
		
		my @data = split("\t",$_);
		my $pos = $data[1];
		my $ref = $data[3];
		my $alt = $data[4];
		my $key = $data[1] . $data[3] . $data[4];
		my $format = $data[8];
		my $gen = $data[9];
		my @data = split(":",$gen);
		my @format_list = split(":",$format);
		my %temp = ();
		for (my $a = 0; $a < scalar(@format_list); $a++)
		{
			$temp{$format_list[$a]} = $data[$a];
		}
		my $final = "";
		foreach(@formats)
		{
			if ($temp{$_} eq "") {$temp{$_} = ".";}			
			$final .= ":" . $temp{$_};
		}
		$final = substr($final,1);

		$db{$key}{$sample} = $final;
	}
	close (IN);
}


open (OUT, ">$outfolder/output/output.vcf");

foreach (@head)
{
	print OUT $_ . "\n";
}

print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

my @sample_list = sort keys %samples;
foreach (@sample_list)
{
	print  OUT "\t" . $_;
}
print OUT "\n";



#my @snps = sort {$a <=> $b} keys %db;
foreach (@snp_list)
{
	print OUT $db_snp{$_};
	my $pos = $_;
	foreach (@sample_list)
	{
		if ($db{$pos}{$_} eq "") {$db{$pos}{$_} = "./.:.:.:.:.:.:.";}
		print OUT "\t" . $db{$pos}{$_};
	}
	print OUT "\n";
}

close (OUT);
exit;








sub help
{
	print "-v VCF_file\n";
	print "-o Output_folder\n";
	print "-t Number_of_threads [default: 6]\n";
	print "-g Path_to_GATK_ReadBackedPhasing\n";
	print "-m Memory (in Gb) per threads [default: 2]\n";
	print "-r Reference file (.fa, .fas)\n";
	print "-b BAM folder\n";
	print "-p ReadBackedPhasing parameters [default: --phaseQualityThresh 2000\n";
	print "-a Path_to_mthreadscmd.pl\n\n";
	print "-x only recreate output\n\n";

	exit;
}
