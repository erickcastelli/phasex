use strict;
use threads;
use Thread::Semaphore;
use Getopt::Std;
use FindBin '$Bin';

#***** System Variables ***************************************
my $program = "mthreadscmd";
my $version = "1.0";
my $os = $^O;
#****************************************************************

our ($max_segments);
our ($opt_c, $opt_t);

getopts('c:t:');

my $max_threads = $opt_t;
if ($max_threads eq "") {$max_threads = 4;}
if (! -e $opt_c) {print "\nCommand list not indicated or not detected\n"; help();}

my @cmds = ();
open (IN, $opt_c);
while (<IN>)
{
	chomp $_;
	if ($_ eq ""){next;}
	push(@cmds,$_);
}
close (IN);
$max_segments = scalar(@cmds);

my $sem = Thread::Semaphore->new($max_threads); # max num of threads
my $current_segment = 0;

my @threads = map {
	# request a thread slot, waiting if none are available:
	$current_segment++; 
	$sem->down;
    
	threads->create(\&run_cmd, $cmds[($current_segment-1)], $current_segment)
	
} 1..$max_segments;
$_->join for @threads;





my %data = ();
my $sample_index = 0;






sub run_cmd {
	 	my ($cmd, $run) = @_;
		
		print "Starting command number $run...\n";

		system ($cmd);
		$sem->up;
	}




sub help {
	print "\n\n";
	print "$program (version $version)\n";
	print "    -c [command list in txt file] \n";
	print "    -t [number of threads] (optional, default: 4)\n";
	print "    example: perl $0 -t 6 -c commands.txt\n";

	exit;	




}






