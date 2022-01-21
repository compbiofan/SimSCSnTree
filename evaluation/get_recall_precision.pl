use warnings;
use strict;
use comparison; 

if(@ARGV < 4){
    die "Function: This code takes the Ginkgo, HMMcopy or Copynumber results, and a ground truth file, output the recall and precision. \nUsage: perl $0 <inferred_cn_file> <groundtruth_cn_file> <inferred method> <ground_truth_type: short for those w/o 2, long for those w/ 2 lines> <tolerance of breakpoint in bp>\n";
}

my ($inferred_f, $gt_f, $inferred_method, $gt_type, $tolerance) = @ARGV;

print "inferred method is " . $inferred_method . "\n"; 
my ($recall, $precision) = comparison::get_recall_precision($inferred_f, $gt_f, $tolerance, $inferred_method, $gt_type);

print "recall is " . $recall . ", and precision is " . $precision . "\n";


