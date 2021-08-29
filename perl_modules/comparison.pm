use warnings;
use strict;
package comparison;

sub get_ginkgo_h{
    my ($ginkgo) = @_;
# read the rising and falls for hmm, save to hmm_h
    my $prev_chr = "NA";
    my $prev_pos;
    my $prev_cn;
    my $prev_e;
    my $ginkgo_h;
    my ($chr, $s, $e, $cn);
    open fh_, "<$ginkgo" or die $!;
    while(<fh_>){
        next if($_ =~ /^CHR/ && $_ =~ /START/);
        #next if($_ =~ /chr\,/);
        my @a = split(/\t/, $_);
        ($chr, $s, $e, $cn) = @a[0 .. 3];
        if($chr ne $prev_chr){
            $prev_pos = 0;
            $prev_cn = 2;
# dealing with the previous chromosome
            if($prev_chr ne "NA"){
                if($prev_cn > 2){
                    $ginkgo_h->{f}->{$prev_chr}->{$prev_pos} = 1;
                }
                elsif($prev_cn < 2){
                    $ginkgo_h->{r}->{$prev_chr}->{$prev_pos} = 1;
                }
            }
        }
        if($cn > $prev_cn){
            $ginkgo_h->{r}->{$chr}->{$s} = 1;
        }
        elsif($cn < $prev_cn){
            $ginkgo_h->{f}->{$chr}->{$s} = 1;
        }
        $prev_chr = $chr;
        $prev_pos = $e;
        $prev_cn = $cn;
    }
    close fh_;
# deal with the last end in the last chromosome
    if($prev_chr ne "NA"){
        if($prev_cn > 2){
            $ginkgo_h->{f}->{$prev_chr}->{$prev_pos} = 1;
        }
        elsif($prev_cn < 2){
            $ginkgo_h->{r}->{$prev_chr}->{$prev_pos} = 1;
        }
    }
    return $ginkgo_h;
}

sub get_hmm_h{
    my ($hmm) = @_;
# read the rising and falls for hmm, save to hmm_h
    my $prev_chr = "NA";
    my $prev_pos;
    my $prev_cn;
    my $prev_e;
    my $hmm_h;
    my ($chr, $s, $e, $cn);
    open fh_, "<$hmm" or die $!;
    while(<fh_>){
        next if($_ =~ /chr\,/);
        my @a = split(/\,/, $_);
        ($chr, $s, $e, $cn) = @a[0 .. 3];
        if($chr ne $prev_chr){
            $prev_pos = 0;
            $prev_cn = 2;
# dealing with the previous chromosome
            if($prev_chr ne "NA"){
                if($prev_cn > 2){
                    $hmm_h->{f}->{$prev_chr}->{$prev_pos} = 1;
                }
                elsif($prev_cn < 2){
                    $hmm_h->{r}->{$prev_chr}->{$prev_pos} = 1;
                }
            }
        }
        if($cn > $prev_cn){
            $hmm_h->{r}->{$chr}->{$s} = 1;
        }
        elsif($cn < $prev_cn){
            $hmm_h->{f}->{$chr}->{$s} = 1;
        }
        $prev_chr = $chr;
        $prev_pos = $e;
        $prev_cn = $cn;
    }
    close fh_;
# deal with the last end in the last chromosome
    if($prev_chr ne "NA"){
        if($prev_cn > 2){
            $hmm_h->{f}->{$prev_chr}->{$prev_pos} = 1;
        }
        elsif($prev_cn < 2){
            $hmm_h->{r}->{$prev_chr}->{$prev_pos} = 1;
        }
    }
    return $hmm_h;
}

# updated 08292021, for a list of CNs of ground truth without any 2's in the list. chr, start, end, cn, each line is a CN.
sub get_gt_h_clean{
    my ($gt) = @_;
    my $gt_h;
    my $stat;
    open fh_, "<$gt" or die $!;
    while(<fh_>){
        chomp;
        my @a = split(/\t/, $_);
        if($a[2] > 2){
            $stat = "r";
        }
        elsif($a[2] < 2){
            $stat = "f"
        }
        $gt_h->{$stat}->{$a[0]}->{$a[1]} = 1;
        $gt_h->{$stat}->{$a[0]}->{$a[2]} = 1;
    }
    close fh_;
    return $gt_h;
}


sub get_gt_h{
    my ($gt) = @_;
# read the rising and falling for ground truth, save to gt_h
    my $gt_h;
    my $prev_chr = "NA";
    my $prev_pos;
    my $prev_cn;
    my $prev_e;
    my $gt_tmp;
    open fh_, "<$gt" or die $!;
    while(<fh_>){
        chomp;
        my @a = split(/\t/, $_);
        $gt_tmp->{$a[0]}->{$a[1]}->{cn} = $a[3];
        $gt_tmp->{$a[0]}->{$a[1]}->{e} = $a[2];
    }
    close fh_;
    $prev_chr = "NA";
    foreach my $chr (keys %$gt_tmp){
        my $tag = 0;
        foreach my $pos (sort {$a <=> $b} keys %{$gt_tmp->{$chr}}){
            my $cn = $gt_tmp->{$chr}->{$pos}->{cn};
            my $e = $gt_tmp->{$chr}->{$pos}->{e};
            if($tag == 0){
                $tag = 1;
                if($prev_chr ne "NA"){
                    if($prev_cn > 2){
                        $gt_h->{f}->{$prev_chr}->{$prev_e} = 1;
                    }
                    elsif($prev_cn < 2){
                        $gt_h->{r}->{$prev_chr}->{$prev_e} = 1;
                    }
                }
# the first coming to this chromosome
                if($cn > 2){
                    $gt_h->{r}->{$chr}->{$pos} = 1;
                }
                elsif($cn < 2){
                    $gt_h->{f}->{$chr}->{$pos} = 1;
                }
            }
            elsif($prev_e != $pos){
                $prev_cn = 2;
# dealing with the end before this start (before the gap)
                if($prev_cn > 2){
                    $gt_h->{f}->{$chr}->{$prev_e} = 1;
                }
                elsif($prev_cn < 2){
                    $gt_h->{r}->{$chr}->{$prev_e} = 1;
                }
            }
            else{
# dealing with continuous breakpoints
                if($prev_cn > $cn){
                    $gt_h->{f}->{$chr}->{$pos} = 1;
                }
                elsif($prev_cn < $cn){
                    $gt_h->{r}->{$chr}->{$pos} = 1;
                }
            }
            $prev_e = $e;
            $prev_cn = $cn;
        }
        $prev_chr = $chr;
    }
    return $gt_h;
}

sub get_cn_recall_precision{
    my ($cn, $gt, $t) = @_;
    $DB::single = 1;
    my $cn_h = &get_cn_h($cn);
    my $gt_h = &get_cn_h($gt);
    #my $gt_h = &get_gt_h($gt);
    #my $gt_total = &count_total($gt_h, "r") + &count_total($gt_h, "f");
    my $gt_total = &count_total_cn($gt_h);
    $DB::single = 1;
    my $cn_total = &count_total_cn($cn_h);
    my $common = &get_common_cn($gt_h, $cn_h, $t);
    my $recall = $common / $gt_total;
    my $precision = $common / $cn_total; 
    return ($recall, $precision);
}

sub get_cn_h{
    my ($cn_f) = @_;
    # get copynumber's hash table
    my $prev_chr = "NA";
    my $prev_e;
    my $cn_h;
    open fh_, "<$cn_f" or die $!;
    while(<fh_>){
        chomp;
        my ($chr, $s, $e) = split(/\t/, $_);
        if($chr !~ /^chr/){
            $chr = "chr" . $chr;
        }
        if($chr ne $prev_chr){
            $prev_e = $e;
        }
        else{
# in the same chromosome
            $cn_h->{$chr}->{$prev_e} = 1;
            $prev_e = $e;
        }
        $prev_chr = $chr;
    }
    close fh_;
    return $cn_h;
}

# if gt contains all bins, gt_status=long; otherwise gt_status=short
sub get_recall_precision{
    my ($hmm, $gt, $t, $method, $gt_status) = @_;
    my $hmm_h;
    if(defined $method && $method eq "ginkgo"){
        $hmm_h = &get_ginkgo_h($hmm);
    }
    else{
        $hmm_h = &get_hmm_h($hmm);
    }
    # modified the way to read gt on 08292021 for the short version available
    my $gt_h;
    if($gt_status eq "long"){
        $gt_h = &get_gt_h($gt);
    }
    else{
        $gt_h = &get_gt_h_clean($gt);
    }
        
# compare the two
    my $gt_r_total = &count_total($gt_h, "r");
    my $gt_f_total = &count_total($gt_h, "f");
    my $hmm_r_total = &count_total($hmm_h, "r");
    my $hmm_f_total = &count_total($hmm_h, "f");
    my $common_r = &get_common($gt_h, $hmm_h, "r", $t);
    my $common_f = &get_common($gt_h, $hmm_h, "f", $t);

    my $common = $common_r + $common_f;
    my $recall = $common / ($gt_r_total + $gt_f_total);
    my $precision;
    if($hmm_r_total + $hmm_f_total == 0){
        return ("NA", "NA");
    }
    else{
        $precision = $common / ($hmm_r_total + $hmm_f_total); 
    }
    return ($recall, $precision);

#print "recall is " . $recall . ", precision is " . $precision . "\n";
}

sub count_total_cn{
    my ($h) = @_;
    my $total = 0;
    foreach my $chr (keys %$h){
        foreach my $pos (keys %{$h->{$chr}}){
            $total ++;
        }
    }
    return $total;
}
sub count_total{
    my ($h, $stat) = @_;
    my $total = 0;
    foreach my $chr (keys %{$h->{$stat}}){
        foreach my $pos (keys %{$h->{$stat}->{$chr}}){
            $total ++;
        }
    }
    return $total;
}

sub get_common{
# only count those in h1 that are hit by h2 once, no matter how many h2 are in h1
    my ($h1, $h2, $stat, $t) = @_;
    my $tp = 0;
    foreach my $chr (keys %{$h1->{$stat}}){
        foreach my $p1 (keys %{$h1->{$stat}->{$chr}}){
            foreach my $p2 (keys %{$h2->{$stat}->{$chr}}){
                if(abs($p1 - $p2) < $t){
                    $tp ++;
                    last;
                }
            }
        }
    }
    return $tp;
}

sub get_common_cn{
# only count those in h1 that are hit by h2 once, no matter how many h2 are in h1
    my ($h1, $h2, $t) = @_;
    my $tp = 0;
    foreach my $chr (keys %$h1){
        foreach my $p1 (keys %{$h1->{$chr}}){
            foreach my $p2 (keys %{$h2->{$chr}}){
                if(abs($p1 - $p2) < $t){
                    $tp ++;
                    last;
                }
            }
        }
    }
    return $tp;
}



1;
