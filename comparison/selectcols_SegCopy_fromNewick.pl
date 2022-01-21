use warnings;
use strict;

if(@ARGV == 0){
    die "This takes a newick file, a SegCopy file, and a chromosome, and output a new SegCopy file so that only the columns with leaf\${number} mentioned in the newick file, and only the rows with the chromosome will be output. \nUsage: $0 <SegCopy> <Newick> <chromosome> <leaf_prefix: leaf>\n";
}

my ($segcopy, $newick, $chr, $leaf_prefix) = @ARGV;

my $newick_str = "";
open fh_, "<$newick" or die $!;
while(<fh_>){
    chomp;
    $newick_str = $_;
    last;
}
close fh_;

my %nums = map{$leaf_prefix.$_ => 1} split(/\D+/, $newick_str);

my ($h, $header, $cells) = &read_segcopy($segcopy, $chr);

# cells contain the names of all leaves in subsequent order. now come up with a new array that have the index of the leaves that are in nums, the index should be corresponding to those in cells.
my @arr;
my $index = 0;
my @headers = split(/\t/, $header);
print join("\t", @headers[0 .. 2]);
foreach my $cell (@$cells){
    if(defined $nums{$cell}){
        push @arr, $index;
        print "\t" . $cell;
    }
    $index += 1;
}
print "\n";

# output $h if the columns are in nums
# print the leaves that are involved in newick tree
foreach my $chr (keys %$h){
    foreach my $pos (sort {$a <=> $b} keys %{$h->{$chr}}){
        print join("\t", $chr, split(/\./, $pos));
        foreach my $index (@arr){
            print "\t" . $h->{$chr}->{$pos}->{$index};
        }
        print "\n";
    }
}
             


1;


sub read_segcopy{
    my ($segcopy, $chr) = @_;
    my $str_;
    my @cells;
    my $h_;
    open fh_, "<$segcopy" or die $!;
    while(<fh_>){
        chomp;
        my @a = split(/\s+/, $_);
        if($_ =~ /^CHR/){
            $str_ = $_;
            @cells = @a[3 .. $#a];
            next;
        }
        if($a[0] eq $chr){
            foreach my $j (3 .. $#a){
                my $j_ = $j - 3;
                $h->{$a[0]}->{join(".", @a[1 .. 2])}->{$j_} = $a[$j];
            }
        }
    }
    close fh_;
    return ($h, $str_, \@cells);
}


