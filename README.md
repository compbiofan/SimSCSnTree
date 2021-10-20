# SCSim

Note: This repository contains scripts that were developed for "SCSim: a single-cell simulator for multi-purpose intra-tumor heterogeneity benchmark study". 

Authors: Xian Fan (xfan2@fsu.edu), Luay Nakhleh (nakhleh@rice.edu) 

## Table of Contents
- [Installing SCSim.](#install_SCSim)
    * [Software requirements](#software_requirements)
    * [Environment setup](#environment_setup)
- [Usage of SCSim.](#usage_of_single_cell_simulator)
    * [General usage.](#general_usage)
    * [Control of CNA size and rate (#CNA)
    * [Control of whole chromosmoe duplication](#WCD)
    * [Control of SNV rate](#SNV)
    * [Control of tree structure](#tree_structure)
    * [Control of parameters for longitudinal study](#longitudinal)
    * [Control of read depth fluctuation](#read_fluctuation)
- [Examples.](#examples)
    * [Simulating large dataset](#large_dataset)
    * [Simulating reads with different ploidies](#ploidies)
    * [Simulating reads with different levels of fluctuation](#fluctuations)
- [Extract data from npy.](#data_extraction)
    * [Making ground truth from the simulator for comparison](#ground_truth)
    * [Generating a newick formatted tree from .npy file in simulation](#newick)


# <a name="install_SCSim"></a>Installing SCSim.
## <a name="software_requirements"></a>Software requirements ##

1. Python 3 or up.

2. Python modules: numpy, graphviz, anytree. 

## <a name="environment_setup"></a>Environment setup ##

Suppose $this_dir is the path of this package.

1. Add pipeline to your directory. In bash,

    ```if [ -d $this_dir ]; then PATH="$this_dir:$PATH" fi```

2. Make the binary from the revised wgsim. 

    ```cd $this_dir/wgsim-master```

    ```gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm```

    ```chmod u+x wgsim```

3. Python modules: numpy, graphviz, anytree. 
    
    ```pip install numpy```
        
    ```pip install graphviz```
        
    ```pip install anytree```

# <a name="usage_of_single_cell_simulator"></a>Usage of SCSim.
## <a name="general_usage"></a>General usage##

SCSim has two steps. Step 1 generate a tree, each node of which contains a genome and each edge of which CNA(s) and/or SNV(s) are imputed. Step 2 samples reads from the genomes on the nodes selected. The following shows how each step works. 

1. Generate the tree with CNVs/SNVs on the edges. This step generates two npy files, one for the tree (containing all information of CNVs/SNVs on the edges and the tree structure) and one for intermediate files (containing chromosome name, chromosome length and the index of the nodes). The .npy files will be used in step 2 for sampling reads. 
    
    ```python main.par.overlapping.py``` followed by the following parameters grouped by their functions.  

    * Parameters controlling file IO:

        -r (--directory)    Location of simulated data. The program will remove the whole directory if it already exists. Otherwise it will create one. (default: test)
        -t (--template-ref) The reference file to sequence the reads. 
        -o (--outfile)      The standard output file, will be saved in output folder, just give the file name. (default: std.out)

    * Parameters controlling tree structure: -n, -B, -A, -F, -G
    
        -n (--cell-num)     Number of the cells. Always greater than -F treewidth. Treewidth controls the total number of clones whereas cell-num controls the total number of cells sequenced at a certain tree depth. 
        -B (--Beta)         The program uses the Beta-splitting model to generate the phylogenetic tree. Specify a value between [0, 1]. (default: 0.5)
        -A (--Alpha)        The Alpha in Beta-splitting model. Specify a value between [0, 1]. The closer Alpha and Beta, the more balanced the tree. (default: 0.5).
        -F (--treewidth)    The mean of the tree width distribution. The final tree width will be sampled from a Gaussian with this mean and a fixed standard deviation. (default: 8)
        -G (--treedepth)    The mean of the tree depth distribution. The final tree depth will be sampled from a Gaussian with this mean and a fixed standard deviation. (default: 4 counting from the first cancer cell)
        -H (--treewidthsigma)	The standard deviation of the tree width distribution. To get exactly the tree width defined by -F, use a very small standard deviation, e.g., 0.0001. (default: 0.5)
        -K (--treedepthsigma)	The standard deviation of the tree depth distribution. To get exactly the tree depth defined by -F, use a very small standard deviation, e.g., 0.0001. (default: 0.5)
        
    * Parameters controlling CNAs: -c, -d, -m, -e, -a

        -c (--cn-num)       The average number of copy number variations to be added on a branch. (default: 1)
        -d (--del-rate)     The rate of deletion as compared to amplification. (default: 0.5)
        -m (--min-cn-size)  Minimum copy number size. (default: 200,000bp)
        -e (--exp-theta)    The parameter for the Exponential distribution for copy number size, beyond the minimum one. (default: 0.000001)
        -a (--amp-p)        The parameter for the Geometric distribution for the number of copies amplified. (default: 0.5)	
       
    * Parameters controlling whole chromosome duplication on the branch to the root: -X, -W, -C, -E, -J

        -X (--multi-root)   The multiplier of the mean CNV on root. (default: 4)
        -W (--whole-amp)    If there is whole chromosome amplification, 1 as yes. (default: 1) 
        -C (--whole-amp-rate)   Whole amplification rate: rate of an allele chosen to be amplified (default: 0.2)
        -E (--whole-amp-num)    Whole amplification copy number addition, which occurs to one allele at a time. (default: 1)
        -J (--amp-num-geo-par)  Whole amplification copy number distribution (geometric distribution parameter: the smaller, the more evenly distributed). (default: 1)
    
    * Parameters controlling SNVS: -R

        -R (--snv-rate)     The rate of the snv. snv-rate * branch-length = # snvs. (default: 1)
        
    * Parameters controlling bulk sequencing:

        -U (--bulk-levels)	The levels of the bulk sequencing separated by semicolon. The definition of the levels is the same as in -L. The default for this option is NA, meaning no bulk sequencing. 
        
    * Parameters for longitudinal study:

        -L (--levels)	This is for both tree inference and longitidunal study. For multiple levels, use semicolon to separate them. The first tumor cell has level 1. If counting from the bottom (leaf) of the tree, use minus before the number. For example, -1 is the leaf level. The range of the level should be within [-depth, depth]. Users can specify desired levels according to -G to know which levels are available. If that is the case, use a very small -K to make sure the depth is not smaller than the biggest level you specify. (default: -1) 
        -U (--bulk-levels)	The levels of the bulk sequencing separated by semicolon. The definition of the levels is the same as in -L. The default for this option is NA, meaning no bulk sequencing. 
        
2. Sample reads from specified genomes. 

    ```python main.par.overlapping.py -k 1``` followed by the following parameters grouped by their functions.  
    
    * Parameter controlling which step to run:

        -k (--skip-first-step)  If the alleles for all nodes have been made, the step can be skipped. Make it 1 then. (default: 0)
        
    * Parameters for IO:

        -f (--fa-prefix)    The prefix of the alleles and read names. (default: ref)
        
    * Parameter specifying which read simulator to use:

        -S (--wgsim-dir)    The directory of the binary of wgsim. It is in the same folder of this main.py. (need to specify) 
        
    * Parameters controlling read depth fluctuation:

        -x (--Lorenz-x)     The value on the x-axis of the point furthest from the diagonal on the Lorenz curve imitating the real coverage uneveness. (default: 0.5) 
        -y (--Lorenz-y)     The value on the y-axis of the Lorenz curve imitating the real coverage unevenness. x > y. The closer (x, y) to the diagonal, the better the coverage evenness. (default: 0.4) 
        -v (--coverage)     The average coverage of the sequence. (default: 0.02)
        -l (--readlen)      Read length for each read sequenced. (default: 35bp)
        -w (--window-size)  Within a window, the coverage is according to a Gaussian distribution. Neighboring windows' read coverage is according to a Metropolis Hasting process. (default: 200000bp)
        -u (--acceptance-rate)  The probability to accept a proposal in Metropolis Hasting. (default: 0.5)
    
    * Parameters controlling bulk sequencing:
    
        -V (--cov-bulk)	The coverage of the bulk sequencing. The same for all levels. This parameter is needed when -U is identified. (default: 30) 
    
    * Parameters for parallel job submissions:

        -p (--processors)   Numbers of processors available.
        -Y (--leaf-index-range) For parallele job submission. >= min, < max leaf index will be processed. min.max. This counts leaf nodes from 0. (default: -1)
        -I (--leaf-ID-range) For parallele job submission. >= min, < max leaf ID will be processed. min.max. This counts from root 0-based so that internal nodes sequencing can be parallelized. (default: -1). When both -Y and -I are -1, all leaves will be processed.
        
    * Parameters for clonality study:  
          	
        -M (--single-cell-per-node)	If this is on, each node represents one cell and there is no clonality in the node. In this case tree_width will be the same as n (leaf num). 1 is on. (default: 0)
    
    For a complete list of options, type
        
    ```python main.par.overlapping.py --help``` 
    
## <a name="CNA"></a>Control of CNA size and rate.

On a branch, the number of the CNA imputed follows a Poisson distribution, the mean of which follows an exponential distribution with p specified by -c (--cn-num). The deletion rate as compared to copy number gain follows a binomial distribution with p specified by -d (--del-rate). The CNA size follows an exponential distribution with p specified by -e (--exp-theta) plus a minimum CNA size specified by -m (--min-cn-size). If it is a copy number gain, the numbers of gain follows a Geometric distribution with p specified by -a (--amp-p). 

        -c (--cn-num)       The average number of copy number variations to be added on a branch. (default: 1)
        -d (--del-rate)     The rate of deletion as compared to amplification. (default: 0.5)
        -m (--min-cn-size)  Minimum copy number size. (default: 200,000bp)
        -e (--exp-theta)    The parameter for the Exponential distribution for copy number size, beyond the minimum one. (default: 0.000001)
        -a (--amp-p)        The parameter for the Geometric distribution for the number of copies amplified. (default: 0.5)	

## <a name="WCD"></a>Control of whole chromosome duplication.

The whole chromosome duplications are imputed in the trunk branch connecting the normal cell and the first tumor cell if -W (--whole-amp) is 1. For each chromosome, the probability that it is amplified equals -C (--whole-amp-rate) and the number of copies amplified follows a geometric distribution with p specified by -J (--amp-num-geo-par) multiplied by a number specified by -E (--whole-amp-num). 

        -X (--multi-root)   The multiplier of the mean CNV on root. (default: 4)
        -W (--whole-amp)    If there is whole chromosome amplification, 1 as yes. (default: 1) 
        -C (--whole-amp-rate)   Whole amplification rate: rate of an allele chosen to be amplified (default: 0.2)
        -E (--whole-amp-num)    Whole amplification copy number addition, which occurs to one allele at a time. (default: 1)
        -J (--amp-num-geo-par)  Whole amplification copy number distribution (geometric distribution parameter: the smaller, the more evenly distributed). (default: 1)

## <a name="SNV"></a>Control of SNV rate.

On a branch, the number of the SNV imputed follows a Poisson distribution, the mean of which equals to snv-rate (specified by -R) multiplied by branch length which is sampled from an exponential distribution with p=1.

         -R (--snv-rate)     The rate of the snv. snv-rate * branch-length = # snvs. (default: 1)

## <a name="tree_structure"></a>Control of tree structure.

The binary tree's branch splitting follows Beta-splitting model so that the splitting of the cells between the left and right branches for each split follows a Beta distribution, whose alpha and beta parameters are specified by -A (--Alpha) and -B (--Beta). 

The splitting ends when the number of cells / subclones on the leaf level reaches -n (--cell-num). The tree witdth and depth can be controlled by Gaussian distributions, the mean and standard deviation of which are specified by -F (--treewidth), -H (--treewidthsigma), -G (--treedpeth) and -K (--treedepthsigma). 

        -n (--cell-num)     Number of the cells. Always greater than -F treewidth. Treewidth controls the total number of clones whereas cell-num controls the total number of cells sequenced at a certain tree depth. 
        -B (--Beta)         The program uses the Beta-splitting model to generate the phylogenetic tree. Specify a value between [0, 1]. (default: 0.5)
        -A (--Alpha)        The Alpha in Beta-splitting model. Specify a value between [0, 1]. The closer Alpha and Beta, the more balanced the tree. (default: 0.5).
        -F (--treewidth)    The mean of the tree width distribution. The final tree width will be sampled from a Gaussian with this mean and a fixed standard deviation. (default: 8)
        -G (--treedepth)    The mean of the tree depth distribution. The final tree depth will be sampled from a Gaussian with this mean and a fixed standard deviation. (default: 4 counting from the first cancer cell)
        -H (--treewidthsigma)	The standard deviation of the tree width distribution. To get exactly the tree width defined by -F, use a very small standard deviation, e.g., 0.0001. (default: 0.5)
        -K (--treedepthsigma)	The standard deviation of the tree depth distribution. To get exactly the tree depth defined by -F, use a very small standard deviation, e.g., 0.0001. (default: 0.5)
        
## <a name="longitudinal"></a>Control of parameters for longitudinal study:

        -L (--levels)	This is for both tree inference and longitidunal study. For multiple levels, use semicolon to separate them. The first tumor cell has level 1. If counting from the bottom (leaf) of the tree, use minus before the number. For example, -1 is the leaf level. The range of the level should be within [-depth, depth]. Users can specify desired levels according to -G to know which levels are available. If that is the case, use a very small -K to make sure the depth is not smaller than the biggest level you specify. (default: -1) 
        -U (--bulk-levels)	The levels of the bulk sequencing separated by semicolon. The definition of the levels is the same as in -L. The default for this option is NA, meaning no bulk sequencing. 

## <a name="read_fluctuation"></a>Control of read depth fluctuation.



        -x (--Lorenz-x)     The value on the x-axis of the point furthest from the diagonal on the Lorenz curve imitating the real coverage uneveness. (default: 0.5) 
        -y (--Lorenz-y)     The value on the y-axis of the Lorenz curve imitating the real coverage unevenness. x > y. The closer (x, y) to the diagonal, the better the coverage evenness. (default: 0.4) 
        -v (--coverage)     The average coverage of the sequence. (default: 0.02)
        -l (--readlen)      Read length for each read sequenced. (default: 35bp)
        -w (--window-size)  Within a window, the coverage is according to a Gaussian distribution. Neighboring windows' read coverage is according to a Metropolis Hasting process. (default: 200000bp)
        -u (--acceptance-rate)  The probability to accept a proposal in Metropolis Hasting. (default: 0.5)

    
# <a name="examples"></a>Examples. 

## <a name="large_dataset"></a>Simulating large dataset. 
##  

The following lists the command to simulate the large dataset. Step 2 of the simulator is the same as the general one described in "Usage". 

```python main.par.py -S $wgsim-master -r $dir -n 10000 -p 1 -X 8 -t $ref -W 1 -C 0.3 -E 1 -l 36 -m 2000000 -e 5000000```

## <a name="ploidies"></a>Simulating reads with different ploidies. 
## 

The following lists the command to simulate the tree and the alternative alleles (step 1 of the simulator) for different ploidies. Step 2 of the simulator is the same as the general one described in "Usage". 

* Ploidy 1.55

```python main.par.py -S $wgsim-master -r $dir -n 100 -p 1 -X 25 -t $ref -W 0 -l 36 -m 2000000 -e 5000000 -d 1 -c 3```

* Ploidy 2.1

```python main.par.py -S $wgsim-master -r $dir -n 100 -p 1 -X 8 -t $ref -W 1 -C 0.05 -l 36 -m 2000000 -e 5000000```

* Ploidy 3.0

```python main.par.py -S $wgsim-master -r $dir -n 100 -p 1 -X 8 -t $ref -W 1 -C 0.5 -l 36 -m 2000000 -e 5000000 -E 1```

* Ploidy 3.8

```python main.par.py -S $wgsim-master -r $dir -n 100 -p 1 -X 8 -t $ref -W 1 -C 0.9 -l 36 -e 5000000 -E 1 -m 10000000```

* Ploidy 5.26

```python main.par.py -S $wgsim-master -r $dir -n 100 -p 1 -X 8 -t $ref -W 1 -C 0.9 -l 36 -e 5000000 -E 1 -m 10000000 -J 0.55```

## <a name="fluctuations"></a>Simulating reads with different levels of fluctuation ##

The following lists the command to simulate the reads  (step 2 of the simulator) for different fluctuations. 

* MALBAC

```python main.par.py -S $wgsim-master -r $dir -l 36 -x 0.5 -y 0.27 -k 1```

* DOP-PCR

```python main.par.py -S $wgsim-master -r $dir -l 36 -x 0.5 -y 0.28 -k 1```

* TnBC

```python main.par.py -S $wgsim-master -r $dir -l 36 -x 0.5 -y 0.33 -k 1```

* Bulk

```python main.par.py -S $wgsim-master -r $dir -l 36 -x 0.5 -y 0.38 -k 1```

# <a name="commands_3methods"></a>Commands to run HMMcopy, Ginkgo and CopyNumber. #

## <a name="hmmcopy"></a>HMMcopy ##

### Installing the software and preparing the files. These steps are general for all datasets, but they should be done only once. ###

1. Download new version of HMMcopy from https://github.com/shahcompbio/single_cell_pipeline/tree/master/single_cell/workflows/hmmcopy. Install according to the instruction. 

2. Download scripts to preprocess files for HMMcopy from https://shahlab.ca/projects/hmmcopy_utils. Make the binaries according to the instruction.

3. Prepare mappability file. 

    * Download the bigWig file from http://genome.ucsc.edu/cgi-bin/ hgFileUi?db=hg19&g=wgEncodeMapability. Choose kmer size 36. Note this is for hg19 only, i.e., bams generated from simulation or real data should be in hg19 as well. The downloaded file's name is wgEncodeCrgMapabilityAlign36mer.bigWig. 

    * Make the mappability file.

        ```mapCounter -w 200000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY wgEncodeCrgMapabilityAlign36mer.bigWig > $ref.mp.seg```

        mapCounter can be found in hmmcopy_utils/bin. 

        $ref: the reference fa file in absolute path.

4. Prepare the gc file. 

    ```gcCounter -w 200000 $ref > $ref.gc.seg```

    $ref: the reference fa file in absolute path. 

    gcCounter can be found in hmmcopy_utils/bin. 

### Preparing the read count file and run HMMcopy. These steps are specific to a given bam. ###

1. Generate the read count file given a bam file.
    * Make the raw read count as a wig file ($bam.wig). 
    
        ```readCounter -w $window_size $bam > $bam.wig```

        $bam: the bam file which is the input. 

        readCounter can be found in hmmcopy_utils/bin.

    * Correct read count by mappability and gc files. This step generates a csv file ($bam.wig.csv).
 
        ```python correct_read_count.py $ref.gc.seg $ref.mp.seg $bam.wig $bam.wig.csv```

        correct_read_count.py can be found in the new HMMcopy package.  

2. Run HMMcopy. 

    ```Rscript hmmcopy.R --corrected_data=$bam.wig.csv --outdir=$dir --sample_id=$id --param_multiplier="1,2,3,4,5,6" --param_str=$s --param_e=$e --param_mu="0,1,2,3,4,5,6,7,8,9,10,11" --param_l=20 --param_nu=$nu --param_k="100,100,700,100,25,25,25,25,25,25,25,25" --param_m="0,1,2,3,4,5,6,7,8,9,10,11" --param_eta=50000 --param_g=3 --param_s=1```

    Rscript version should be R-3.5.1 or above. We tried lower versions of R and found they do not work. 

    hmmcopy.R can be found in the new HMMcopy package. 

    $dir is the output directory for this cell. 

    Input is $bam.wig.csv. Output file containing segmentation and absolute copy number is in $dir/0/segs.csv. Predicted ploidy can be found in the second row, first column in $dir/0/metrics.csv.

    The explanation of the parameters can be found in https://rdrr.io/bioc/HMMcopy/man/HMMsegment.html. We suggest using a large $s (e.g., 10,000,000). For $e and $nu, we found the best combination is >0.999999 and 4, respectively in terms of F1 score (harmonic mean of recall and precision).   

## <a name="ginkgo"></a>Ginkgo ##

### Installing the software. ###

1. Download command line version of Ginkgo from https://github.com/robertaboukhalil/ginkgo. Install according to the instruction. 

    In particular, to use our modified version of ginkgo (more details in #3), install the R packages, ctc, DNAcopy and inline by the following R command in R (we use R/4.0.0):

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("ctc")
    BiocManager::install("DNAcopy")
    install.packages("inline")
    install.packages("gplots")
    install.packages("scales")
    install.packages("plyr")
    install.packages("ggplot2")
    install.packages("gridExtra")
    install.packages("fastcluster")
    install.packages("heatmap3")

    In the directory of ginkgo, type "make". This will build all the executable files from the C++ files. 

2. In the directory of Ginkgo, in the folder of /genomes, make a subdirectory with name hg19. Inside /genomes/hg19 run the following command:

    ```wget http://qb.cshl.edu/ginkgo/uploads/hg19.original.tar.gz -O - | tar -xz```

3. To tune alpha, which controls the significance level to accept a change point, use our modified version of ginkgo.sh and process.R, which can be found in run_ginkgo/ folder in this repository. Put ginkgo.sh in cli/, and process.R in scripts/, respectively, in the ginkgo command line folder you downloaded. If you don't have to tune alpha, use the original ginkgo.sh and process.R, but eliminate the --alpha option.  

### Run Ginkgo. Inputs are bam files. Output is one big file holding all copy number changes for all cells. ###

1. Prepare the bed files. 

    For each bam file you have, run

    ```bedtools bamToBed -i $dir/$id.bam > $dir/$id.bed```

    You can download and install bedtools in https://bedtools.readthedocs.io/en/latest/. 

    $dir is the directory holding all the bam files.

2. Gather the bed files and output them into cell.list file.

    ```ls $dir/*.bed > $dir/cell.list``` 

3. Run Ginkgo. 

    ```cli/ginkgo.sh --input $dir/ --genome hg19 --binning variable_175000_48_bwa --maskbadbins --cells $dir/cell.list --alpha $alpha```

    We use variable_175000_48_bwa after testing the performance on different options. We found this gives us good accuracy. We use bwa because our alignment was done by bwa. Users can choose bowtie if they align with bowtie.  

    Output file can be found in $dir/SegCopy. Each row is a genomic region. The columns are chromosome, start, end, followed by the absolute copy number inference for each cell for this region. 

## <a name="copynumber"></a>CopyNumber ##

### Installing the software. ###

1. Install CopyNumber package on your R, as described in https://bioconductor.org/packages/release/bioc/html/copynumber.html.  

2. Download CopyNumber.R from run_CopyNumber/ in this depository. This is the script we wrote that calls winsorize and multipcf.segments, two key steps in CopyNumber package.  

3. CopyNumber does not take bam file direclty as input. Instead, it takes the read counts. We wrote a script that converts the csv file generated by readCounter and correct_read_count.py in HMMcopy to a file whose format is compatible to CopyNumber's input. User can find this script in run_CopyNumber/collect_rc_4copynumber.py in this depository.  Put this file together with CopyNumber.R.

### Run CopyNumber. ###

1. Prepare the input file.  

    For each bam file, run readCounter and correct_read_count.py steps in HMMcopy (see the instruction above in HMMcopy section). Put all "*.wig.csv" files together in one folder $dir. Note: the suffix of the files generated has to be wig.csv as our script uses this to recognize the files to be converted to the input file for CopyNumber. 

    Run

    ```python collect_rc_4copynumber.py $dir``` 

    This will create a file called copynumber.input.csv in $dir. It will be the input file for CopyNumber.R.

2. Run CopyNumber.R.

    ```Rscript $this_dir/CopyNumber.R $dir/copynumber.input.csv $output_f $gamma```
    
    $output_f is the output file, which will be put in $dir/. 

    $gamma is the weight of the penalty on changing a state. For more details about how to tune $gamma, refer to our ROC analysis in Fig. 1. 

    If using default gamma, which is 40, run

    ```Rscript $this_dir/CopyNumber.R $dir/copynumber.input.csv $output_f```

## <a name="aneufinder"></a>AneuFinder ##

### Installing the software. ###

    In R, type 

        ```
        if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

        BiocManager::install("AneuFinder")
        BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
        ```

### Running AneuFinder. ###

    In R, type

        ```
        library(AneuFinder)

        library(BSgenome.Hsapiens.UCSC.hg19)

        var.width.ref <- system.file("extdata", "hg19_diploid.bam.bed.gz", package="AneuFinderData")
        blacklist <- system.file("extdata", "blacklist-hg19.bed.gz", package="AneuFinderData")
        datafolder<-"ploidy${ploidy}_1"
        outputfolder<-"ploidy${ploidy}_1_out"
        Aneufinder(inputfolder = datafolder, outputfolder = outputfolder, assembly = 'hg19', numCPU = 1, binsizes = c(2e5), chromosomes = c(1:22,'X','Y'), blacklist = blacklist, correction.method = 'GC', GC.BSgenome = BSgenome.Hsapiens.UCSC.hg19, refine.breakpoints=FALSE, method = 'edivisive')
        ```

    Here ${ploidy} refers to the dataset being processed. In ploidy${ploidy}_1 are the bam files of all the cells to be processed. Folder ploidy${ploidy}_1_out will contain the results from AneuFinder. In particular, the result file is in ploidy${ploidy}_1_out/BROWSERFILES/method-edivisive/binsize_2e+05_stepsize_2e+05_CNA.bed.gz. This file will be further processed in section "Generating violin plots". 


### Run AneuFinder. ###
    
# <a name="misc"></a>Miscellaneous. #
## <a name="reference_file"></a>Referece file ##

All reference files used in this paper is hg19. The scripts to download it and to process the fa files for each chromosome are below.

```wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz```

```tar -vxf chromFa.tar.gz```

```for i in `seq 1 22` X Y; do cat chr$i.fa >> hg19.fa; done```

## <a name="mapping"></a>Mapping the reads to the reference ##

For both simulated and real data, we use bwa to align the reads to the reference. We eliminated reads with mapping quality score < 40 in creating the bam file. The following are the commands we used to generate the bam files.

For simulated data, we have an extra step that merges the fastq files corresponding to the paternal and maternal alleles for each end of the paired end reads. Command as follows.

```cat $dir/leaf${n}_allele0_1.fq $dir/leaf${n}_allele1_1.fq > $dir/leaf${n}_1.fq```

```cat $dir/leaf${n}_allele0_2.fq $dir/leaf${n}_allele1_2.fq > $dir/leaf${n}_2.fq```

${n} is the index of the leaf. Here we use "leaf" as the prefix of the file names generated by the simulator. $dir is the directory that contains the simulated fastq files.

The following step is the same for both simulated and real data.

gen_bam/make_bam_from_fq.sh $dir/leaf${n} $hg19 $processor $mapping_qual_t 

This script requires installing bwa and samtools.

$hg19 is the reference fasta file in the absolute path. $processor is the number of processor to run bwa. $mapping_qual_t is the threshold of mapping quality. We set it to be 40 for all experiments in this paper. 

The outputs of this step are the sorted bam (duplication removal step was also performed in this script) and the bai file, with the names $dir/leaf${n}.sorted.bam[.bai].  

## <a name="ground_truth"></a>Making ground truth from the simulator for comparison. ##

1. Read the from_first_step.tree.npy file generated in the first step of the simulator and convert it to a csv file. 

    ```python read_tree.py -s -f from_first_step.tree.npy > gt.all.csv```

    This step generates a file in bed format that contains all ground truth CNAs for each cell. The fourth column (1-based) is the cell ID.

2. Generate the ground truth for each individual cell from gt.all.csv.

    ```mkdir gt_sep```

    ```python comparison/sep_groundtruth.py gt.all.csv gt_sep/gt``` 

    This script will generate ground truth file for each cell in gt_sep folder, with the prefix gt. 

    For Ginkgo and HMMcopy, we compare the result with the ground truth for each cell separately. For CopyNumber, the ground truth is the combination of all cells that are involved in the study. We eliminate the CNAs in ground truth whose supporting cells are less than five, for CopyNumber alone. The commands are as follows.

    ```perl extract_gt.pl gt.all.csv selected.leaves > gt.selected.csv```

    ```perl -ane 'print join("\t", @F[0 .. 2]) . "\n"' gt.selected.csv | sort | uniq -c | perl -ane 'print join("\t", @F[1 .. 3]) . "\n" if($F[0] > 5)' > gt.forCopyNumber.csv``` 

    This step generates gt.forCopyNumber.csv as the ground truth to be compared to CopyNumber's results.

3. If you want to generate segcopy formatted ground truth file for comparison, use the following command. 

    ```python bin_groundtruth.py -a segcopy_f -b gt.all.csv --leafonly(optional) > gt.all.segcopyformatted```

    segcopy_f is a file you generated from Ginkgo (see Ginkgo under commands to run HMMcopy, Ginkgo and CopyNumber (#commands_3methods) for details. gt.all.csv is the file you generated in step #1 in this section. If you want to include all ancestral nodes along with the leaves, do not put the option --leafonly in the command line. Otherwise use --leafonly in your command. 

    The format of the output is "chr, start, end, leaf#1, leaf#2, ...". From the fourth column (1-based), the entries are integer copy numbers. 

## <a name="newick"></a>Generate a newick-formatted tree from .npy file in simulation. ##

    The simulator automatically stores the tree structure in from_first_step.tree.npy. To generate the newick string of the tree, use the following command.

    ```python gen_newick.py from_first_step.tree.npy > newicktree.txt```

# <a name="plots"></a>Generating plots (violin plots, plots with Lorenz and Beta distribution, ROCs, venn diagram and flip count histograms). #

## <a name="violin_plots"></a>Generating violin plots ##

1. Prepare the files.

    * In ~/.bashrc, add one line

        ```export PERL5LIB=$this_dir/perl_modules```

        ```source ~/.bashrc```

        in which $this_dir is the directory of this repository. This step will take the perl module file automatically.  

    * Prepare HMMcopy files. In a directory $dir/HMMcopy,

        ```for i in `seq 1 $n`; do mkdir $dir/HMMcopy/p$i; done```

        $dir="data/$exp/rep_${rep}"

        in which $exp is "fluc" for the fluctuation experiment and "ploidy" for the ploidy experiment. ${rep} is the index of the repitition. 
    
        In each $dir/HMMcopy/p[1-$n] folder, put the corresponding HMMcopy's result whose prefix is leaf, and the name of such a file is leaf$ID.csv for the cell whose ID is $ID. Here $n represents the number of ploidies/fluctuations. In our paper, $n = 5 for ploidy; $n = 4 for fluctuations. 

        Next, convert from csv to bed format. 
    
        ```for i in `seq 1 $n`; do python scripts_plots/general_scripts/convertHMMcopySeg2bed_pop.py $dir/HMMcopy/p$i; done``` 

    * Prepare Ginkgo files. In a directory $dir/Ginkgo, 

        ```for i in `seq 1 $n`; do mkdir $dir/Ginkgo/p$i; done```

        In each $dir/Ginkgo/p[1-$n] folder, put the corresponding SegCopy file from Ginkgo's output. 

        Convert SegCopy to bed file and separate the results for each cell.

        ```for i in `seq 1 $n`; do python scripts_plots/general_scripts/convertGinkgoSegCopy2bed_pop.py $dir/Ginkgo/p$i/SegCopy; done``` 

        $n has the same meaning as described in HMMcopy's section.

    * Prepare CopyNumber files.

        In a folder $dir/CopyNumber, 

        ```for i in `seq 1 $n`; do mkdir $dir/CopyNumber/p$i; done```

        In each $dir/CopyNumber/p$i folder, put the corresponding CopyNumber output file.  

    * Prepare AneuFinder files.

        In a folder $dir/AneuFinder,

        ```for i in `seq 1 $n`; do mkdir $dir/AneuFinder/p$i; done```

        In each $dir/Ginkgo/p[1-$n] folder, put the corresponding binsize_2e+05_stepsize_2e+05_CNV.bed file from AneuFinder's output (which can be found in $output_folder/BROWSERFILES/method-edivisive/ in .gz format).

        Convert the bed file and separate the results for each cell.

        ```for i in `seq 1 $n`; do python scripts_plots/general_scripts/convertAneuFinder2bed_pop.py $dir/AneuFinder/p$i/binsize_2e+05_stepsize_2e+05_CNV.bed leaf; done```

        $n has the same meaning as described in HMMcopy's section.

    * Prepare ground truth file.

        In a folder $dir/gt, 

        ```for i in `seq 1 $n`; do mkdir $dir/gt/p$i; done```

        In each $dir/gt/p$i folder, put the corresponding gt.all.csv and gt_sep folder. For more details about the ground truth file, please read miscellaneous section. 

2. To generate violin plots, organize the data in such as a way that $dir is in parallel to the folder holding the scripts generating the plots. The scripts include fluc_qual.py, fluc_quan.py, ploidy_qual.py and ploidy_quan.py. They can be found in script_plots/ folder in this repository. For an example of how to organize the files, see the structure in script_plots/.  

    * To generate violin plots for qualitative comparison for ploidy experiment, run

        ```python ploidy_qual.py ploidy $rep 5 5```

        The first 5 represents the number of ploidies being investigated. The second 5 is the threshold of the cells supporting a ground truth CNA, the result of which will be used to evaluate CopyNumber's performance.

    * To generate violin plots for quantitative comparison for ploidy experiment, run

        ```python ploidy_quan.py ploidy $rep 5```

        In quantitative analysis we don't evaluate CopyNumber as it does not output absolute copy number.

    * To generate violin plots for qualitative comparison for fluc experiment, run

        ```python fluc_qual.py fluc $rep 4 5```

        The first number 4 represents the number of levels of fluctuation being investigated. The second number 5 is the threshold of the cells supporting a ground truth CNA, the result of which will be used to evaluate CopyNumber's performance.

    * To generate violin plots for quantitative comparison for fluc experiment, run

        ```python fluc_quan.py fluc $rep 4```

        In quantitative analysis we don't evaluate CopyNumber as it does not output absolute copy number.

    * To generate violin plots for AneuFinder, run

        ```python scripts_plots/violin/scripts/ploidy_qua[l|n]_AneuFinder.py ploidy 1 1 5```

## <a name="scatter_plots"></a>Generating scatter plots of ploidies with histogram. ##

1. Prepare a meta file, in which each line represents a ploidy, and there are three columns. The first column is the file name (with the path) of the ground truth ploidy file, the second the file name (with the path) of the inferred ploidy file, and the third a color denoting this ploidy. The aforementioned two files should have the same format. Both will have two columns, the first the cell name, the second the ploidy corresponding to the cell. The number of rows is the number of cells in the study. The two files should have the same set of cells. The following is an example of the meta file.

    ```
    $path/gt_ploidy_p1.csv $path/ginkgo_ploidy_p1.csv    purple 
    $path/gt_ploidy_p2.csv $path/ginkgo_ploidy_p2.csv    green 
    $path/gt_ploidy_p3.csv $path/ginkgo_ploidy_p3.csv    orange 
    $path/gt_ploidy_p4.csv $path/ginkgo_ploidy_p4.csv    red 
    $path/gt_ploidy_p5.csv $path/ginkgo_ploidy_p5.csv    dodgerblue
    ```

    The following is a few lines to show the format of actual or inferred ploidy file. 

        leaf11   1.56
        leaf14   1.65

    These files can be obtained by running scripts in 

        ```python utilities/calc_ploidy_[gt.py|segcopy.py|hmmcopy.sh]```

2. Plot the scatterplot and histograms. 

    ```python scripts_plots/ploidy_comparison/scatterplot_hmmcopy.py $meta_file $output_fig.png```

    $output_fig.png is the output figure file. (a) is the summary of the five ploidies. (b)-(e) are the scatter plots for each ploidy, each having a histogram showing the percentage of each predicted ploidy value on the right. 

## <a name="read_count_copy_number"></a>Generating plots on the whole genome with read count and absolute copy number from ground truth or a method. ##

1. Prepare the file for plot.

    * For plotting ground truth together with HMMcopy's inference of absolute copy number, run

        ```perl scripts_plots/wgplots/bam_gt_est_forplot.pl $fai $w $gt_f $bam $hmmcopy_segs > $file_for_plot```
        
        $fai is the reference fai file in absolute path. 
        
        $w is the window size for each dot representing the read count (the primary Y axis). We use 200,000 in our paper.  
        
        $gt_f is the ground truth file generated by the simulator (see section "Making ground truth from the simulator for comparison."). Specifically, the file to be used contains the ground truth for this cell only, i.e., the one produced by the script comparison/sep_groundtruth.py.
        
        $bam is the bam file of the cell to be plotted.
        
        $hmmcopy_segs is HMMcopy's segs.csv file.  

        $file_for_plot is the file that combines the read count, ground truth and HMMcopy's inference for each bin. 

    * For plotting the inferred absolute copy number from two candidate ploidies of HMMcopy, run 

        ```perl scripts_plots/wgplots/bam_gt_twohmmcopy_forplot.pl $fai $w $bam $hmmcopy_segs_1 $hmmcopy_segs_2 > $file_for_plot```
        
        $hmmcopy_segs_1 and $hmmcopy_segs_2 represent the segs.csv from two candidate ploidies of HMMcopy.
        
        $file_for_plot is the file that combines the read count, HMMcopy's inference for candidate ploidy 1 and 2 for each bin. 

2. Plot the read count on the whole genome together with the absolute copy numbers. 

    ```Rscript scripts_plots/wgplots/plot_rc_cn.r $file_for_plot $pdf_output $legend1 $legend2```

    $file_for_plot is the output from step 1.

    $pdf_output is the name of the pdf file containing the plot.

    $legend1 and $legend2 represent the label for the two columns with the absolute copy number in $file_for_plot. For the first case discussed in step 1, they are "Ground Truth" and "HMMcopy", respectively. For the second case discussed in step 1, they are "HMMcopy (candidate ploidy = x)" and "HMMcopy (candidate ploidy = y)". Please fill in x and y according to the data you use. Note that the scale of secondary Y axis needs to be tuned to align with the primary Y axis scale.  

## <a name="lorenz_beta"></a>Generating the plot including both Lorenz curve and the corresponding Beta distribution ##

1. Generate separate plots in one PDF for different levels of fluctuation.

    ```Rscript scripts_plots/lorenz_beta/plot_lorenz_beta_separate.r $pdf_output```

    $pdf_output is the output pdf file.

2. Generate a summary plot in one PDF for different levels of fluctuation.

    ```Rscript scripts_plots/lorenz_beta/plot_lorenz_beta_summary.r $pdf_output```

    $pdf_output is the output pdf file.

## <a name="venn_diagram"></a>Generating venn diagram ##
1. Generate the sets of breakpoints from the three tools.

    ```python scripts_plots/venn/Venn_data.py $thr $hmmcopy_result $ginkgo_result $copynumber_result > $output.txt```

    $thr is the threshold which determines whether two breakpoints from two tools are the same or not. We have used $thr=400000 in our analysis in the paper. 

    $hmmcopy_result is the output of HMMcopy.

    $ginkgo_result is the output of Ginkgo.

    $copynumber_result is the output of CopyNumber.

    $output.txt is the input file for the next step.

2. Plot the Venn diagram in Python.

    ```python scripts_plots/venn/Venn_diagram_plot.py $output.txt $fig_name```

    $fig_name is the name of the output figure.

# <a name="PAUP"></a>Using PAUP to generate a maximum parsimony tree and count flips. #

## Installing and running PAUP. ##

1. Download PAUP from http://phylosolutions.com/paup-test/. In the website, under Command-line binaries, you can find the command-line versions of PAUP according to your OS. In case you need the GUI version of PAUP, download the file named PAUP_dev_icc.zip.

2. Unzip the package, open terminal and go the directory of the unzipped packaged. To make the command-line executable, run

    ```chmod a+x paup4a166_osx```

## Inferring the Maximum Parsimony tree using PAUP 

1. Run the script named Nexus_gen.py to generate the input file of PAUP.

    ```python Nexus.py $infile $nexus_file.nex```

    $infile could be the output of HMMcopy, Ginkgo, or CopyNumber

    $nexus_file.nex is the output of Nexus.py which will be given as input to PAUP

2. Open PAUP (in the following case, the binary version is that of Macintosh OS X)

    ```./paup4a166_osx```

3. In the shell, type 

    ```Execute $nexus_file.nex```

    $outfile is the path to the output of Nexus.py in the last section. The output of this step including the information of the tree will be saved in a file named logfile.log

4. Parse the output of PAUP (logfile.log).

    ```perl PAUP/process_PAUP.pl logfile.log $parsed```

    $parsed is the output of this script.

## Counting the number of copy number changes (or flips) across the tree ##
    
    Run the script named flip_counter.py to generate the histogram of the flip-counts across the tree

    ```python flip_counter.py $parsed $hist_name```

    $hist_name is the name of the histogram plot.
