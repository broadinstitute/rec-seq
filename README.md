## Rec-seq
Rec-seq is a method for determining the DNA specificity and potential off-target substrates of site-specific 
recombinases via high-throughput sequencing of recombined substrates from a pool of partially randomized DNA 
sequences. The Rec-seq data analysis script, rec-seq.py, is used to quantify 
the enzyme’s substrate specificity. Post-recombination sequencing reads that contain the matched target core 
sequence are aligned to the native target sequence, with no gaps allowed. After alignment, reads with excessive 
numbers of mismatches are considered to be the result of sequencing errors and are filtered out of subsequent 
analysis (determined by max_mismatch_count parameter). For the remaining sequences, at each position in the 
recombinase target, the abundance of the canonical base (A<sub>i</sub>) and the sum of the non-canonical bases 
(B<sub>i</sub>) are calculated. The same analysis is performed for the sequencing reads of the input library, 
except the abundance of the canonical base and the sum of the non-canonical bases are expressed as fractions 
α<sub>i</sub> and β<sub>i</sub>. The enrichment score for each position is then calculated as the ratio 
r<sub>i</sub> = (A<sub>i</sub>/B<sub>i</sub>)/(α<sub>i</sub>/β<sub>i</sub>). Analysis is performed separately 
for the left and right half-sites, using as input the sequencing reads from experiments with either left- or 
right-randomized half-sites.

### Input
1. **index file** - tab-delimited text file defining analysis of individual experiments. Each row in the index file specifies an enzyme variant used in the experiment, a date of the experiment, and fastq files containing pre- and post-selection sequencing reads for both left and right half sites. Index file contains the following columns:
   * `Enzyme_variant` - enzyme variant used in the experiment
   * `Date` - date of the experiment
   * `Left_library_file` - fastq file containing post-selection sequencing reads for the left half site
   * `Right_library_file` - fastq file containing post-selection sequencing reads for the right half site
   * `Left_control_file` - fastq file containing pre-selection sequencing reads for the left half site
   * `Right_control_file` - fastq file containing pre-selection sequencing reads for the left half site

2. **substrate file** - tab-delimited text file defining canonical substrate sequences and their formats. Each file listed in the index file must have an entry in the substrates file. Library files and their corresponding control files must have the same site layout. Substrate file contains the following columns:
   * `Library_file` - fastq file containing sequencing reads
   * `Substrate` – nucleotide sequence of the canonical substrate
   * `Site_layout` - layout of the substrate in the format <left half-site length>;<core length>;<right half-site length>

3. **fastq files** (listed in index/substrate files) containing sequencing reads.

4. **max_mismatch_count** - parameter that controls which sequencing reads that are included in the analysis, only reads that have at most max_mismatch_count mismatches when compared to left and right half sites of the canonical substrate are considered.

### Output
Analysis output, printed to stdout, is a tab-delimited text with the following columns:
* `Enzyme_variant` - enzyme variant used in analyzed experiment
* `Date` - date of analyzed experiment
* `Left_library_file` - fastq file containing post-selection sequencing reads for the left half site 
* `Right_library_file` - fastq file containing post-selection sequencing reads for the right half site
* `Total_read_count` - total number of reads present in post-selection library file
* `Control_read_count` - total number of reads present in pre-selection control file
* `Core_count` - number of reads containing core sequence in post-selection library file
* `Control_core_count` - number of reads containing core sequence in pre-selection control file
* `Position` - position for which enrichment was computed, negative for the left half-site and positive for the right half-site.
* `Nucleotide` - nucleotide for which enrichment was computed
* `Match` - true when the nucleotide for which enrichment was computed matches a canonical substrate base at a given position 
* `Library_count` - number of reads (containing core sequence and having at most max_mismatch_count mismatches) matching a given nucleotide at a given position in post-selection library file
* `Control_count` - number of reads (containing core sequence and having at most max_mismatch_count mismatches) matching a given nucleotide at a given position in pre-selection control file
* `Enrichment` - computed value of enrichment of a given nucleotide at a given position

### Usage
    rec-seq.py [-h] [-mmc MAX_MISMATCH_COUNT] index_file substrate_file
    
    positional arguments:
      index_file            input index file
      substrate_file        input substrates file
    
    optional arguments:
      -h, --help            show help message and exit
      -mmc MAX_MISMATCH_COUNT, --max_mismatch_count MAX_MISMATCH_COUNT
                        default max mismatch count value 5

### Requirements
Python 3.2 or later
