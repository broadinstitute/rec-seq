#-------------------------------------------------------------------------------
# Name:        rec-seq.py
# Purpose:     Computation of Rec-seq enrichment scores
#
# Authors:      Vlado Dancik and Jeffrey L. Bessen
#
# Copyright:   (c) 2019 Broad Institute
# License:     https://github.com/broadinstitute/rec-seq/blob/master/LICENSE
#-------------------------------------------------------------------------------

import sys
import argparse
import os.path
import pandas as pd

TAB = '\t'

HEADER = str('Enzyme_variant' + TAB
            +'Date' + TAB
            +'Left_library_file' + TAB
            +'Right_library_file' + TAB
            +'Total_read_count'+TAB
            +'Control_read_count' + TAB
            +'Core_count' + TAB
            +'Control_core_count' + TAB
            +'Position' + TAB
            +'Nucleotide' + TAB
            +'Match' + TAB
            +'Library_count' + TAB
            +'Control_count' + TAB
            +'Enrichment')


def log(string):
    """ Print a string to error stream. """
    print(string, file=sys.stderr)


def parse_substrate_format(substrate_format, sep=';'):
    """
        Parse site-layout string
    """
    parsed_format = substrate_format.split(sep)
    if len(parsed_format) == 3:
        return [int(x) for x in parsed_format]
    else:
        return None


def substrate_map(substrate_file):
    """
        Load a substrate file, parse substrate formats, identify core sequence,
        and return a map indexed by a substrate file name
    """

    map = {}
    substrates =  pd.read_csv(substrate_file, sep='\t')
    if 'Library_file' not in substrates.columns:
        log('ERROR: required column "Library_file" not found in input substrates file')
    elif 'Substrate' not in substrates.columns:
        log('ERROR: required column "Substrate" not found in input substrates file')
    elif 'Site_layout' not in substrates.columns:
        log('ERROR: required column "Site_layout" not found in input substrates file')
    else:
        for i in range(len(substrates)):
            library_file = substrates.Library_file[i]
            format = parse_substrate_format(substrates.Site_layout[i], sep=';')
            substrate = substrates.Substrate[i]
            core = substrate[format[0] : (format[0]+format[1])]
            map[library_file] = (substrate, core, format)
        return map
    return None


def next_read(f):
    """
       Get a next sequencing read from the input file f.
    """

    f.readline()
    line = f.readline()
    f.readline()
    f.readline()
    return line


def find_core(read, core, core_position_sum, core_position_count, start = -1):
    """
       Find the core sequence, trying "average" position first for efficiency.
    """

    if start < 0 and core_position_count > 0:
        core_position = round(core_position_sum/core_position_count)
        if len(read) > core_position+len(core):
            if read[core_position:core_position+len(core)]==core:
                return core_position
    return read.find(core, start+1)



def match(substrate, target):
    """
        Compare two sequences of the same length and return indexes of the
        positions where the sequences differ.
    """

    if len(substrate) != len(target):
        log("WARNING: len(substrate) != len(target) "+str(len(substrate)) +'!='+ str(len(target)))
        return None
    diff = []

    for i,(x,y) in enumerate(zip(substrate, target)):
        if x!=y:
            diff.append(i)
    return diff


def find_substrate(read, target, position, max_mismatch_count):
    """
       Return a substrate at a given position if it has fewer then max_mismatch_count
       mismatches from the target sequence.
    """

    if position > 0 and len(read) >= position+len(target):
        substrate = read[position : position+len(target)]
        if len(match(substrate, target)) <= max_mismatch_count:
            return substrate

    return None


def nucleotide_index(nucleotide):
    """
        Convert nucleotide character to a number.
    """

    if nucleotide == 'C':
         return 2
    if nucleotide == 'G':
         return 3
    if nucleotide == 'A':
         return 1
    if nucleotide == 'T':
         return 4
    return 0


def nucleotide(nucleotide_index):
    """
        Convert nucleotide index to a character.
    """

    nucleotides = ['?','A','C','G','T']
    if 1 <= nucleotide_index and nucleotide_index <= 4:
        return nucleotides[nucleotide_index]
    return '?'


def nucleotide_distribution(filename, substrate, core, offset, max_mismatch_count):
    """
        Compute nucleotide distribution
    """

    # we keep track of "average" core position to find it efficiently
    core_position_sum = 0
    core_position_count = 0

    # accumulate nucleotide distribution per position
    # start with one count for each nucleotide, aka "uniform prior"
    nucleotide_counts = [[1]*5 for i in range(len(substrate))]

    read_count = 0
    core_count = 0
    with open(filename) as f:
        read = next_read(f)
        read_count += 1

        while read:
            core_position = find_core(read, core, core_position_sum, core_position_count)
            if core_position > 0:
                core_count += 1

            while core_position > 0:
                #loop through multiple possible core positions
                target = find_substrate(read, substrate, core_position-offset, max_mismatch_count)
                if target == None:
                    core_position = find_core(read, core, core_position_sum, core_position_count, core_position)
                else:
                    for pos, nucleotide in enumerate(target):
                        ni = nucleotide_index(nucleotide)
                        nucleotide_counts[pos][ni] += 1
                    core_position_sum += core_position
                    core_position_count += 1
                    core_position = -1

            read = next_read(f)
            read_count += 1

    log(filename+': read_count = '+str(read_count)+', core_count = '+str(core_count)
        +', match_count = '+str(core_position_count)+', ratio = '+str(core_count/read_count))
    return (read_count, core_count, nucleotide_counts)


def get_control(controls, control_file, substrates, max_mismatch_count, root):
    """
        Get control distribution from controls map, computing it from
        a control file if necessary
    """

    control_distribution = controls.get(control_file)
    if control_distribution == None:
        (substrate,core,format) = substrates[control_file]
        control_distribution = nucleotide_distribution(root+control_file, substrate, core, format[0], max_mismatch_count)
        controls[control_file] = control_distribution
    return control_distribution


def enrichment_score(filename, control, substrate, format, max_mismatch_count):
    """
        Compute enrichment scores
    """

    core = substrate[format[0]:(format[0]+format[1])]
    (read_count, core_count, nucleotide_counts) = nucleotide_distribution(filename, substrate, core, format[0], max_mismatch_count)
    (control_read_count, control_core_count, control_nucleotide_counts) = control
    enrichment = [[0]*5 for i in range(len(substrate))]

    for i in range(len(substrate)):
        total = sum(nucleotide_counts[i])
        control_total = sum(control_nucleotide_counts[i])

        for j in range(1,5):
            frequency = nucleotide_counts[i][j]/total
            odds_ratio = frequency / (1-frequency)
            control_frequency = control_nucleotide_counts[i][j]/control_total
            control_odds_ratio = control_frequency / (1-control_frequency)
            enrichment[i][j] = odds_ratio/control_odds_ratio

    return (read_count, core_count, nucleotide_counts, enrichment)


def print_enrichments(print_range, row, enrichments, control_distribution, substrate, format):
    """
        Print computed enrichments to stdout
    """

    (read_count, core_count, nucleotide_counts, enrichment) = enrichments
    (control_read_count, control_core_count, control_nucleotide_counts) = control_distribution

    for k in print_range:
        for j in range(1,5):
            print(row.Enzyme_variant + TAB
                + row.Date + TAB
                + row.Left_library_file + TAB
                + row.Right_library_file + TAB
                + str(read_count) + TAB
                + str(control_read_count) + TAB
                + str(core_count) + TAB
                + str(control_core_count) + TAB
                + str(k - int(format[0]+format[1]/2)) + TAB
                + nucleotide(j) + TAB
                + str(substrate[k]==nucleotide(j)) + TAB
                + str(nucleotide_counts[k][j]) + TAB
                + str(control_nucleotide_counts[k][j]) + TAB
                + str(enrichment[k][j]))


def process(row, substrates, controls, max_mismatch_count, root):
    """
        Compute enrichment scores for an experiment defined by a row in the
        input index file
    """

    log('Processing '+row.Date+' '+row.Enzyme_variant)
    # check that all files have substrate info available
    if row.Left_control_file not in substrates:
        log('WARNING: substrate info not available for '+row.Left_control_file)
    elif row.Left_library_file not in substrates:
        log('WARNING: substrate info not available for '+row.Left_library_file)
    elif row.Right_control_file not in substrates:
        log('WARNING: substrate info not available for '+row.Right_control_file)
    elif row.Right_library_file not in substrates:
        log('WARNING: substrate info not available for '+row.Right_library_file)
    else:

        # enrichments for left half-site
        control_distribution = get_control(controls, row.Left_control_file, substrates, max_mismatch_count, root)
        (substrate,core,format) = substrates[row.Left_library_file]
        enrichments = enrichment_score(root+row.Left_library_file, control_distribution, substrate, format, max_mismatch_count)
        print_enrichments(range(0,format[0]), row, enrichments, control_distribution, substrate, format)

        # enrichments for right half-site
        control_distribution = get_control(controls, row.Right_control_file, substrates, max_mismatch_count, root)
        (substrate,core,format) = substrates[row.Right_library_file]
        enrichments = enrichment_score(root+row.Right_library_file, control_distribution, substrate, format, max_mismatch_count)
        print_enrichments(range(format[0]+format[1],len(substrate)), row, enrichments, control_distribution, substrate, format)


def read_index_file(index_file):
    """
        Read index file and check that all required columns are present
    """
    input = pd.read_csv(index_file, sep='\t')
    if 'Enzyme_variant' not in input.columns:
        log('ERROR: required column "Enzyme_variant" not found in input index file')
    elif 'Date' not in input.columns:
        log('ERROR: required column "Date" not found in input index file')
    elif 'Left_library_file' not in input.columns:
        log('ERROR: required column "Left_library_file" not found in input index file')
    elif 'Right_library_file' not in input.columns:
        log('ERROR: required column "Right_library_file" not found in input index file')
    elif 'Left_control_file' not in input.columns:
        log('ERROR: required column "Left_control_file" not found in input index file')
    elif 'Right_control_file' not in input.columns:
        log('ERROR: required column "Right_control_file" not found in input index file')
    else:
        return input
    return []


def main():

    (index_file, substrate_file, max_mismatch_count) = parse_command_line()
    input = read_index_file(index_file)
    substrates = substrate_map(substrate_file)
    if len(input) == 0 or substrates == None:
        return

    root = os.path.dirname(index_file)
    if len(root) > 0:
        root = root + '/'
    print(HEADER)
    controls = {}

    for i in range(len(input)):
        try:
            process(input.iloc[i], substrates, controls, max_mismatch_count, root)
        except KeyboardInterrupt:
            raise
        except:
            type, value, traceback = sys.exc_info()
            log('WARNING: '+str(type)+' '+str(value))


def parse_command_line():
    """
        Parse command-line arguments
    """
    parser = argparse.ArgumentParser(description='compute rec-seq enrichment scores')
    parser.add_argument('index_file', help='input index file')
    parser.add_argument('substrate_file', help='input substrates file')
    parser.add_argument('-mmc', '--max_mismatch_count', type=int, default=5,
                        help='default max mismatch count value 5')
    args = parser.parse_args()
    log(args)
    return (args.index_file, args.substrate_file, args.max_mismatch_count)


if __name__ == '__main__':
    assert sys.version_info >= (3, 2)
    main()
