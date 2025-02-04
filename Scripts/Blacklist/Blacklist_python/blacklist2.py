"""
blacklist2.py
improved blacklist generation

Copyright (c) 2024 Brydon P. G. Wall

Original Blacklist by Alan Boyle

This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions.
"""
import math
import os
import argparse
import pysam
import statistics
import re
import numpy as np

class SequenceData():
    def __init__(self, bam_file: str = None, bam_index_file: str = None) -> None:
        self.bins_input: np.ndarray[np.uint32]
        self.bins_multimapping: np.ndarray[np.uint32]
        self.bins_temp: np.ndarray[np.uint32]
        self.total_reads: int = 0
        self.bam_file: str = bam_file
        self.bam_index_file: str = bam_index_file

    def get_input_bins_per_region(
            self, 
            mappability: np.ndarray[np.uint8], 
            bin_size: int, 
            bin_overlap: int, 
            referece_name: str
            ) -> None:

        temp_count_list = np.zeros(len(mappability), dtype=np.uint32) # [0 for i in range(len(mappability))]
        multi_count_list = np.zeros(len(mappability), dtype=np.uint32) # [0 for i in range(len(mappability))]
        z_multi_center = 0

        bam = pysam.AlignmentFile(filename=self.bam_file, index_filename=self.bam_index_file)
        
        # This restricts output to a specific chromosome
        if bam.get_tid(referece_name) != -1:

            # Keep vector of locations and counts of reads
            # filter on mappability based on each read length
            for al in bam.fetch(region=referece_name):
                
                # Filter mappability
                if mappability[al.reference_start] > 0 and mappability[al.reference_start] <= al.query_alignment_length:
                    
                    # Count reads at each position that are uniquely mappable
                    temp_count_list[al.reference_start] = temp_count_list[al.reference_start] + 1

                # A Multimapping read
                else:
                    
                    # Count reads at each position that are not uniquely mappable
                    multi_count_list[al.reference_start] = multi_count_list[al.reference_start] + 1
                    z_multi_center += 1
                
                self.total_reads += 1

        bam.close()
        
        array_size = len(temp_count_list) - bin_size
        self.bins_input = np.zeros(array_size, dtype=np.uint32)
        self.bins_multimapping = np.zeros(array_size, dtype=np.uint32)

        k = 0
        for i in range(0, array_size, bin_overlap):
            
            read_center = np.sum(temp_count_list[i:i+bin_size])
            multi_center = np.sum(multi_count_list[i:i+bin_size])

            self.bins_input[k] = read_center
            self.bins_multimapping[k] = multi_center
         
            k += 1

def get_mappability_bins(
        mappability: np.ndarray[np.uint8],
        bin_size: int,
        bin_overlap: int,
        unique_length: int = 36
        ) -> np.ndarray[np.uint32]:
    """
    Takes as input the mappability at each base and returns binned counts of mappability
    """
    num_bins = (len(mappability) - bin_size) // bin_overlap + 1
    mappability_bins = np.zeros(num_bins, dtype=np.uint32)

    for i in range(0, (len(mappability) - bin_size) + 1, bin_overlap):
        bin_values = mappability[i:i+bin_size]
        unique_center = np.count_nonzero((bin_values > 0) & (bin_values <= unique_length))
        mappability_bins[i // bin_overlap] = unique_center

    return mappability_bins

def rankify(bin_list):
    sorted_indices = np.argsort(bin_list)
    ranks = np.empty_like(bin_list)
    
    rank = 1
    i = 0

    while i < len(bin_list):
        j = i

        # Get elements of the same rank
        while j < len(bin_list) - 1 and bin_list[sorted_indices[j]] == bin_list[sorted_indices[j + 1]]:
            j += 1

        m = j - i + 1
        average_rank = rank + (m - 1) * 0.5

        for j in range(i, i + m):
            index = sorted_indices[j]
            ranks[index] = average_rank

        # Increment rank and index
        rank += m
        i += m

    return ranks

def quantile_normalize(data) -> None:
    cell_count: int = len(data)
    bin_count: int = len(data[0])

    # First calculate rank means
    ranked_mean = np.zeros((bin_count))
    for cell_ID in range(cell_count):
        x = np.zeros((bin_count))
        for bin_ID in range(bin_count):
            x[bin_ID] = data[cell_ID][bin_ID]
        np.ndarray.sort(x)
        for bin_ID in range(bin_count):
            ranked_mean[bin_ID] += x[bin_ID]
    for bin_ID in range(bin_count):
        ranked_mean[bin_ID] /= cell_count

    # calculate half value for ties
    ranked_mean_tie = np.zeros((bin_count - 1))
    for bin_ID in range(bin_count - 1):
        ranked_mean_tie[bin_ID] = (ranked_mean[bin_ID] + ranked_mean[bin_ID + 1]) / 2

    # Iterate through each cell line
    for c_count in range(cell_count):
        bins = np.zeros((bin_count))
        for b_count in range(bin_count):
            bins[b_count] = data[c_count][b_count]

        bins = rankify(bins)

        bins_q_normalized = np.zeros((bin_count))
        for b_count in range(bin_count):
            if bins[b_count] % 1 != 0:
                bins_q_normalized[b_count] = ranked_mean_tie[int(math.floor(bins[b_count])) - 1]
            else:
                bins_q_normalized[b_count] = ranked_mean[int(bins[b_count] - 1)]
            data[c_count][b_count] = bins_q_normalized[b_count]

    return data

def get_abnormal_regions(
        input_data: np.ndarray,
        bins_map: np.ndarray,
        region_type: str,
        normalize: bool,
        constant: float = 1000000.0
        ) -> np.ndarray:

    norm_temp = 0.0
    result = np.zeros(len(bins_map))
    data = np.zeros((len(input_data), len(bins_map)))
    
    # process by column
    for j in range(len(input_data)):

        if region_type == 'reads': # reads
            input_data[j].bins_temp = input_data[j].bins_input
        elif region_type == 'multi': # multimapping
            input_data[j].bins_temp = input_data[j].bins_multimapping

        # track all rows in this column
        for i in range(len(bins_map)):
            norm_temp = 0.0
            if normalize:
                if region_type == 'multi' and input_data[j].bins_input[i] != 0: # multimapping
                    norm_temp = input_data[j].bins_temp[i] / input_data[j].bins_input[i]
                elif bins_map[i] != 0:
                    norm_temp = input_data[j].bins_temp[i] / bins_map[i]

            elif input_data[j].total_reads != 0 and constant != 0:
                norm_temp = (input_data[j].bins_temp[i] / input_data[j].total_reads) * constant

            data[j][i] = norm_temp

    data = quantile_normalize(data)

    # Now we collapse rows
    means = np.zeros(len(input_data))
    
    # over each row
    for i in range(len(bins_map)):

        # over each column
        for j in range(len(input_data)):
            means[j] = data[j][i]

        # This is median signal
        result[i] = statistics.median(means)

    return result

def get_bed(
        read_norm_list: np.ndarray[np.float64],
        multi_list: np.ndarray[np.float64],
        bins_map: np.ndarray[np.uint32],
        bin_size: int,
        bin_overlap: int,
        region: str,
        bridge_constant:int = 200,
        output_original: bool = False,
        original_only: bool = False,
        weak_percentile: float = 0.99,
        strong_percentile: float = 0.999
        ) -> str:

    # Output the raw count information
    temp_first: int
    temp_last: int
    in_region: bool = False
    temp_hit: int = 0
    hit_code: int = 0
    hit_counter: int = 0
    miss: int = 0
    region_class: str = ""

    # calculate quantiles parameters
    weak_n = 10 ** len(str(weak_percentile).split('.')[-1])
    weak_i = int(weak_percentile * weak_n) - 1
    strong_n = 10 ** len(str(strong_percentile).split('.')[-1])
    strong_i = int(strong_percentile * strong_n) - 1

    # Generate the threshold levels for weak and strong hits
    read_weak_thresh = statistics.quantiles(
        read_norm_list, 
        n=weak_n, 
        method='inclusive'
        )[weak_i]
    read_strong_thresh = statistics.quantiles(
        read_norm_list, 
        n=strong_n, 
        method='inclusive'
        )[strong_i]
    multi_weak_thresh = statistics.quantiles(
        multi_list, 
        n=weak_n, 
        method='inclusive'
        )[weak_i]
    multi_strong_thresh = statistics.quantiles(
        multi_list, 
        n=strong_n, 
        method='inclusive'
        )[strong_i]
    
    # We can bridge regions that have 0 signal as well.
    min_thresh = min(read_norm_list)

    current_thresholds = {
        'read_weak_met': False,
        'multi_weak_met': False,
        'min_thresh_met': False,
        'read_strong_met': False,
        'multi_strong_met': False,
    }
    current_region_size = 0
    bed_list = []
    for i in range(len(bins_map)):
        current_region_size += 1
        previous_thresholds = {key: value for key, value in current_thresholds.items()}
        current_thresholds['read_weak_met'] = read_norm_list[i] >= read_weak_thresh
        current_thresholds['multi_weak_met']  = multi_list[i] >= multi_weak_thresh
        current_thresholds["min_thresh_met"]  = read_norm_list[i] <= min_thresh
        
        weak_thresh_met = (
            current_thresholds['read_weak_met'] or 
            current_thresholds['multi_weak_met'] or 
            current_thresholds["min_thresh_met"]
            )
        
        if weak_thresh_met:
            
            # check strong thresholds
            current_thresholds["read_strong_met"] = read_norm_list[i] >= read_strong_thresh
            current_thresholds["multi_strong_met"] = multi_list[i] > multi_strong_thresh

            # tracking for overlapping regions
            miss = 0
            temp_last = i

            # If this is a new region, record it
            if not in_region:
                temp_first = i
                in_region = True

            # Check to see if this bin passes a threshold
            if current_thresholds["read_strong_met"]:
                temp_hit = 1
                hit_code = hit_code | 1
                hit_counter += 1

            if current_thresholds["multi_strong_met"]:
                temp_hit = 1
                hit_code = hit_code | 2
                hit_counter += 1
            
        else:
            if miss < (bin_size // bin_overlap) + bridge_constant: # bridge over adjacent bins plus 100 * 200 = 20kb
                miss += 1                           # recommend 5k for the smaller genomes
            else: # nothing in this distance
                in_region = False

                # If we hit a threshold we output the whole region
                if temp_hit == 1:
                    if hit_code == 2:
                        region_class = "Low Mappability"
                    else:
                        region_class = "High Signal Region"
                    
                    if output_original or original_only:
                        if output_original and not original_only:
                            region_class = f'*{region_class}*'
                        bed_list.append(
                            (region, temp_first*bin_overlap, temp_last*bin_overlap + bin_size, region_class)
                            )
                        
                    temp_hit = 0
                    hit_code = 0
                    hit_counter = 0
        
        # new outputs
        region_change = previous_thresholds != current_thresholds

        if (region_change or i == len(bins_map) - 1) and not original_only:
            
            def get_region_class(thresholds: dict[str: bool]) -> str:
                region_class = ''
                if thresholds['min_thresh_met']:
                    region_class += 'Min Threshold Met, '
                if thresholds['read_strong_met']:
                    region_class += 'Very High Signal, '
                elif thresholds['read_weak_met']:
                    region_class += 'High Signal, '
                if thresholds['multi_strong_met']:
                    region_class += 'Very Low Mappability, '
                elif thresholds['multi_weak_met']:
                    region_class += 'Low Mappability, '
                return region_class

            # end of file
            if i == len(bins_map) - 1:
                
                # get region class
                region_class = get_region_class(current_thresholds)
                
                # if any threshold is met
                if region_class != '':
                    region_class = region_class.rstrip(', ')
                    bed_list.append(
                        (region, (i - current_region_size) * bin_overlap, i * bin_overlap + bin_size, region_class)
                        )
            
            # not end of file
            else:

                # get region class
                region_class = get_region_class(previous_thresholds)
                
                # if any threshold was met in previous region
                if region_class != '':
                    region_class = region_class.rstrip(', ')
                    bed_list.append(
                        (region, (i - current_region_size) * bin_overlap, i * bin_overlap + bin_size, region_class)
                        )
            
            # reset region size
            current_region_size = 0

    # If we were in a region when when hit the end of the chromosome, output region
    # This may go past chromosome end!!
    if temp_hit == 1:
        if hit_code == 2:
            region_class = "Low Mappability"
        else:
            region_class = "High Signal Region"
        if output_original or original_only:
            if output_original and not original_only:
                region_class = f'*{region_class}*'
            bed_list.append(
                (region, temp_first*bin_overlap, temp_last*bin_overlap + bin_size, region_class)
                )
    
    out_string = ''
    for line in bed_list:
        line: tuple
        out_string += line[0] + '\t'
        out_string += str(line[1]) + '\t'
        out_string += str(line[2]) + '\t'
        out_string += line[3] + '\n'

    return out_string.strip()

def get_args() -> list:
    parser = argparse.ArgumentParser(
        prog='blacklist2.py',
        description='''improved blacklist generation\n
					Copyright (c) 2023 Brydon Wall. \
                    Original Blacklist by Alan Boyle. \
                    This program comes with ABSOLUTELY \
                    NO WARRANTY. This is free software, \
                    and you are welcome to redistribute \
                    it under certain conditions.'''
	)
    parser.add_argument(
        '-b', '--bams',
        default='./input',
        nargs='?',
        help='A comma separated list of bam files or \
            directories containing bam files to process. \
            Default is ./input.'
    )
    parser.add_argument(
        '-m', '--mappability',
        default='./mappability',
        nargs='?',
        help='This is the directory of mappability file(s); \
            Default is ./mappability.'
    )
    parser.add_argument(
        '-r', '--regions',
        default='all',
        nargs='?',
        help='A comma separated list to tell \
            the program to look for and output only the \
            specified regions. If mappability file for \
            specific region is not found, region will be \
            skipped; Default is all.'
    )
    parser.add_argument(
        '-i', '--bin', 
        default='1000',
        nargs='?',
        help='Sets the size of bin to use; default is 1000.'
        )
    parser.add_argument(
        '-p', '--overlap',
        default='100',
        nargs='?',
        help='Sets the overlap between bins; default is 100.'
    )
    parser.add_argument(
        '-g', '--bridge',
        default='200',
        nargs='?',
        help='Sets the bridge size. Does not affect new output; \
            default is 200.'
    )
    parser.add_argument(
        '-u', '--unique',
        default='36',
        nargs='?',
        help='Sets the unique length within the get_mappability_bins \
            function. According to Boyle, "This is arbitraty and \
            defines how long a read needs to be to be considered \
            unique Should be set to something actually calculated in \
            the uint8 files."'
    )
    parser.add_argument(
        '-l', '--original',
        action='store_true',
        help='Will output original regions\' annotations marked with \
            asterisks: *region*.'
    )
    parser.add_argument(
        '-L', '--originalOnly',
        action='store_true',
        help='Will output with original algorithm from Boyle only.'
    )
    parser.add_argument(
        '-w', '--weak',
        default='0.99',
        nargs='?',
        help='Sets the weak percentile for calculating thresholds; \
            default is 0.99.'
    )
    parser.add_argument(
        '-s', '--strong',
        default='0.999',
        nargs='?',
        help='Sets the strong percentile for calculating thresholds; \
            default is 0.999.'
    )
    parser.add_argument(
        '-o', '--output',
        default='./b2output.bed',
        nargs='?',
        help='Path to output .bed file; Default is ./b2output.bed.'
    )
    parser.add_argument(
        '-v', '--view',
        action='store_true',
        help='Whether to print output to terminal and not to file.'
    )
    parser.add_argument(
        '-n', '--noMerge',
        action='store_true',
        help='Do not merge regions if they are overlaping'
    )

    args = parser.parse_args()
    
    # parse bams
    bams = args.bams.split(',')
    mapp = args.mappability
    regions = args.regions.split(',')
    bin_ = int(args.bin)
    overlap = int(args.overlap)
    bridge = int(args.bridge)
    unique = int(args.unique)

    
    # parse bools
    original = False
    original_only = False
    view = False
    no_merge = False
    if args.original:
        original = not original
    if args.originalOnly:
        original_only = not original_only
        original = True
    if args.view:
        view = not view
    if args.noMerge:
        no_merge = not no_merge

    output = args.output

    weak = float(args.weak)
    strong = float(args.strong)
    
    out_list = [
        bams,
        mapp,
        regions,
        bin_,
        overlap,
        bridge,
        unique,
        original,
        original_only,
        weak,
        strong,
        output,
        view,
        no_merge
    ]

    return out_list

def get_bams(
        bam_dirs: list[str]
        ) -> np.ndarray[SequenceData]:
    
    out_list = []
    bam_files = []
    for bam in bam_dirs:

        # if file
        if os.path.isfile(bam) and bam.endswith('.bam'):
            bam_files.append(f'{bam}')

        # if directory
        elif os.path.isdir(bam):

            try:
                for entry in os.scandir(bam):
                    if entry.is_file() and entry.name.endswith('.bam'):
                        bam_file = os.path.normpath(bam + '/' + entry.name)
                        bam_files.append(bam_file)
            except OSError:
                raise FileNotFoundError(f'Error in reading files from: {bam}')
        else:
            raise FileNotFoundError(f'{bam} not found.')
        
        for bam_file in bam_files:

            # look for missing index file
            bai = bam_file + '.bai'
            if not os.path.exists(bai):
                raise FileNotFoundError(f'BAM index {bai} not found.')
            out_list.append(SequenceData(bam_file, bai))

    return np.array(out_list, dtype=SequenceData)

def get_mappability(
        mappability_dir: str, 
        regions: str, 
        view: bool = False
        ) -> dict[str, int]:

    out_dict = {}

    for region in regions:
        
        # regex pattern
        pattern = region + r'(?!\d)' # not a digit

        # Get mappability file
        mappability_file = False
        for entry in os.scandir(mappability_dir):
            if re.search(pattern, entry.name):
                mappability_file = os.path.normpath(mappability_dir + '/' + entry.name)
        
        # No mappability found
        if not mappability_file:
            if not view:
                print(f'Could not find mappability file for "{region}" in {mappability_dir} - skipping')
        else:
            # Parse mappability to list of ints
            try:
                with open(mappability_file, "rb") as in_file:
                    mappability = np.fromfile(in_file, dtype=np.uint8)
            except:
                raise FileNotFoundError(f"Error in reading mappability file. Expected uint8 at: {mappability_file}")
            
            out_dict[region] = mappability
    
    return out_dict

def get_regions(
        regions: list[str], 
        bams: np.ndarray[SequenceData]
        ) -> set[str]:
    if len(regions) == 1 and regions[0] == 'all':
        out_set = set()
        for bam in bams:
            with pysam.AlignmentFile(bam.bam_file, 'rb') as aln:
                regions = aln.references
                for ref in regions:
                    out_set.add(ref)
        return out_set
    else:
        return regions

def merge_annotations(a: str, b: str) -> str:
    
    # Make set 
    out_set = set()
    a = a.strip().split(', ')
    b = b.strip().split(', ')
    for x in a:
        out_set.add(x)
    for y in b:
        out_set.add(y)
    
    # Combine Signals and Mappabilities
    if 'High Signal' in out_set and 'Very High Signal' in out_set:
        out_set.remove('High Signal')
    if 'Low Mappability' in out_set and 'Very Low Mappability' in out_set:
        out_set.remove('Low Mappability')
    
    out_set = list(out_set)

    # Sort while ignoring the word 'Very '
    out_set.sort(
        key=lambda a: a[5:] if a.startswith('Very ') else a
        )
    out_str = ', '.join(out_set) + '\n'
    return out_str.strip(', ')

def merge_bed(bed: str) -> str:
    
    out_str = ''
    bed_lines = bed.splitlines(keepends=True)

    merged_lines = []

    # Chop overlapping regions
    chr_dict: dict[list] = {}
    for line in bed_lines:
        line = line.split('\t')
        line[1] = int(line[1])
        line[2] = int(line[2])
        if not line[0] in chr_dict:
            chr_dict[line[0]] = [line]
        else:
            chr_dict[line[0]].append(line)
    
    chopped_set: set[int] = set()
    i = 0
    for chr in chr_dict:
        lines = chr_dict[chr]
        for line in lines:
            chopped_set.add(line[1])
            chopped_set.add(line[2])
            lines[i] = line
            i += 1
        chopped_list: list[int] = list(chopped_set)
        chopped_list.sort()
        
        # Merge

        for i in range(len(chopped_list) - 1):
            lines_to_merge = []
            start = chopped_list[i]
            stop = chopped_list[i + 1]

            # loop through 
            for line in lines:

                # break out of loop
                if line[1] > stop:
                    break
                
                # to merge found
                elif line[1] <= start and line[2] >= stop:
                    lines_to_merge.append(line)
                
                # remove lines that are no longer needed
            
            # merge annotations
            annotation = ''
            for line in lines_to_merge:
                annotation = merge_annotations(annotation, line[3])
            
            if len(lines_to_merge) >= 1:
                merged_lines.append(
                    [
                        lines_to_merge[0][0],
                        start,
                        stop,
                        annotation
                    ]
                )

        # merge touching 
        for i in range(len(merged_lines) - 1):
            line1 = merged_lines[i]
            line2 = merged_lines[i + 1]
            if line1[2] == line2[1] and line1[3] == line2[3]:
                line2[1] = line1[1]
                line1[3] = 'Merged\n'

        for line in merged_lines:
            if line[3] != 'Merged\n':
                out_str += '\t'.join([str(item) for item in line])

    return out_str.strip()

def main():

    # Parameters
    args = get_args()
    bams: list[str] = args[0]
    mappability_dir: str = args[1]
    regions: list[str] = args[2]
    bin_size: int = args[3]
    bin_overlap: int = args[4]
    bridge_size: int = args[5]
    unique_length: int = args[6]
    output_original: bool = args[7]
    original_only: bool = args[8]
    weak_percentile: float = args[9]
    strong_percentile: float = args[10]
    output_path: str = args[11]
    view: bool = args[12]
    no_merge: bool = args[13]
    out_str = ''

    # Perform Calculations per BAM and region
    bam_data = get_bams(bams)

    # Process regions
    regions = get_regions(regions, bam_data)

    # Process mappability files
    mappability_dict = get_mappability(mappability_dir, regions, view)

    # Loop through regions
    for region in regions:
        
        # Get mappability
        try:
            mappability = mappability_dict[region]
        
        # Skip region
        except KeyError:
            continue

        # Get input bins
        for bam in bam_data:
            bam.get_input_bins_per_region(
                mappability,
                bin_size,
                bin_overlap,
                region
            )

        bins_map = get_mappability_bins(mappability, bin_size, bin_overlap, unique_length)

        read_norm_list = get_abnormal_regions(bam_data, bins_map, 'reads', True)
        multi_list = get_abnormal_regions(bam_data, bins_map, 'multi', False)

        # get .bed
        bed = get_bed(
            read_norm_list,
            multi_list,
            bins_map,
            bin_size,
            bin_overlap,
            region,
            bridge_size,
            output_original,
            original_only,
            weak_percentile,
            strong_percentile
            )
        if not no_merge:
            bed = merge_bed(bed)

        out_str += bed

    # write output
    if not view:
        with open(output_path, 'w') as out_file:
            out_file.write(out_str)
            print(f'Blacklist written to {output_path}')

    # print
    if view:
        print(out_str)

if __name__ == '__main__':
    main()
