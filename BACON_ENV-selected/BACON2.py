'''
BArcode COuNter (BACON) v0.1.1
DEVELOPED BY ANDREW CHINN
THANKS TO THE LAB OF DR. PAMELA HANSON
AT FURMAN UNIVERSITY, SC, USA

***
BASED ON BARCOUNT v0.92
DEVELOPED BY stephan.kamrad@crick.ac.uk
FOUND AT https://github.com/Bahler-Lab/barcount
***

FOR QUESTIONS OR HELP CONTACT
chinan3@furman.edu

***
TO DO:
- Make log (fix it lol)
- Account for mutations (check but not with primers)
- Add exception errors (please)
- Formalize i/o process (make better output stats)
- Deicde on a new name
'''

#IMPORT LIBRARIES
import ntpath
import sys
import argparse
import configparser
import re
from warnings import warn
import pandas as ps
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.CheckSum import seguid
from time import sleep
from tqdm import tqdm
from pathlib import Path
import multiprocessing
import math
import glob
import os
import Levenshtein

#MAIN FUNCTION
def run_bacon(barcode_table, output, fastq_input, FWD_Primer, REV_Primer, barcode_table_alt, barcode_rev_comp, output_name, position, lock, max_read_length, max_barcode_length, max_lev_distance):

    #INITIALIZE BARCODE KEY FILES
    barcode_data = ps.read_csv(barcode_table, header=None, names=['gene',], index_col=1)
    barcode_data.loc['unassigned'] = np.nan

    barcode_data_alt = ps.read_csv(barcode_table_alt, header=None, names=['gene_alt',], index_col=1)
    barcode_data_alt.loc['unassigned_alt'] = np.nan

    #INITIALIZE OUTPUT STATS
    stats = ps.Series()
    stats['fastq_entries'] = 0
    stats['primer_lev_match'] = 0
    stats['matched_norm'] = 0
    stats['matched_alt'] = 0
    stats['matched_lev'] = 0
    stats['unmatched'] = 0
    stats['percent_matched'] = 0
    stats['failed_max_read'] = 0
    stats['failed_max_barcode'] = 0
    stats['percent_reads_lost'] = 0
    barcode_data['count'] = 0

    #INITIALIZE DEBUG LOG
    debug_file = open(output+'{}_debug.csv'.format(output_name), 'w')
    debug_log = ['fastq_entry',#0
                'fastq_header',#1
                'sequence',#2
                'sequence_length',#3
                'failed_max_read',#4
                'fwd_primer_input',#5
                'fwd_primer_mut',#6
                'rev_primer_input',#7
                'rev_primer_mut',#8
                'right_matching',#9
                'left_matching',#10
                'barcode_region',#11
                'barcode_matched',#12
                'match_type',#13
                'matching_quality',#14
                'matched_gene',#15
                ]

    #INITIALIZE MATCHING QUALITY TRACKERS
    #matchingqs = {bc : [] for bc in barcode_data.index}
    #barcodelengths = {bc : [] for bc in barcode_data.index}

    #OPEN FASTQ FILE
    read_seq = open(fastq_input, 'r')

    #CREATE LIST OF BARCODES FROM KEYS
    barcodes_aslist = barcode_data.index.tolist()
    barcodes_alt_aslist = barcode_data_alt.index.tolist()

    #INITIALIZE READ SEQUENCE COUNTER AT 0
    seq_count = 0

    #SET TOTAL NUMBER OF READS FOR PROGRESS BAR
    with open(fastq_input, 'r') as f:
        reads = []
        for line in f:
            if line.startswith("@"):
                reads.append(line)
    total_reads = len(reads)

    #OUTPUT_NAME REASSIGN
    #if math.isnan(output_name) == True:
        #output_name = str(ntpath.basename(fastq_input).split('.fastq')[0])

    #TQDM STUFF WITH MULTIPROCESSING
    with lock:
        bar = tqdm(
            desc=f'Process {output_name}',
            total=total_reads,
            position=position,
            leave=False
        )

    #MAIN LOOP (RUNS ONCE FOR EVERY READ/SEQUENCE IN FASTQ)
    for r in SeqIO.parse(read_seq, 'fastq'):
        seq_count += 1

        debug_file.write(','.join(map(str, debug_log))+'\n')
        debug_log = ['NA']*16


        debug_log[0] = str(seq_count)
        debug_log[1] = str(r.description)
        debug_log[2] = str(r.seq)
        debug_log[3] = len(str(r.seq))

        #TQDM STUFF WITH MULTIPROCESSING
        with lock:
            bar.update(1)

        #READLENGTH FILTER
        if max_read_length and (len(r.seq) > max_read_length):
            stats['failed_max_read'] += 1
            debug_log[4] = str(True)
            continue
        debug_log[4] = str(False)

        #ASSIGN VAR "barcode_temp" AS OUTPUT FROM find_target_region FUNCTION
        barcode_temp, mxb, debug_log[8] = find_target_region(str(r.seq), str(FWD_Primer), str(REV_Primer), max_barcode_length)
        stats['failed_max_barcode'] += mxb
        debug_log[7] = str(REV_Primer)
        debug_log[5] = str(FWD_Primer)
        debug_log[11] = str(barcode_temp)
        
        if barcode_rev_comp == True:
            barcode_temp = str(Seq(barcode_temp).reverse_complement())
        else:
            pass

        matched_barcode, key_type = match_barcode(barcode_temp, barcodes_aslist, barcodes_alt_aslist)

        debug_log[12] = str(matched_barcode)
        debug_log[13] = str(key_type)
        debug_log[14] = 'NA'

        if key_type == 'norm':
            gene = barcode_data.at[matched_barcode,'gene']
            #barcode_data.loc[matched_barcode, 'count'] += 1
            barcode_data.loc[barcode_data.gene == gene, 'count'] += 1
            stats['matched_norm'] += 1
            debug_log[15] = str(gene)
            continue
        elif key_type == 'alt':
            gene = barcode_data_alt.at[matched_barcode,'gene_alt']
            barcode_data.loc[barcode_data.gene == gene, 'count'] += 1
            stats['matched_alt'] += 1
            debug_log[15] = str(gene)
            continue
        elif key_type == 'lev':
            gene = barcode_data.at[matched_barcode, 'gene']
            barcode_data.loc[barcode_data.gene == gene, 'count'] += 1
            stats['matched_lev'] += 1
            debug_log[15] = str(gene)
            continue
        elif key_type == 'none':
            stats['unmatched'] += 1
            debug_log[15] = 'NA'

            #barcodelengths[matched_barcode].append(len(barcode))

    #EXPORT DATA TO CSV AND SORT BY MOST HITS
    debug_file.write(','.join(map(str, debug_log))+'\n')
    
    stats['fastq_entries'] = seq_count
    stats['percent_matched'] = ((stats['matched_norm'] + stats['matched_alt'] + stats['matched_lev']) / (stats['fastq_entries'] - (stats['failed_max_read'] + stats['failed_max_barcode']))) * 100
    stats['percent_reads_lost'] = ((stats['failed_max_read'] + stats['failed_max_barcode']) / (stats['fastq_entries'])) * 100

    if not os.path.isdir(output+'stats/'):
        os.makedirs(output+'stats/')

    if not os.path.isdir(output+'counts/'):
        os.makedirs(output+'counts/')

    barcode_data = barcode_data.sort_values('count', ascending=False)
    barcode_data.to_csv(output+'counts/count_{}.csv'.format(output_name))

    read_seq.close()
    debug_file.close()
    stats.to_csv(output+'stats/stats_{}_SINGLE.csv'.format(output_name))

    #RETURN FILES
    #return stats
    #with lock:
        #bar.close()

#FUNCTION TO FIND REGION PRECEEDING REVERSE PRIMER AND FOLLOWING FORWARD PRIMER (IF PRESENT)
#CURRENTLY ONLY SUPPORT LEV DISTANCES FOR SUBSTITUTION MUTATIONS
def find_target_region(input_seq, leftmost, rightmost, max_region_len):
    p_levmatch = 'FALSE'
    if rightmost != 'nan' and rightmost in input_seq:
        target = input_seq.split(rightmost)[0]
    else:
        test_reg = input_seq[-len(rightmost):]
        tdistance = Levenshtein.distance(test_reg, rightmost)
        if tdistance > max_lev_distance:
            target = input_seq.split(test_reg)[0]
            p_levmatch = 'TRUE'
        else:
            return 'failedfindrevprimer', 1, p_levmatch
    if leftmost != 'nan' and leftmost == target[:len(leftmost)]:
        target = target.split(leftmost, 1)[-1]
    else:
        pass

    if len(target) >= max_region_len:
        return 'failedmaxregionlen', 1, p_levmatch
    else:
        return target, 0, p_levmatch


#FUNCTION TO MATCH A POTENTIAL BARCODE SEQUENCE TO THE BARCODE KEY TABLES
def match_barcode(barcode_region, barcode_key, barcode_key_alt):
    if barcode_region in barcode_key:
        return barcode_region, 'norm'

    elif barcode_region in barcode_key_alt:
        return barcode_region, 'alt'

    else:
        distances = [Levenshtein.distance(barcode_region, str(v)) for v in barcode_key]
        best = min(distances)
        #check if there's a tie or the distance is too great
        if best > max_lev_distance or distances.count(best) > 1:
            return 'unassigned', 'none'
        else:
            #otherwise, return the best one
            return barcode_key[distances.index(best)], 'lev'

#CONFIG FILE SETUP
script_path = Path(__file__).parent.resolve()

config = configparser.RawConfigParser()
config.read('{}/bacon_config.cfg'.format(script_path))

config_details = dict(config.items('bacon_config'))

#GET DICT VARS
barcode_table_up = config_details['uptag_barcodes_path']
barcode_table_up_alt = config_details['uptag_alt_barcodes_path']

barcode_table_dn = config_details['dntag_barcodes_path']
barcode_table_dn_alt = config_details['dntag_alt_barcodes_path']

output = config_details['output_path']
fastq_inputs_path = config_details['fastq_inputs_path']
max_read_length = int(config_details['max_read_length'])
max_barcode_length = int(config_details['max_barcode_length'])
max_lev_distance = int(config_details['max_levenshtein_distance'])

fastq_input_table = ps.read_csv(fastq_inputs_path)

#RUN BACON
if __name__ == '__main__': 

    jobs = []
    lock = multiprocessing.Manager().Lock()

    for b in range(len(fastq_input_table)):
        if fastq_input_table.loc[b, 'up_dn'] == 'up':
            barcode_table = barcode_table_up
            barcode_table_alt = barcode_table_up_alt
        elif fastq_input_table.loc[b, 'up_dn'] == 'dn':
            barcode_table = barcode_table_dn
            barcode_table_alt = barcode_table_dn_alt
        else:
            print('Error: no up nor dn detected for fastq input row {} in the up_dn input column.'.format(b))
            raise SystemExit(0)

        fastq_input = fastq_input_table.loc[b, 'fastq_path']
        FWD_Primer = fastq_input_table.loc[b, 'fwd_primer']
        REV_Primer = str(Seq(fastq_input_table.loc[b, 'rev_primer']).reverse_complement())
        barcode_rev_comp = fastq_input_table.loc[b, 'rev_comp_bc']
        output_name = fastq_input_table.loc[b, 'fastq_cond']
        run = fastq_input_table.loc[b, 'run']

        if run == 'y':
            process = multiprocessing.Process(target=run_bacon, args=(barcode_table, output, fastq_input, FWD_Primer, REV_Primer, barcode_table_alt, barcode_rev_comp, output_name, b+1, lock, max_read_length, max_barcode_length, max_lev_distance))
            jobs.append(process)

    for j in jobs:
        j.start()

    for j in jobs:
        j.join()

    #MASTER STATS LIST
    stats_path = r'{}stats/'.format(output)
    stats_list = glob.glob(os.path.join(stats_path, "*_SINGLE.csv"))

    master_stats = ps.concat((ps.read_csv(f) for f in stats_list), keys=stats_list, axis=1)

    master_stats = master_stats[master_stats.columns[master_stats.nunique() > 1]]

    master_stats.to_csv(output+'stats/master_stats.csv')
