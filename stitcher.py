#!/usr/bin/env python
# V 3.0
# anton jm larsson anton.larsson@basic-genomics.com
import argparse
import pysam
import warnings
import numpy as np
import pygtrie
import portion as P
import itertools
import sys
import time
import os
import json
from scipy.special import logsumexp
from joblib import delayed,Parallel
from multiprocessing import Process, Manager
from collections import Counter
import faulthandler
__version__ = '3.0'
nucleotides = ['A', 'T', 'C', 'G']
nuc_dict = {'A':0, 'T':1, 'C':2, 'G':3, 'N': 4}
np.seterr(divide='ignore')
ll_this_correct = {i:np.log(1-10**(-float(i)/10)) for i in range(1,94)}
ln_3 = np.log(3)
ll_other_correct = {i:-(float(i)*np.log(10))/10 - ln_3 for i in range(1,94)}
ll_N = -np.log(4)

from scipy.sparse import csc_matrix
def make_ll_array(e):
    y = np.array([e[0]/3,e[0]/3,e[0]/3,e[0]/3])
    if e[1] != 4:
        y[e[1]] = 1-e[0]
    return np.log10(y)

# taken from https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def intervals_extract(iterable): 
    iterable = sorted(set(iterable)) 
    for key, group in itertools.groupby(enumerate(iterable), 
    lambda t: t[1] - t[0]): 
        group = list(group) 
        yield [group[0][1], group[-1][1]] 

def interval(t):
    return P.from_data([(True,i[0],i[1], True) for i in t])

def get_time_formatted(time):
    day = time // (24 * 3600)
    time = time % (24 * 3600)
    hour = time // 3600
    time %= 3600
    minutes = time // 60
    time %= 60
    seconds = time
    s = ''.join(['{} day{}, '.format(day, 's'*(1 != day))*(0 != day), 
                 '{} hour{}, '.format(hour,'s'*(1 != hour))*(0 != hour), 
                 '{} minute{}, '.format(minutes,'s'*(1 != minutes))*(0 != minutes), 
                 '{:.2f} second{}, '.format(seconds,'s'*(1 != seconds))*(0 != seconds)])
    s = s[:-2]
    s = s + '.'
    return s

def get_insertions_locs(cigtuples):
    insertion_locs = []
    l = 0
    for c in cigtuples:
        if c[0] == 0:
            l += c[1]
        elif c[0] == 1:
            for i in range(c[1]):
                insertion_locs.append(l)
                l += 1
    return insertion_locs

def get_skipped_tuples(cigtuples, ref_positions):
    skipped_locs = []
    l = 0
    for c in cigtuples:
        if c[0] == 0:
            l += c[1]
        elif c[0] == 3:
            skipped_locs.append((ref_positions[l-1]+1, ref_positions[l]-1))
    return skipped_locs

def using_indexed_assignment(x):
    "https://stackoverflow.com/a/5284703/190597 (Sven Marnach)"
    result = np.empty(len(x), dtype=int)
    temp = x.argsort()
    result[temp] = np.arange(len(x))
    return result



def stitch_reads(read_d, cell, gene, umi, UMI_tag):
    faulthandler.enable()
    master_read = {}
    nreads = len(read_d)
    exonic_list = [0]*nreads
    intronic_list = [0]*nreads
    orientation_counts = {'+': 0, '-': 0, 'NA': 0, 'no': 0}
    read_type_counts = {'TP_read': 0, 'internal': 0, 'FP_read': 0}
    seq_list = []
    qual_list = []
    ref_pos_set = set()
    ref_pos_list = []
    threeprime_start = Counter()
    fiveprime_start = Counter()
    reference_pos_counter = Counter()
    for i,read in enumerate(read_d):
        if read.has_tag('GE'):
            exonic = True
        else:
            exonic = False
        if read.has_tag('GI'):
            intronic = True
        else:
            intronic = False
        try:
            #with open('debug_output.txt', 'w') as f:
                #f.write(read+'\n')
            Q_list = list(read.query_alignment_qualities)
        except TypeError:
            Q_list = [read.query_alignment_qualities]
        
        seq = read.query_alignment_sequence
        cigtuples = read.cigartuples
        insertion_locs = get_insertions_locs(cigtuples)
        
        if len(insertion_locs) > 0:
            seq = "".join([char for idx, char in enumerate(seq) if idx not in insertion_locs])
            Q_list = [qual for idx, qual in enumerate(Q_list) if idx not in insertion_locs]

        ref_positions = read.get_reference_positions()
        skipped_intervals = get_skipped_tuples(cigtuples, ref_positions)

        exonic_list[i] = exonic
        intronic_list[i] = intronic

        seq_list.append(seq)

        qual_list.append(Q_list)

        ref_pos_list.append(ref_positions)
        reference_pos_counter.update(ref_positions)

        ref_pos_set = ref_pos_set | set(ref_positions)

        orientation = read.get_tag('YZ')
        orientation_counts[orientation] += 1
        
        read_type = read.get_tag('XX')
        read_type_counts[read_type] += 1

        if orientation == '+' and read_type == 'TP_read' and read.is_read1 and read.is_reverse:
            threeprime_start.update({read.reference_end: 1})
        if orientation == '-' and read_type == 'TP_read' and read.is_read1 and not read.is_reverse:
            threeprime_start.update({read.reference_start: 1})
        
        if orientation == '+' and read_type == 'FP_read' and read.is_read1 and not read.is_reverse:
            fiveprime_start.update({read.reference_start: 1})
        if orientation == '-' and read_type == 'FP_read' and read.is_read1 and read.is_reverse:
            fiveprime_start.update({read.reference_end: 1})

        if len(master_read) == 0:
            master_read['skipped_interval_list'] = [skipped_intervals]
        else:
            master_read['skipped_interval_list'].append(skipped_intervals)

    sparse_row_dict = {b:[] for b in nucleotides}
    sparse_col_dict = {b:[] for b in nucleotides}
    sparse_ll_dict = {b:[] for b in nucleotides}

    ref_pos_set_array = np.array(list({p for p in ref_pos_set}))

    master_read['ref_intervals'] = interval(intervals_extract(np.sort(ref_pos_set_array)))
    master_read['ref_pos_counter'] = reference_pos_counter
    master_read['skipped_intervals'] = interval(list(set([item for sublist in master_read['skipped_interval_list'] for item in sublist])))


    ref_and_skip_intersect = master_read['ref_intervals'] & master_read['skipped_intervals']
    reference_positions = []
    skipped_positions = []

    if not ref_and_skip_intersect.empty:
        conflict_pos_list = P.iterate(ref_and_skip_intersect, step=1)
        skip_pos_counter = Counter()
        for skip_tuples in master_read['skipped_interval_list']:
            skip_interval = interval(skip_tuples) 
            skip_pos_counter.update([p for p in conflict_pos_list if p in skip_interval])
        for pos in conflict_pos_list:
            if master_read['ref_pos_counter'][pos] > skip_pos_counter[pos]:
                reference_positions.extend(pos)
            else:
                skipped_positions.extend(pos)

        reference_keep_intervals = interval(intervals_extract(reference_positions))
        skip_keep_intervals = interval(intervals_extract(skipped_positions))
        ### No refskip where there is reference sequence
        master_read['skipped_intervals'] = master_read['skipped_intervals'] - reference_keep_intervals
        ### No ref coverage where there is refskip
        master_read['ref_intervals'] = master_read['ref_intervals'] - skip_keep_intervals

    master_read['skipped_intervals'] = master_read['skipped_intervals'] & P.closed(master_read['ref_intervals'].lower, master_read['ref_intervals'].upper)

    ref_pos_set_array = np.array(list({p for p in P.iterate(master_read['ref_intervals'], step=1)}))

    if len(ref_pos_set_array) == 0:
        return (False, ':'.join([gene,cell,umi]))

    ref_to_pos_dict = {p:o for p,o in zip(ref_pos_set_array, using_indexed_assignment(ref_pos_set_array))}

    for i, (seq, Q_list, ref_positions) in enumerate(zip(seq_list, qual_list, ref_pos_list)):
        for b1, Q, pos in zip(seq,Q_list, ref_positions):
            if pos not in master_read['ref_intervals']:
                continue
            for b2 in nucleotides:
                sparse_row_dict[b2].append(i)
                sparse_col_dict[b2].append(ref_to_pos_dict[pos])
                if b1 == b2:
                    sparse_ll_dict[b2].append(ll_this_correct[Q])
                elif b1 == 'N':
                    sparse_ll_dict[b2].append(ll_N)
                else:
                    sparse_ll_dict[b2].append(ll_other_correct[Q])
    sparse_csc_dict = {b:csc_matrix((sparse_ll_dict[b], (sparse_row_dict[b],sparse_col_dict[b])), shape=(i+1,len(ref_pos_set_array))) for b in nucleotides}

    ll_list = [m.sum(axis=0) for m in sparse_csc_dict.values()]
    try:
        ll_sums = np.stack(ll_list)
    except ValueError:
        return (False, ':'.join([gene,cell,umi]))
    full_ll = logsumexp(ll_sums, axis=0)

    prob_max = np.asarray(np.exp(np.amax(ll_sums, axis=0) - full_ll)).ravel()
    nuc_max = np.asarray(np.argmax(ll_sums, axis=0)).ravel()

    master_read['seq'] = ''.join([nucleotides[x] if p > 0.3 else 'N' for p, x in zip(prob_max, nuc_max)])
    master_read['phred'] = np.nan_to_num(np.rint(-10*np.log10(1-prob_max+1e-13)))
    master_read['SN'] = read.reference_name
    master_read['is_reverse'] = 1 if orientation_counts['-'] >= orientation_counts['+'] else 0
    master_read['del_intervals'] =  ~(master_read['ref_intervals'] | master_read['skipped_intervals'])
    master_read['NR'] = nreads
    master_read['IR'] = np.sum(intronic_list)
    master_read['ER'] = np.sum(exonic_list)
    master_read['TC'] = read_type_counts['TP_read']
    master_read['IC'] = read_type_counts['internal']
    master_read['FC'] = read_type_counts['FP_read']
    master_read['cell'] = cell
    master_read['gene'] = gene
    master_read['umi'] = umi
    sam = convert_to_sam(master_read, UMI_tag)
    if sam is None:
        return (False, ':'.join([gene,cell,umi]))
    else:
        return (True, sam)

def assemble_reads(bamfile,gene_to_stitch, cell_set, cell_tag, UMI_tag, only_molecules, q):
    readtrie = pygtrie.StringTrie()
    bam = pysam.AlignmentFile(bamfile, 'rb')
    gene_of_interest = gene_to_stitch['gene_id']
    for read in bam.fetch(gene_to_stitch['seqid'], gene_to_stitch['start'], gene_to_stitch['end']):
        if read.has_tag(cell_tag):
            cell = read.get_tag(cell_tag)
        else:
            continue
        if cell_set is not None:
            if cell not in cell_set:
                continue
        if not read.has_tag(UMI_tag):
            continue
        umi = read.get_tag(UMI_tag)
        if umi == '':
            continue
        else:
            if only_molecules:
                if umi[0] == '_':
                    continue
                if len(umi) >= 5:
                    if umi[:5] == 'Unass':
                        continue
            if read.has_tag('GE'):
                gene_exon = read.get_tag('GE')
            else:
                gene_exon = 'Unassigned'
            if read.has_tag('GI'):
                gene_intron = read.get_tag('GI')
            else:
                gene_intron = 'Unassigned'
            # if it maps to the intron or exon of a gene
            if gene_intron != 'Unassigned' or gene_exon != 'Unassigned':
                # if it is a junction read
                if gene_intron == gene_exon:
                    gene = gene_intron
                    # if it's an only intronic read
                elif gene_intron != 'Unassigned' and gene_exon == 'Unassigned':
                    gene = gene_intron
                    # if it's an only exonic read
                elif gene_exon != 'Unassigned' and gene_intron == 'Unassigned':
                    gene = gene_exon
                    # if the exon and intron gene tag contradict each other
                else:
                    continue
            else:
                continue
        

        if gene == gene_of_interest:
            node = '{}/{}/{}'.format(cell,gene,umi)
            if readtrie.has_node(node):
                readtrie[node].append(read)
            else:
                readtrie[node] = [read]
    mol_list = []
    mol_append = mol_list.append
    for node, mol in readtrie.iteritems():
        info = node.split('/')
        mol_append(stitch_reads(mol, info[0], info[1], info[2], UMI_tag))
    del readtrie
    if len(mol_list) == 0:
        return gene_of_interest
    if len(mol_list) > 50000:
        for m_list in chunks(mol_list, 50000):
            q.put((True, m_list))
    else:
        q.put((True, mol_list))
    return gene_of_interest


def make_POS_and_CIGAR(stitched_m):
    CIGAR = ''
    ref_tuples = [(i[1] if i[0] else i[1]+1, i[2] if i[3] else i[2]-1) for i in P.to_data(stitched_m['ref_intervals'])]
    if stitched_m['skipped_intervals'].empty:
        skipped_tuples = []
    else:
        skipped_tuples = [(i[1] if i[0] else i[1]+1, i[2] if i[3] else i[2]-1) for i in P.to_data(stitched_m['skipped_intervals'])]
    if stitched_m['del_intervals'].empty:
        del_tuples = []
    else:
        del_tuples = [(i[1] if i[0] else i[1]+1, i[2] if i[3] else i[2]-1) for i in P.to_data(stitched_m['del_intervals'])[1:-1]]
    POS = ref_tuples[0][0] + 1
    tuple_dict = {'M': ref_tuples, 'N': skipped_tuples, 'D': del_tuples}
    while sum(len(t) for t in tuple_dict.values()) > 0:
        pos_dict = {k:v[0][0] for k,v in tuple_dict.items() if len(v) > 0}
        c = min(pos_dict, key=pos_dict.get)
        n_bases = np.int_(tuple_dict[c[0]][0][1]-tuple_dict[c[0]][0][0])+1
        if n_bases == 0:
            del tuple_dict[c[0]][0]
            continue
        CIGAR += '{}{}'.format(n_bases,c[0])
        del tuple_dict[c[0]][0]
    return POS, CIGAR

def get_len_from_cigar(cigarstring):
    total_len = 0
    len_string = ''
    for s in cigarstring:
        if s not in ['M', 'N', 'D']:
            len_string += s
        elif s == 'M':
            total_len += int(len_string)
            len_string = ''
        else:
            len_string = ''
    return total_len

def convert_to_sam(stitched_m, UMI_tag):
    sam_dict = {}
    POS, CIGAR = make_POS_and_CIGAR(stitched_m)
    len_from_cigar = get_len_from_cigar(CIGAR)
    if len_from_cigar != len(stitched_m['seq']):
        return None
    sam_dict['QNAME'] = '{}:{}:{}'.format(stitched_m['cell'],stitched_m['gene'],stitched_m['umi'])
    sam_dict['FLAG'] = str(16*stitched_m['is_reverse'])
    sam_dict['RNAME'] = stitched_m['SN']
    sam_dict['POS'] = str(POS)
    sam_dict['MAPQ'] = str(255)
    sam_dict['CIGAR'] = CIGAR
    sam_dict['RNEXT'] = '*'
    sam_dict['PNEXT'] = str(0)
    sam_dict['TLEN'] = str(0)
    sam_dict['SEQ'] = stitched_m['seq']
    sam_dict['QUAL'] = "".join([chr(int(p)) for p in np.clip(stitched_m['phred'],0,126-33)+33])
    sam_dict['NR'] = 'NR:i:{}'.format(stitched_m['NR'])
    sam_dict['ER'] = 'ER:i:{}'.format(stitched_m['ER'])
    sam_dict['IR'] = 'IR:i:{}'.format(stitched_m['IR'])
    sam_dict['TC'] = 'TC:i:{}'.format(stitched_m['TC'])
    sam_dict['IC'] = 'IC:i:{}'.format(stitched_m['IC'])
    sam_dict['FC'] = 'FC:i:{}'.format(stitched_m['FC'])
    sam_dict['BC'] = 'BC:Z:{}'.format(stitched_m['cell'])
    sam_dict['XT'] = 'XT:Z:{}'.format(stitched_m['gene'])
    sam_dict[UMI_tag] = '{}:Z:{}'.format(UMI_tag, stitched_m['umi'])
    return '\t'.join(list(sam_dict.values()))

def yield_reads(read_dict):
    for cell in read_dict:
        for gene in read_dict[cell]:
            #print('\t', gene)
            for umi in read_dict[cell][gene]:
                #print('\t\t', umi)
                yield read_dict[cell][gene][umi], None, cell, gene, umi

def get_length_from_cigar(read):
    total_len = 0
    for (op, len) in read.cigartuples:
        if op == 0:
            total_len += len
    return total_len
def create_write_function(filename, bamfile, version):
    bam = pysam.AlignmentFile(bamfile, 'rb')
    header = bam.header
    
    def write_sam_file(q):
        error_file = open('{}_error.log'.format(os.path.splitext(filename)[0]), 'w')
        stitcher_bam = pysam.AlignmentFile(filename,'wb',header={'HD':header['HD'], 'SQ':header['SQ'], 'PG': [{'ID': 'stitcher.py','VN': '{}'.format(version)}]})
        while True:
            good, mol_list = q.get()
            if good is None: break
            if good:
                g = ''
                for success, mol in mol_list:
                    if success:
                        read = pysam.AlignedRead.fromstring(mol,header)
                        if g == '':
                            g = read.get_tag('XT')
                        stitcher_bam.write(read)
                    else:
                        error_file.write(mol+'\n')
                if g != '':
                    error_file.write('Gene:{}\n'.format(g))
            q.task_done()
        q.task_done()
        error_file.close()
        stitcher_bam.close()
        return None
    return write_sam_file

def extract(d, keys):
    return dict((k, d[k]) for k in keys if k in d)
    
def construct_stitched_molecules(infile, gtffile, cells, gene_file, contig, threads, cell_tag, UMI_tag, gene_identifier,only_molecules, q, version):
    if cells is not None:
        cell_set = set([line.rstrip() for line in open(cells)])
    else:
        cell_set = None
    print('Reading gene info from {}'.format(gtffile))
    gene_list = []
    if gene_identifier == 'gene_id':
        n = 1
    elif gene_identifier == 'gene_name':
        n = 5
    else:
        n = 1
    with open(gtffile, 'r') as f:
        for line in f:
            l = line.split('\t')
            if len(l) < 8:
                continue
            if l[2] == 'gene':
                if contig is not None:
                    if l[0] == contig:
                        try:
                            gene_list.append({'gene_id': l[8].split(' ')[n].replace('"', '').strip(';\n'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4])})
                        except:
                            gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';\n'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4])})
                    else:
                        continue
                else:
                    try:
                        gene_list.append({'gene_id': l[8].split(' ')[n].replace('"', '').strip(';\n'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4])})
                    except:
                        gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';\n'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4])})
    gene_dict = {g['gene_id']: g for g in gene_list}
    
    if gene_file is not None and gene_file != 'None':
        gene_set = set([line.rstrip() for line in open(gene_file)])
        gene_dict = {k:v for k,v in gene_dict.items() if k in gene_set}
    bam = pysam.AlignmentFile(infile, 'rb')
    contig_set = set([d['SN'] for d in bam.header['SQ']])
    prev_l = len(gene_dict)
    gene_dict = {k:v for k,v in gene_dict.items() if v['seqid'] in contig_set}
    new_l = len(gene_dict)
    diff_l = prev_l - new_l
    print(prev_l, new_l)
    if diff_l > 0:
        warnings.warn('Warning: removed {diff_l} genes with contig not present in bam file'.format(diff_l=diff_l))
    bam.close()
    
    
    res = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(assemble_reads)(infile, gene, cell_set, cell_tag, UMI_tag, only_molecules, q) for g,gene in gene_dict.items())


    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stitch together molecules from reads sharing the same UMI')
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output .bam file')
    parser.add_argument('-g','--gtf', metavar='gtf', type=str, help='gtf file with gene information')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--UMI-tag', type=str, default='UB', help='UMI tag to group reads')
    parser.add_argument('--cell-tag', type=str, default='BC', help='cell barcode tag to group reads')
    parser.add_argument('--cells', default=None, metavar='cells', type=str, help='List of cell barcodes to stitch molecules')
    parser.add_argument('--genes', default=None, metavar='genes', type=str, help='List of gene,  one per line.')
    parser.add_argument('--contig', default=None, metavar='contig', type=str, help='Restrict stitching to contig')
    parser.add_argument('--gene-identifier', default='gene_id', metavar='gene_identifier', type=str, help='Gene identifier (gene_id or gene_name)')
    parser.add_argument('--only-molecules', action='store_true', help='Only reconstruct molecues')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()
    infile = args.input
    if infile is None:
        raise Exception('No input file provided.')
    outfile = args.output
    if outfile is None:
        raise Exception('No output file provided.')
    gtffile = args.gtf  
    if gtffile is None:
        raise Exception('No gtf file provided.')
    threads = int(args.threads)
    cells = args.cells
    gene_file = args.genes
    contig = args.contig
    UMI_tag = args.UMI_tag
    cell_tag = args.cell_tag
    gene_identifier = args.gene_identifier
    only_molecules = args.only_molecules
    m = Manager()
    q = m.JoinableQueue()
    p = Process(target=create_write_function(filename=outfile, bamfile=infile, version=__version__), args=(q,))
    p.start()
    
    print('Stitching reads for {}'.format(infile))
    
    start = time.time()
    construct_stitched_molecules(infile, gtffile, cells, gene_file, contig, threads,cell_tag, UMI_tag,gene_identifier, only_molecules, q, __version__ )
    q.put((None,None))
    p.join()
    end = time.time()
    
    print('Finished writing stitched molecules from {} to {}, took {}'.format(infile, outfile, get_time_formatted(end-start)))
