"""
Adjust GWAS summary statistics.
"""
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# import matplotlib

import pandas
import scipy as sp
import gzip
import pylab
import itertools as it

import h5py

def get_sid_pos_map(sids):
    sids = set(sids)
    sid_map = {}
    chrom_pos_dict = {}
    for chrom_i in range(1,23):
        chrom_pos_dict[chrom_i] = [] 
    for chrom_i in range(1,23):
        fn = '/project/PCMA/faststorage/1_DATA/1k_genomes/ALL_1000G_phase1integrated_v3_chr%d_impute.legend.gz'%chrom_i
        with gzip.open(fn) as f:
            f.next()
            for line in f:
                l = line.split()
                sid = l[0]
                if sid in sids:
                    pos = int(l[1])
                    sid_map[l[0]]={'pos':pos, 'chrom':chrom_i}
                    chrom_pos_dict[chrom_i].append(pos)
    return {'sid_map':sid_map, 'chrom_pos_dict':chrom_pos_dict}


def plot_manhattan(result_file,fig_filename='/project/PCMA/faststorage/2_RESULTS/figures/manhattan_all.png'):
    """
    Generates a Manhattan plot for the PCMA results...
    """
    res = pandas.read_csv(result_file)
    sids = list(res.SID)
    print 'Getting SNP positions from 1K Genomes data'
    d = get_sid_pos_map(sids)
    sid_map = d['sid_map']
    chrom_pos_dict = d['chrom_pos_dict']
    print 'Calculating X-axis offsets'
    chrom_offset_dict = {}
    x_tick_pos = []
    x_tick_lab = []
    x_offset = 0
    for chrom_i in range(1,23):
        chrom_offset_dict[chrom_i]=x_offset
        old_x_offset = x_offset
        x_offset += max(chrom_pos_dict[chrom_i])
        x_tick_pos.append((old_x_offset+x_offset)/2.0)
        x_tick_lab.append(str(chrom_i))
    
    print 'Calculating X-axis positions'
    ps = sp.array(res.pvCHI2)
#         ps = sp.array(res.pvCPC)
    x_positions=sp.empty(len(ps))
    chromosomes=sp.empty(len(ps))
    for i, sid in enumerate(sids):
        if sid=='SID': 
            continue
        chrom_i=sid_map[sid]['chrom']
        pos=sid_map[sid]['pos']
        x_offset = chrom_offset_dict[chrom_i]
        x_positions[i]=x_offset+pos
        chromosomes[i]=chrom_i
    
    neg_log_ps = -sp.log10(ps)
    ps_filter = neg_log_ps>3
    filtered_log_ps = neg_log_ps[ps_filter]   
    filtered_pos = x_positions[ps_filter] 
    filtered_chroms = chromosomes[ps_filter]
    
    color_map = {1:{'x_pos':[],'ps':[]}, 2:{'x_pos':[],'ps':[]},
                 3:{'x_pos':[],'ps':[]}, 4:{'x_pos':[],'ps':[]}}
    for lps,pos,chrom in it.izip(filtered_log_ps,filtered_pos,filtered_chroms):
        if chrom%2==0:
            if lps<7.301029:
                color_map[1]['x_pos'].append(pos)
                color_map[1]['ps'].append(lps)
            else:
                color_map[3]['x_pos'].append(pos)
                color_map[3]['ps'].append(lps)
        else:
            if lps<7.301029:
                color_map[2]['x_pos'].append(pos)
                color_map[2]['ps'].append(lps)
            else:
                color_map[4]['x_pos'].append(pos)
                color_map[4]['ps'].append(lps)
        
    print 'Filtering and plotting'
    with pylab.style.context('fivethirtyeight'):
        pylab.figure(figsize=(14,5))
        pylab.plot(color_map[1]['x_pos'],color_map[1]['ps'],'.',color='#1199EE',alpha=0.2)
        pylab.plot(color_map[2]['x_pos'],color_map[2]['ps'],'.',color='#11BB00',alpha=0.2)
        pylab.plot(color_map[3]['x_pos'],color_map[3]['ps'],'.',color='#AA99EE',alpha=0.7)
        pylab.plot(color_map[4]['x_pos'],color_map[4]['ps'],'.',color='#AABB00',alpha=0.7)
        pylab.ylabel('-log(P-value)')
        pylab.xlabel('Chromosomes')
        pylab.xticks(x_tick_pos,x_tick_lab)
        pylab.tight_layout()
        pylab.savefig(fig_filename)
    pylab.clf()
