import os, sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

'''
Check if a given IBDMap/IBDReduce analysis has a genome-wide significant signal
'''
def check_result(ibdreduce_file, significance_threshold=0.05, write_output=True):

    df = pd.read_csv(ibdreduce_file, delimiter='\t', skiprows=3)
    
    num_values = (df['PAdj'] <= significance_threshold).sum()
    if num_values > 0:

        if write_output:

            o = []
            dff = df[df['PAdj'] <= significance_threshold].reset_index(drop=True)
            for x in range(len(dff)):
                minpf = dff.loc[x, 'PAdj']
                adjpf = dff.loc[x, 'PVal']
                chromf = dff.loc[x, 'CHROM']
                posf = dff.loc[x, 'POS']
                o.append({"chrom": chromf, "pos": posf, "pval":adjpf})
            
            odf = pd.DataFrame(o)
            odf.to_csv(ibdreduce_file+'.results', sep='\t', index=None)

        return True
    
    else:

        return False

'''
Build signal "windows" from consecutive genome-wide significant breakpoints 
with one non-significant flanking breakpoint on each side (if they exist)
'''
def build_signal_windows(ibdreduce_file, significance_threshold=0.05):

    df = pd.read_csv(ibdreduce_file, sep='\t', skiprows=3)
    df['CHROM'] = pd.to_numeric(df['CHROM'])
    df = df.sort_values(['CHROM', 'POS'])
    df.reset_index(drop=True, inplace=True)
    df['idx'] = df.index
    data = df.values

    i = 0
    in_group = False
    overall_windows = []
    sig_breakpoints = 0

    while i < len(data):

        current_window = []
        current_lowest_pval = 1.0
        
        while i < len(data) and float(data[i][5]) < significance_threshold:

            if data[i][5] < current_lowest_pval:
                current_lowest_pval = data[i][5]
            chrom = data[i][0]

            try:
                left_flank = data[i-1][1]
                if int(left_flank) > int(data[i][1]):
                    left_flank = data[i][1]
            except:
                left_flank = data[i][1]
            try:
                right_flank = data[i+1][1]
                if int(right_flank) < int(data[i][1]):
                    right_flank = data[i][1]
            except:
                right_flank = data[i][1]

            tup = (i, chrom, left_flank, data[i][1], right_flank)
            current_window.append(tup)
            i += 1
            in_group = True

        if i + 1 < len(data) and float(data[i + 1][5]) < significance_threshold:
            i += 1

        while i < len(data) and float(data[i][5]) < significance_threshold:

            if data[i][5] < current_lowest_pval:
                current_lowest_pval = data[i][5]
            chrom = data[i][0]

            try:
                left_flank = data[i-1][1]
                if int(left_flank) > int(data[i][1]):
                    left_flank = data[i][1]
            except:
                left_flank = data[i][1]
            try:
                right_flank = data[i+1][1]
                if int(right_flank) < int(data[i][1]):
                    right_flank = data[i][1]
            except:
                right_flank = data[i][1]

            tup = (i, chrom, left_flank, data[i][1], right_flank)
            current_window.append(tup)
            i += 1
            in_group = True
        
        if in_group:
            sig_breakpoints += 1
            overall_chr = current_window[0][1]
            start_idx = current_window[0][0]
            start_pos = current_window[0][2]
            end_idx = current_window[-1][0]
            end_pos = current_window[-1][-1]
            window_obj = {"chrom": int(overall_chr), "start": int(start_pos), "end": int(end_pos), "low_pval": current_lowest_pval}
            overall_windows.append(window_obj)

            current_lowest_pval = 1.0
            in_group = False

        i += 1

    df=pd.DataFrame(overall_windows)
    df.to_csv(ibdreduce_file+'.windows', sep='\t', index=None)
    
    return(len(overall_windows))

'''
Generate a Manhattan plot from IBDReduce results file.
'''

def manhattan_plot(ibdreduce_file):

    df = pd.read_csv(ibdreduce_file, sep='\t', skiprows=3)
    df['log_p'] = -np.log10(df['PVal'])
    sig_line = df['PAdjCutoff'].iloc[0]
   
    plt.figure(figsize=(10, 5))
    
    # Convert chromosomes to numeric and sort
    df['CHROM'] = pd.to_numeric(df['CHROM'])
    chromosomes = sorted(df['CHROM'].unique())
    colors = ['#1f77b4', '#2ca02c']
    
    x_pos = 0
    x_ticks = []
    x_labels = []
   
    for idx, chrom in enumerate(chromosomes):
        chrom_data = df[df['CHROM'] == chrom]
        plt.scatter(chrom_data['POS'] + x_pos, 
                    chrom_data['log_p'],
                    c=colors[idx % 2],
                    s=1,
                    alpha=0.7)
        
        x_ticks.append(x_pos + (chrom_data['POS'].max() - chrom_data['POS'].min())/2)
        x_labels.append(str(int(chrom)))
        x_pos += chrom_data['POS'].max()
    
    if sig_line:
        plt.axhline(y=-np.log10(sig_line), color='red', linestyle='--')
    
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p-value)')
    plt.title('Example Manhattan plot')
    plt.xticks(x_ticks, x_labels)
    plt.ylim(0,7)
    plt.grid(False)
    
    plt.tight_layout()
    plt.savefig('manhattan_plot.png', dpi=300, bbox_inches='tight')
    plt.show()
