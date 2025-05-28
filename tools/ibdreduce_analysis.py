import os, sys
import pandas as pd
import matplotlib.pyplot as plt

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