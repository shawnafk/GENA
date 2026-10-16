#0224, deleted one

def gen_name(prefix,rawname,prob):
    Para_file = prefix + rawname + '_parameter.csv'
    L1A_file = prefix + rawname + '_L1A_' + prob + '.csv'
    Para_tab = prefix + rawname + '_ptab_' + prob + '.csv'
    return (L1A_file,Para_file,Para_tab)

prefix = r"/Users/shawn/Desktop/genaL1/0224_0301_PartialL1/"
rawname = r"CE7ZJ-ZP-PLD-CMENAVCID-LVDS-20250304T162102-20250304T201843.raw"
#load libs for generating parameter tab
from zhengjs import df as gdf
import pandas as pd
import numpy as np

def generates_para_tab(prefix,rawname,prob):
    fnames = gen_name(prefix,rawname,prob)
    df = gdf.read_csv(fnames[1])
    df['real_BX1'] = df['BX1'] * 2 + df['BG']
    df['real_BX2'] = df['BX2'] * 2 + df['BG']
    df['real_BY1'] = df['BY1'] * 2 + df['BG']
    df['real_BY2'] = df['BY2'] * 2 + df['BG']
    #df['rounded_HV_MCP1'] = df['HV_MCP1'].apply(round_down_to_interval, args=(voltage_range,))
    #df['rounded_HV_MCP2'] = df['HV_MCP2'].apply(round_down_to_interval, args=(voltage_range,))
    if prob == 'a':
        plist = ['MCP1', 'AX', 'AY']
        #plist = ['MCP1', 'AX', 'AY','rounded_HV_MCP1']
        df['combination'] = df[plist].apply(tuple, axis=1)
    if prob == 'b':
        plist = ['MCP2','real_BX1', 'real_BY1', 'real_BY2']
        #plist = ['MCP2','real_BX1', 'real_BX2', 'real_BY1', 'real_BY2','rounded_HV_MCP2']
        df['combination'] = df[plist].apply(tuple, axis=1)

    # 计算相邻行的时间差（以分钟为单位）
    time_differences = df['Date'].diff().dt.total_seconds() / 60

    # 判断 combination 是否变化或者时间差是否大于 30 分钟
    changes = (df['combination'].ne(df['combination'].shift())) | (time_differences > 30)

    # 对变化情况进行累积求和，标记不同的组
    combination_changes = changes.cumsum()
    # 统计每个组合的长度
    combination_lengths = df.groupby(combination_changes).size()
    # 筛选出长度大于等于 300 的组合
    valid_combinations = combination_lengths[combination_lengths >= 300]
    # 生成测试 ID
    test_id = 1
    test_results = []
    if prob == 'a':
        for comb_num, length in valid_combinations.items():
            start_index = df[combination_changes == comb_num].index.min()
            end_index = df[combination_changes == comb_num].index.max()
            combination = df.loc[start_index, 'combination']
            test_results.append({
                'test_id': test_id,
                'MCP': df['MCP1'][start_index],
                'AX': df['AX'][start_index],
                'AY': df['AY'][start_index],
                'start_date': df['Date'][start_index].strftime('%Y-%m-%d-%H:%M:%S'),
                'end_date': df['Date'][end_index].strftime('%Y-%m-%d-%H:%M:%S'),
                'start_ts': df['Timestamp'][start_index],
                'end_ts': df['Timestamp'][end_index],
                'length': length
            })
            test_id += 1
    if prob == 'b':
        for comb_num, length in valid_combinations.items():
            start_index = df[combination_changes == comb_num].index.min()
            end_index = df[combination_changes == comb_num].index.max()
            combination = df.loc[start_index, 'combination']
            test_results.append({
                'test_id': test_id,
                'MCP': df['MCP2'][start_index],
                'BX1': df['real_BX1'][start_index],
                'BX2': df['real_BX2'][start_index],
                'BY1': df['real_BY1'][start_index],
                'BY2': df['real_BY2'][start_index],
                'start_date': df['Date'][start_index].strftime('%Y-%m-%d-%H:%M:%S'),
                'end_date': df['Date'][end_index].strftime('%Y-%m-%d-%H:%M:%S'),
                'start_ts': df['Timestamp'][start_index],
                'end_ts': df['Timestamp'][end_index],
                'length': length
            })
            test_id += 1

    # 将结果转换为 DataFrame
    result_df = pd.DataFrame(test_results)

    # 保存结果到新的 CSV 文件
    result_df.to_csv(fnames[2], index=False)

generates_para_tab(prefix,rawname,'a')
generates_para_tab(prefix,rawname,'b')

#gen channel distribution
import matplotlib
matplotlib.use('Agg')
from zhengjs import quicklook as ql
from zhengjs import df as gdf
from zhengjs import filters as flt
from zhengjs import physical as phy
phy.init_param(1)
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

def show_paras2(ax,r,flag,gap=10,mk='',msize=0):
    #rdate = pd.to_datetime(r.iloc[:, 0],format="%Y-%m-%d %H:%M:%S.%f")
    #ddate = pd.to_datetime(d.iloc[:, 0],format="%Y-%m-%d %H:%M:%S.%f")
    rdate = pd.to_datetime(r.iloc[:, 0])
    if flag == 'a':
        ax[0].plot(rdate[::gap], abs(r['HV_MCP1'][::gap]),marker=mk,markersize=msize)
        ax[1].plot(rdate[::gap], r['TA'][::gap],marker=mk)
    elif flag == 'b':
        ax[0].plot(rdate[::gap], abs(r['HV_MCP2'][::gap]),marker=mk,markersize=msize)
        ax[1].plot(rdate[::gap], r['TB'][::gap],marker=mk,markersize=msize)
    else:
        print('name not correct')
        exit(1)
    ax[0].set_ylabel('Volt')
    ax[0].set_ylim(2600,3000)
    
    #temp
    ax[1].set_ylabel('Temp')
    ax[1].set_ylim(-20,20)
    minor_locator = ticker.MultipleLocator(1)
    ax[1].yaxis.set_minor_locator(minor_locator)
    minor_locator = ticker.MultipleLocator(25)
    ax[0].yaxis.set_minor_locator(minor_locator)
    for i in range(2):
        ax[i].grid(which='minor', linestyle='--', alpha=0.5)
        ax[i].grid(which='major', linestyle='-', alpha=0.8)
    return 0

#A
def gen_plots(prefix,rawname,prob,folder='./'):
    if prob == 'a':
        xlim = phy.A_NS_XLIM
        ylim = phy.A_NS_YLIM
    elif prob == 'b':
        xlim = phy.B_NS_XLIM
        ylim = phy.B_NS_YLIM
    fnames = gen_name(prefix,rawname,prob)
    ptab = pd.read_csv(fnames[2])
    if(len(ptab)) == 0:
        return 1
    pdf = gdf.read_csv(fnames[1])
    df = gdf.read_csv(fnames[0])
    for row in ptab.itertuples():
        subdf = df.split_time(row.start_date,row.end_date)
        if(len(subdf) == 0):
            print("this one is empty")
            continue
        subpdf = pdf.split_time(row.start_date,row.end_date)
        df4 = flt.start_eff(subdf)
        df8 = flt.start_stop_eff(df4)
        ns = flt.start_ns_eff(df8,xlim,ylim)
        match = flt.match_eff(ns)
        level2_df = flt.l1_2_l2(match,prob,1)
        onboard = flt.onboard_eff(level2_df)
        if prob =='a':
            table_data = [
                ['MCP', row.MCP],
                ['AX', row.AX],
                ['AY', row.AY],
                ['StartDate', row.start_date],
                ['EndDate', row.end_date],
            ]
        if prob =='b':
            table_data = [
                ['MCP', row.MCP],
                ['BX1', row.BX1],
                ['BX2', row.BX2],
                ['BY1', row.BY1],
                ['BY2', row.BY2],
                ['StartDate', row.start_date],
                ['EndDate', row.end_date],
            ]

        height_ratios = [1] * 2 + [3] * 8 + [4] * 2
        fig, ax = plt.subplots(12, 1, sharex='col', sharey='row', figsize=(15, 80), gridspec_kw={'height_ratios': height_ratios})
        table = ax[0].table(cellText=table_data, loc='top')
        show_paras2([ax[0],ax[1]],subpdf,'a')
        ql.show_hists_8(ax[2:10],subdf)
        ql.show_counts(ax[10:12],subdf,df4,df8,ns,match,onboard)
        ql.format_date(fig,ax[-1],subdf)


        if prob == 'a':
            fig.savefig(folder + prob
                    + '_' + str(row.MCP) + '_' + str(row.AX) + '_' + str(row.AY)
                    + '_' + pd.to_datetime(row.start_date).strftime('%Y_%m_%d_%H_%M_%S') 
                    + '_to_' + pd.to_datetime(row.end_date).strftime('%Y_%m_%d_%H_%M_%S') 
                    +'.jpg',
                    dpi=300, bbox_inches='tight')
        if prob == 'b':
            fig.savefig(folder + prob 
                    + '_' + str(row.MCP) + '_' + str(row.BX1) + '_' + str(row.BX2) 
                    + '_' + str(row.BY1) + '_' + str(row.BY2)
                    + '_' + pd.to_datetime(row.start_date).strftime('%Y_%m_%d_%H_%M_%S') 
                    + '_to_' + pd.to_datetime(row.end_date).strftime('%Y_%m_%d_%H_%M_%S') 
                    +'.jpg',
                    dpi=300, bbox_inches='tight')
            
#gen_plots(prefix,rawname,'a')
gen_plots(prefix,rawname,'b')
