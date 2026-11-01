#part 1 packets define and seek method
import numpy as np
# 定义结构体
#will use array spped up? 
ep1_struct = np.dtype([
    ('HV_N', np.dtype('>i2')),
    ('HV_P', np.dtype('>u2')),
    ('HV_C', np.dtype('>i2')),
    ('HV_G', np.dtype('>i2')),
    ('HV_A', np.dtype('>i2')),
    ('HV_B', np.dtype('>i2')),
    # Monitor ADC bytes are unsigned; int8 turns codes >127 into negative values.
    ('monitors', (np.uint8, 8)),
    ('temps', (np.uint8, 3)),
    ('stat', np.uint8),  
    ('Evt_Num_Tot', np.dtype('>u2')),
    ('Evt_Num_Eff', np.dtype('>u2')),
    ('Evt_Num_Aux', np.dtype('>u2')),
    ('instruction_counts', (np.uint8, 2)),
])

ep2_struct = np.dtype([
    ('version', np.uint8),  
    ('err', np.uint8),      
    ('stat', np.uint8),     
    ('flex', np.uint8, 24),     
    ('args_addr', np.uint8),     
    ('args_value', (np.uint8, 4)),     
])

event_struct = np.dtype([
    ('flag', np.uint8),     
    ('time_s', np.dtype('>u2')),     
    ('time_us', np.dtype('>u2')),     
    ('se_amp', (np.uint8, 3)),
    ('c0', (np.dtype('>u2'), 4)),       
    ('c1', (np.dtype('>u2'), 2)),       
    ('c3', (np.dtype('>u2'), 2)),       
    ('c5', (np.dtype('>u2'), 2)),       
    ('c7', (np.dtype('>u2'), 2)),       
])

packet_struct = np.dtype([
    ('header', (np.uint8, 2)),
    ('seq1', np.dtype('>u2')),
    ('seq2', np.dtype('>u2')),
    ('dataLen', np.dtype('>u2')),
    ('time_s', np.dtype('>u4')),  
    ('time_ms', np.dtype('>u2')), 
    ('ep1', ep1_struct),        
    ('ep2', ep2_struct),        
    ('events', (event_struct, 62)),
    ('verification_sum', np.dtype('>u2')), 
])


pdf_type = {
        'PacketID': np.uint32,
        'BM': np.uint8,
        'BX': np.uint16,
        'Timestamp': np.float64,
    
        'HV_N': np.float32,
        'HV_P': np.float32,
        'HV_C': np.float32,
        'HV_G': np.float32,
        'HV_A': np.float32,
        'HV_B': np.float32,
        
        'A_N12' : np.float32,
        'A_N5'  : np.float32,
        'A_3_3' : np.float32,
        'D_3_3' : np.float32,
        'HV_N_I': np.float32,
        'HV_C_I': np.float32,
        #'TDC_V' : np.float32,
        #'D_1_5' : np.float32,
        
        'Temp_A': np.float32,
        'Temp_B': np.float32,
        'Temp_PCB': np.float32,

        'Evt_Num_Tot': np.uint16,
        'Evt_Num_Eff': np.uint16,
        'Evt_Num_Aux': np.uint16,
    
        'MJ2000_q1': np.float32,
        'MJ2000_q2': np.float32,
        'MJ2000_q3': np.float32,
        'MJ2000_Star_X': np.int32,
        'MJ2000_Star_Y': np.int32,
        'MJ2000_Star_Z': np.int32,
        
        'Ant_Ang': np.float32,
        'Ant_Azi': np.float32,
        'Star_L' : np.float32,
        'Star_NZ': np.float32,
        'Shadow' : np.uint8,
        'PlatformGroup': np.uint8,
        'PosePresent': np.bool_,
        'ArgsAddress': np.uint8,
        'ArgsValue': np.uint32,
        
        
                
        'Th_A_MCP': np.int16,
        'Th_A_X':   np.int16,
        'Th_A_Y':   np.int16,
        'Th_B_MCP': np.int16,       
        'Th_B_X1':  np.int16,     
        'Th_B_X2':  np.int16,   
        'Th_B_Y1':  np.int16,   
        'Th_B_Y2':  np.int16,   
        'Th_B_G':   np.int8,
    }

df0_type ={
        'UUID':          np.dtype('S36'),
        'PID':           np.uint8,
        'BM':            np.uint8,
        'BX':            np.uint16,
        'Timestamp':     np.float64,
        # Internal seconds counter can reach ~1e4 s; float32 loses the 20 us tick.
        'TimeInternal':  np.float64,
        'EvtLevel':      np.uint8,
        'EP':            np.uint16,
        'prob_flag':     np.uint8,
        'ch_tof_0':      np.uint16,
        'ch_tof_1':      np.uint16,
        'ch_tof_2':      np.uint16,
        'ch_tof_3':      np.uint16,
        'ch_1s':         np.uint16,
        'ch_1e':         np.uint16,
        'ch_3s':         np.uint16,
        'ch_3e':         np.uint16,
        'ch_5s':         np.uint16,
        'ch_5e':         np.uint16,
        'ch_7s':         np.uint16,
        'ch_7e':         np.uint16,
        'is_empty_slot': np.bool_,
        'PacketID':      np.uint32,
        'EventID':       np.uint64,
}


df1_type ={
        'UUID':          np.dtype('S36'),
        'PacketID':      np.uint32,
        'EventID':       np.uint64,
        'PID':           np.uint8,
        'BM':            np.uint8,
        'BX':            np.uint16,
        'Timestamp':     np.float64,
        # Keep the same precision as L0; float32 discards 20 us ticks at large values.
        'TimeInternal':  np.float64,
        'EvtLevel':      np.uint8,
        'EP':            np.float32,
        'Tof_0':         np.float32,
        'Tof_1':         np.float32,
        'Tof_2':         np.float32,
        'Tof_3':         np.float32,
        'TX_1s':         np.float32,
        'TX_1e':         np.float32,
        'TX_2s':         np.float32,
        'TX_2e':         np.float32,
        'TY_1s':         np.float32,
        'TY_1e':         np.float32,
        'TY_2s':         np.float32,
        'TY_2e':         np.float32,
}

#1a is more like level 2 (without TOF, SPEC, and ENEGRY);
#no need to keey all channel time
#Entities should not be multiplied unnecessarily
#df1A_type ={
#        'BM':        np.uint8,
#        'BX':        np.uint16,
#        'Timestamp': np.int64,
#        'EvtLevel':  np.uint8,
#        'EP':        np.float32,
#        'X1':        np.float32,
#        'X2':        np.float32,
#        'Y1':        np.float32,
#        'Y2':        np.float32,
#        'Tof_M':       np.float32,
#        'Tof_X':       np.float32,
#        'Tof_Y':       np.float32,
#}

#Below are reduced data, only keep 8 channel recorded data
#markd as 1A
#!also interplate all parameter for further demostration (previous pdf and data are stored separatly, now no need to load them twice)
#also include interplated parameters for -
#- calculate energy/tof
#- statellite attitude and position, incident angle

df1a_type ={
        'UUID':          np.dtype('S36'),
        'PacketID':      np.uint32,
        'EventID':       np.uint64,
        'PID':           np.uint8,
        'BM':            np.uint8,
        'BX':            np.uint16,
        
        'HV_C': np.float32,
        'HV_G': np.float32,
        'MJ2000_q1': np.float32,
        'MJ2000_q2': np.float32,
        'MJ2000_q3': np.float32,
        'MJ2000_Star_X': np.float32,
        'MJ2000_Star_Y': np.float32,
        'MJ2000_Star_Z': np.float32,
        
        'Timestamp':     np.float64,
        # L1A inherits the L0/L1 internal-time precision.
        'TimeInternal':  np.float64,
        'EvtLevel':      np.uint8,
        'EP':            np.float32,
        'Tof_0':         np.float32,
        'Tof_1':         np.float32,
        'Tof_2':         np.float32,
        'Tof_3':         np.float32,
        'TX_1s':         np.float32,
        'TX_1e':         np.float32,
        'TX_2s':         np.float32,
        'TX_2e':         np.float32,
        'TY_1s':         np.float32,
        'TY_1e':         np.float32,
        'TY_2s':         np.float32,
        'TY_2e':         np.float32,
}

'''
Type (数据诊断与修复等级) 定义：

[最高置信度：完美物理数据]
0: Golden Data (完全正常，无任何问题，直接采用)

[中高置信度：标签时序错位 (狸猫换太子)]
1: Swap_X1 (TX_1s与TX_1e错位，修复后完美)
2: Swap_X2 (TX_2s与TX_2e错位，修复后完美)
3: Swap_Y1 (TY_1s与TY_1e错位，修复后完美)
4: Swap_Y2 (TY_2s与TY_2e错位，修复后完美)

[中等置信度：系统性补偿]
10: Rollover_Global (全局溢出翻转，起始正常，所有截止信号加上翻转周期后完美)
20: Time_Shift_Global (全局时间平移，扣除起始平移量后完美)

[低置信度：降维抢救 (好轴救坏轴/一维能谱保留)]
30: Reconstruct_X1 (X1完全丢失/被噪声填充，利用完美的Y轴和X2重构X1)
31: Reconstruct_X2 (利用Y轴和X1重构X2)
40: 1D_Y_Only (X轴彻底损毁无法重构，剥离X轴，仅保留Y轴坐标和整体能量)
41: 1D_X_Only (Y轴彻底损毁无法重构，剥离Y轴，仅保留X轴坐标和整体能量)

[零置信度：彻底报废]
255: Trash (纯噪声或严重污染，无法修复，直接丢弃)
'''
# Level 2 derived event products. GRP was never populated and has been removed.

df2_type ={
        'UUID':          np.dtype('S36'),
        'PacketID':      np.uint32,
        'EventID':       np.uint64,
        'PID':           np.uint8,
        'BM':        np.uint8,
        'BX':        np.uint16,
        'HV_C': np.float32,
        'HV_G': np.float32,
        'MJ2000_q1': np.float32,
        'MJ2000_q2': np.float32,
        'MJ2000_q3': np.float32,

        'MJ2000_Star_X': np.float32,
        'MJ2000_Star_Y': np.float32,
        'MJ2000_Star_Z': np.float32,

        'Timestamp': np.float64,
        # L2 inherits the L1A internal-time precision.
        'TimeInternal':  np.float64,
        'EvtLevel':  np.uint8,
        'EP':        np.float32,
        'X1':        np.float32,
        'X2':        np.float32,
        'Y1':        np.float32,
        'Y2':        np.float32,
        'PX':        np.float32,
        'PY':        np.float32,
        'Tof_M':     np.float32,
        'Tof_X':     np.float32,
        'Tof_Y':     np.float32,
        'Tof':       np.float32,
        'Spec':      np.uint8,
        'Energy':    np.float32,
}
