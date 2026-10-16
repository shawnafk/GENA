import numpy as np

def neff(func):
    """Return a negated version of an eff function."""
    def wrapper(*args, **kwargs):
        return ~func(*args, **kwargs)
    wrapper.__name__ = f"n_{func.__name__}"
    wrapper.__doc__ = f"Negated {func.__name__}"
    return wrapper

#filter level 1
def triggered_eff(df): 
    sum = np.array((df['TX_1s'] + df['TX_2s'] + df['TY_1s'] + df['TY_2s'] + df['TX_1e'] + df['TX_2e'] + df['TY_1e'] + df['TY_2e']))
    return (sum!=0)

def start_eff(df): 
    multi = np.array((df['TX_1s'] * df['TX_2s'] * df['TY_1s'] * df['TY_2s']))
    return (multi!=0)

def start_only_eff(df): 
    multi = np.array((df['TX_1s'] * df['TX_2s'] * df['TY_1s'] * df['TY_2s'])) 
    sum = np.array((df['TX_1e'] + df['TX_2e'] + df['TY_1e'] + df['TY_2e']))
    return (multi!=0) & (sum==0)

def start_stop_eff(df): 
    multi = np.array((df['TX_1s'] * df['TX_2s'] * df['TY_1s'] * df['TY_2s'] * df['TX_1e'] * df['TX_2e'] * df['TY_1e'] * df['TY_2e']))
    return (multi!=0)

#-------------------------------------------------------------------------------------

#second condition
def length_eff(df,params, prob = 'a', dir = 'X', edge ='1', channel = 's',re_err=0.1):
    if prob == 'a':
        vmax = params[f'A_{dir}_SUM']*(1+re_err)
    elif prob == 'b':
        vmax = params[f'B_{dir}_SUM']*(1+re_err)
    return df[f'T{dir}_{edge}{channel}'] < vmax

#std 1
def dir_length_eff(df,params,prob='a',dir='X'):
    return length_eff(df, params, prob,dir,'1','s') & length_eff(df,params, prob,dir,'2','s')


def start_length_eff(df, params, prob = 'a'):
    return dir_length_eff(df,params,prob,'X') & dir_length_eff(df,params,prob,'Y')


#third condition, or not
def sum_eff(df,params,prob = 'a',dir='X',re_err=0.1):
    if prob == 'a':
        xu = params[f'A_{dir}_SUM'] * (1+re_err)
        xl = params[f'A_{dir}_SUM'] * (1-re_err)
    elif prob == 'b':
        xu = params[f'B_{dir}_SUM'] * (1+re_err)
        xl = params[f'B_{dir}_SUM'] * (1-re_err)
    sumx =  np.array((df[f'T{dir}_1s'] + df[f'T{dir}_2s']))
    return (sumx<xu) & (sumx>xl)

#third condition, or not
def diff_eff(df,params,prob = 'a',dir='X',channel='s',re_err=0.1):
    if prob == 'a':
        xu = params[f'A_{dir}R'] * (1+re_err)
        xl = params[f'A_{dir}L'] * (1-re_err)
    elif prob == 'b':
        xu = params[f'B_{dir}R'] * (1+re_err)
        xl = params[f'B_{dir}L'] * (1-re_err)
    diff_1_2 =  np.array((df[f'T{dir}_1{channel}'] - df[f'T{dir}_2{channel}']))
    return (diff_1_2<xu) & (diff_1_2>xl)


#std 2
def start_sum_eff(df,params,prob = 'a'):
    return  sum_eff(df,params,prob,dir='X') & sum_eff(df,params,prob,dir='Y')

#std 3
def start_diff_eff(df,params,prob = 'a'):
    return  diff_eff(df,params,prob,dir='X') & diff_eff(df,params,prob,dir='Y')

#std 4
def end_diff_eff(df,params,prob = 'a'):
    return  diff_eff(df,params,prob,dir='X',channel='e') & diff_eff(df,params,prob,dir='Y',channel='e')

#for those sum_start_x1_x2 in [55, 100], while sum_start_y1_y2 is good, [35, 55], IF
#tof is exists (and trustable), derive x2 time
#anode time still work, time shift happend, use anode_tof to repair x2
#other situation, just use original x2

#for level filtering
#! convert to difference 1-2


def select_circ(df,center=[0,0],radius=10,layer=1):
    if layer == 1:
        x,y = (df['X1'], df['Y1'])
    elif layer == 2:
        x,y = (df['X2'], df['Y2'])
    elif layer == 3:
        x,y = (df['PX'], df['PY'])
    r = np.sqrt((x-center[0])**2 + (y-center[1])**2)
    return (r<radius) 

def select_rect(df,center=[0,0],width=10,height=10,layer=1):
    if layer == 1:
        x,y = (df['X1'], df['Y1'])
    elif layer == 2:
        x,y = (df['X2'], df['Y2'])
    elif layer == 3:
        x,y = (df['PX'], df['PY'])
    w = width / 2
    h = height / 2
    cond1 = (x > center[0] - w) & (x < center[0] + w)
    cond2 = (y > center[1] - h) & (y < center[1] + h)
    return cond1 & cond2


from . import constants as const
#*flag denotes the layer
def ns_2_xy(df,flag,prob,snid=1,params=None):
    if params is None:
        params = const.init_param(snid)
    if prob == 'a':
        KX,BX = params['A_KX'],params['A_BX']
        KY,BY = params['A_KY'],params['A_BY']
    elif prob == 'b':
        KX,BX = params['B_KX'],params['B_BX']
        KY,BY = params['B_KY'],params['B_BY']
    else:
        raise ValueError("prob must be 'a' or 'b'")
    if flag == 0:
        #! use 1 - 2, the sign of x/y coordinates are matched with satellite!!! final
        x = KX * (df['TX_1s'] - df['TX_2s']) + BX
        y = KY * (df['TY_1s'] - df['TY_2s']) + BY
    elif flag == 1:
        x = KX * (df['TX_1e'] - df['TX_2e']) + BX
        y = KY * (df['TY_1e'] - df['TY_2e']) + BY
    else:
        raise ValueError("flag must be 0: carbon file layer or 1 mcp layer")
    return x,y



#-------------------------- for mismatch data
def sum_selected_eff(df,v1,v2,params,prob = 'a',dir='X',re_err=0.1):
    if prob == 'a':
        xu = params[f'A_{dir}_SUM'] * (1+re_err)
        xl = params[f'A_{dir}_SUM'] * (1-re_err)
    elif prob == 'b':
        xu = params[f'B_{dir}_SUM'] * (1+re_err)
        xl = params[f'B_{dir}_SUM'] * (1-re_err)
    sumx =  np.array((df[f'T{dir}_{v1}'] + df[f'T{dir}_{v2}']))
    return (sumx<xu) & (sumx>xl)


def channel_ok(df,dir,prob,params):
    return dir_length_eff(df,params,dir=dir,prob=prob) & ...

from functools import partial
def golden_std(params,prob = 'a',dir='all'):
    if dir == 'all':
        return [
        partial(start_length_eff, params=params, prob=prob),
        partial(start_sum_eff, params=params, prob=prob),
        partial(start_diff_eff, params=params, prob=prob),
        partial(end_diff_eff, params=params, prob=prob),
        partial(tofxy_eff)
        ]
    elif dir == 'X' or dir == 'Y':
        return [
        partial(length_eff, params=params, prob=prob,dir=dir,edge='1',channel='s'),
        partial(length_eff, params=params, prob=prob,dir=dir,edge='2',channel='s'),
        partial(sum_eff, params=params, prob=prob,dir=dir),
        partial(diff_eff, params=params, prob=prob,dir=dir),
        partial(diff_eff, params=params, prob=prob,dir=dir,channel='e'),
        ]

def soft_std(params,prob = 'a',dir='all'):
    if dir == 'all':
        return [
        partial(start_diff_eff, params=params, prob=prob),
        partial(end_diff_eff, params=params, prob=prob),
        ]
    elif dir == 'X' or dir == 'Y':
        return [
        partial(length_eff, params=params, prob=prob,dir=dir,edge='1',channel='s'),
        partial(length_eff, params=params, prob=prob,dir=dir,edge='2',channel='s'),
        partial(sum_eff, params=params, prob=prob,dir=dir),
        partial(diff_eff, params=params, prob=prob,dir=dir),
        partial(diff_eff, params=params, prob=prob,dir=dir,channel='e'),
        ]


#start_stop_eff(df) & match_eff(df) & start_ns_eff(df)
from functools import reduce
#for i in conds[1:]:
#    idx = idx & i
#what if flt need args
def filter_df(df,flts,op='and'):
    if op == 'and':
        combined_mask = reduce(lambda x, y: x & y, [flt(df) for flt in flts])
    elif op == 'or':
        combined_mask = reduce(lambda x, y: x | y, [flt(df) for flt in flts])
    return combined_mask
