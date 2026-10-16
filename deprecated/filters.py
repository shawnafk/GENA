import numpy as np
from .. import physical as phy
from ...merge import gio as gdf

#select data with some conditions
#return len and sub dataframe
def select_data(df,cond):
    return df[cond]


#bad condition

#condition
def count_eff(df): 
    sum = np.array((df['TX_1s'] + df['TX_2s'] + df['TY_1s (ns)'] + df['TY_2s'] + df['TX_1e'] + df['TX_2e'] + df['TY_1e'] + df['TY_2e']))
    cond = (sum!=0)
    sub_df = select_data(df,cond)
    return sub_df

#all zero data
def all_zero(df):
    sum = np.array((df['TX_1s'] + df['TX_2s'] + df['TY_1s'] + df['TY_2s'] + df['TX_1e'] + df['TX_2e'] + df['TY_1e'] + df['TY_2e']))
    cond = (sum==0)
    sub_df = select_data(df,cond)
    return sub_df

#4 start effective
def start_eff(df): 
    multi = np.array((df['TX_1s'] * df['TX_2s'] * df['TY_1s'] * df['TY_2s']))
    cond = (multi!=0)
    sub_df = select_data(df,cond)
    return sub_df

#one eff
def start_x1_eff(df):
    return select_data(df,df['TX_1s']!=0)

def start_x2_eff(df):
    return select_data(df,df['TX_2s']!=0)

def start_y1_eff(df):
    return select_data(df,df['TY_1s']!=0)

def start_y2_eff(df):
    return select_data(df,df['TY_2s']!=0)

def x1_eff(df):
    return select_data(df,df['TX_1e']!=0)

def x2_eff(df):
    return select_data(df,df['TX_2e']!=0)

def y1_eff(df):
    return select_data(df,df['TY_1e']!=0)

def y2_eff(df):
    return select_data(df,df['TY_2e']!=0)



#2 iff, pair
def start_x_iff(df):
    cond = (df['TX_1s']==0) & (df['TX_2s']==0)
    return select_data(df,cond)

def start_y_iff(df):
    cond = (df['TY_1s']==0) & (df['TY_2s']==0)
    return select_data(df,cond)

#3 iff

def start_stop_eff(df): 
    multi = np.array((df['TX_1s'] * df['TX_2s'] * df['TY_1s'] * df['TY_2s'] * df['TX_1e'] * df['TX_2e'] * df['TY_1e'] * df['TY_2e']))
    cond = (multi!=0)
    sub_df = select_data(df,cond)
    return sub_df

def start_nostop_eff(df): 
    multi = np.array((df['TX_1s'] * df['TX_2s'] * df['TY_1s'] * df['TY_2s'])) 
    sum = np.array((df['TX_1e'] + df['TX_2e'] + df['TY_1e'] + df['TY_2e']))
    cond = (multi!=0 and sum==0)
    sub_df = select_data(df,cond)
    return sub_df



def start_ns_eff(df,xlim=50,ylim=40):
    cond = (df['TX_1s']<xlim) & (df['TX_2s']) & (df['TY_1s']) & (df['TY_2s']<ylim)
    sub_df = select_data(df,cond)
    return sub_df

def start_x1s_iff(df,xlim):
    cond = df['TX_1s'] > xlim
    sub_df = select_data(df,cond)
    return sub_df

def start_x2s_iff(df,xlim):
    cond = df['TX_2s'] > xlim
    sub_df = select_data(df,cond)
    return sub_df

def start_y1s_iff(df,ylim):
    cond = df['TY_1s'] > ylim
    sub_df = select_data(df,cond)
    return sub_df

def start_y2s_iff(df,ylim):
    cond = df['TY_2s'] > ylim
    sub_df = select_data(df,cond)
    return sub_df

def match_eff(df,xl=45,xu=55,yl=35,yu=45):
    sumx =  np.array((df['TX_1s'] + df['TX_2s']))
    sumy =  np.array((df['TY_1s'] + df['TY_2s']))
    cond = (sumx<xu) & (sumx>xl) & (sumy<yu) & (sumy>yl)
    sub_df = select_data(df,cond)
    return sub_df

def match_iff(df,xl=45,xu=55,yl=35,yu=45):
    sumx =  np.array((df['TX_1s'] + df['TX_2s']))
    sumy =  np.array((df['TY_1s'] + df['TY_2s']))
    cond = (sumx>xu) | (sumx<xl) | (sumy>yu) | (sumy<yl)
    sub_df = select_data(df,cond)
    return sub_df

    
#level 2
def l1_2_l2(df,flag,snid):
    phy.init_param(snid)
    if flag == 'a':
        xc = phy.A_X_CENTER_NS
        yc = phy.A_Y_CENTER_NS
        xcoef = phy.A_X_COF
        ycoef = phy.A_Y_COF
    if flag == 'b':
        xc = phy.B_X_CENTER_NS
        yc = phy.B_Y_CENTER_NS
        xcoef = phy.B_X_COF
        ycoef = phy.B_Y_COF
    x1 = (df['TX_2s'] - df['TX_1s'])
    y1 = (df['TY_2s'] - df['TY_1s'])
    x2 = (df['TX_2e'] - df['TX_1e'])
    y2 = (df['TY_2e'] - df['TY_1e'])

    x1 =  (x1 - xc) * xcoef
    y1 =  (y1 - yc) * ycoef
    x2 =  (x2 - xc) * xcoef
    y2 =  (y2 - yc) * ycoef
    px = np.arctan((x2-x1)/50)*180/np.pi
    py = np.arctan((y2-y1)/50)*180/np.pi
    return gdf.ENA_DataFrame({'Date': df['Date'],'Timestamp': df['Timestamp'], 'X1': x1, 'Y1': y1, 'X2': x2, 'Y2': y2, 'Ang_X': px, 'Ang_Y': py})

def onboard_eff(df,center=[0,0],c_size = [90,30],m_size=[150,60]):
    cxl = center[0] - c_size[0]/2
    cxr = center[0] + c_size[0]/2
    cyd = center[1] - c_size[1]/2
    cyu = center[1] + c_size[1]/2
    
    mxl = center[0] - m_size[0]/2
    mxr = center[0] + m_size[0]/2
    myd = center[1] - m_size[1]/2
    myu = center[1] + m_size[1]/2
    
    cond1 = (df['X1']>cxl) & (df['Y1']>cyd) & (df['X1']<cxr) & (df['Y1']<cyu)
    cond2 = (df['X2']>mxl) & (df['Y2']>myd) & (df['X2']<mxr) & (df['Y2']<myu)
    sub_df = select_data(df,cond1&cond2)
    return sub_df
    
def angle_eff(df,Ax,Ay):
    x1,y1 = (df['X1'], df['Y1'])
    x2,y2 = (df['X2'], df['Y2'])
    px = np.arctan((x1-x2)/50)*180/np.pi
    py = np.arctan((y1-y2)/50)*180/np.pi
    cond1 = (px>-Ax) & (px<Ax)
    cond2 = (py>-Ay) & (py<Ay)
    sub_df = select_data(df,cond1&cond2)
    return sub_df

def energy_eff(df,Emin,Emax):
    sub_df = select_data(df,(df['Energy']>Emin) & (df['Energy']<Emax))
    return sub_df
