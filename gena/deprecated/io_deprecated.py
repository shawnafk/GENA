''' Performance
pqt io: gzip, default compress: 1G l1a 500M l1b, 73M pdf
====================SUMMARY====================
Num packets processed 2497030
Load RAW to memory costs: 0.18 s
Searching frame index costs: 8.11 s
Parse and dump parameters costs: 1.29 s
Parse and dump df0 costs: 40.48 s
Total time costs: 50.11 s

hdf v1 io: complib='blosc',complevel=5. 3.4G in total
====================SUMMARY====================
Num packets processed 2497030
Load RAW to memory costs: 0.19 s
Searching frame index costs: 2.63 s
Parse and dump parameters costs: 3.76 s
Parse and dump df0 costs: 203.25 s
Total time costs: 209.86 s

hdf v2 io: gzip, default level (maybe 4) 1.0G
====================SUMMARY====================
Num packets processed 2497030
Load RAW to memory costs: 0.17 s
Searching frame index costs: 2.34 s
Parse and dump parameters costs: 3.73 s
Parse and dump df0 costs: 139.62 s
Total time costs: 145.85 s

level 2: 1.2G in total
====================SUMMARY====================
Num packets processed 2497030
Load RAW to memory costs: 0.16 s
Searching frame index costs: 7.87 s
Parse and dump parameters costs: 2.37 s
Parse and dump df0 costs: 84.62 s
Total time costs: 95.02 s

level 1: 1.2G in total
====================SUMMARY====================
Num packets processed 2497030
Load RAW to memory costs: 0.18 s
Searching frame index costs: 2.22 s
Parse and dump parameters costs: 2.40 s
Parse and dump df0 costs: 83.81 s
Total time costs: 88.61 s

level 0 do not compress: 18G
====================SUMMARY====================
Num packets processed 2497030
Load RAW to memory costs: 0.85 s
Searching frame index costs: 3.55 s
Parse and dump parameters costs: 1.36 s
Parse and dump df0 costs: 54.85 s
Total time costs: 60.62 s
'''
#------------------------------------------------------------------ hdf v1: complib='blosc',complevel=5
def init_hdf(filename):
    pdf = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.pdf_type.items()})
    l1 = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1_type.items()})
    l2 = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df2_type.items()})
    with pd.HDFStore(filename, mode='w',complib='blosc',complevel=5,index=True) as store:
        store.put('/P', pdf, format='table', data_columns=['Timestamp'],index=True)
        for prob in ['/A', '/B']:
            store.put(prob + '/L1', l1, format='table', data_columns=['Timestamp'],index=True)
            store.put(prob + '/L2', l2, format='table', data_columns=['Timestamp'],index=True)


def save_hdf(filename,df,group_name):
    #file_exists = os.path.exists(filename)
    #df = dict_2_df(df_dict,df_type)
    df.to_hdf(filename, key=group_name, mode='a',append=True,complib='blosc',complevel=5, format='table',data_columns=['Timestamp'] )
    
#partial read example
#low = pd.to_datetime('20250312 00:00:00').timestamp()
#high = pd.to_datetime('20250312 12:00:00').timestamp()
#
#df = pd.read_hdf(
#    'gena_lv5.h5',
#    key='P',
#    where=f'Timestamp > {low} & Timestamp < {high}'
#
#eg: pdf = gio.from_hdf('gena_lv5.h5',key='P',where=f'Timestamp > {low} & Timestamp < {high}')
def from_hdf(filename,*args,**kwargs):
    df = pd.read_hdf(filename,*args,**kwargs).drop_duplicates().sort_values(by='Timestamp') 
    edf = ENA_DataFrame(df)
    return edf


#---------------------------- self organized hdf file --------------------------------------
import h5py
from typing import Callable, List

def gena_h5_init(file_path: str):
    if not os.path.exists(file_path):
        #w- will write on the file
        with h5py.File(file_path, 'w-') as f:
            f.create_group('P')
            f.create_group('L1/A/')
            f.create_group('L1/B/')
            f.create_group('L2/A/')
            f.create_group('L2/B/')
    
def gena_h5_write(file_path, df: pd.DataFrame, group_name: str, chunk_id: int) -> None:
    """
    添加一个 DataFrame chunk 到存储
    参数:
        df: 要存储的 DataFrame
    """
    if 'Timestamp' not in df.columns:
        raise ValueError("DataFrame must contain 'Timestamp' columns")
    #chunk id is the id of execution of add_chunk
    chunk_id_str = str(chunk_id).zfill(4)
    with h5py.File(file_path, 'a') as f:
        # create chunk subgroups
        chunk_group = f[group_name].create_group(f'c_{chunk_id_str}')
        
        #Create dataset for every column
        for col in df.columns:
            data = df[col].values
            dtype = df[col].dtype
            chunk_group.create_dataset(
                col,
                data=data,
                compression='gzip',
                #compression='blosc',
                #compression='lzo',
                compression_opts = 1,
                chunks=True,
                dtype=dtype
            )



#default is load all rows, idx_list = None; :
#from memory_profiler import profile
#@profile
def gena_h5_load_cols(file_path, group_name: str, cols: List[str], rows: List[List[bool]] = None) -> pd.DataFrame:
    """
    从指定的多个 key 中提取所有 chunks 的数据，并组合成一个完整的 DataFrame
    
    参数:
        cols: column names of the loading HDF5 data
        rows: Line index of the loading HDF5 data, lists of bools lists,
                first lists denots different chunks;
                seconds denotes rows; its length should match data length in the chunk.
    
    返回:
        pd.DataFrame: 组合后的完整 DataFrame
    """
    # 存储所有列的数据
    #all_data = {}
    # 记录每个列的数据类型（确保一致性）
    #dtypes = {}
    
    store_dict = {col: [] for col in cols}

    with h5py.File(file_path, 'r') as f:
        chunks = f[group_name]
        num_chunks = len(chunks)
        if rows == None:
        #performance: change loop order doesn't affect efficience
        #directly read to dict is 14 sec
        #change to dataframe cost another 15 sec
            for k in cols:
                dat = []
                for c in range(num_chunks):
                    chunk_id_str = 'c_'+str(c).zfill(4)
                    dat.append(chunks[chunk_id_str][k][:])    
                store_dict[k] = np.concatenate(dat)
        else:
            #non_empty_indices = non_empty_indices = [(c, idx) for c, idx in zip(range(num_chunks), rows) if sum(idx)!= 0] 
            #non_empty_rows = [idx for _, idx in non_empty_indices]
            #non_empty_chunks = [c for c, _ in non_empty_indices]
            for k in cols:
                dat = []
                #for c,idx in zip(non_empty_chunks,non_empty_rows):
                for c,idx in zip(range(num_chunks),rows):
                    #tstart = time.time()
                    #if(sum(idx)==0):
                        #print('skiped')
                        #continue
                    chunk_id_str = 'c_'+str(c).zfill(4)
                    dat.append(chunks[chunk_id_str][k][idx])    
                    #tend = time.time()
                    #print('time load partial:', tend-tstart)
                store_dict[k] = np.concatenate(dat)   
    #df = pd.DataFrame(store_dict)
    return store_dict

#all columns; with row id or not
def gena_h5_load(file_path, group_name: str, cols: List[str]=None, rows: List[List]=None):
    if cols == None:
        if group_name == 'P':
            all_cols = [col for col,_ in st.pdf_type.items()]
    
        if group_name == 'L1/A' or group_name == 'L1/B':
            all_cols = [col for col,_ in st.df1_type.items()]
    
        if group_name == 'L2/A' or group_name == 'L2/B':
            all_cols = [col for col,_ in st.df2_type.items()]
        return gena_h5_load_cols(file_path,group_name, all_cols, rows)
    else:
        return gena_h5_load_cols(file_path,group_name, cols, rows)


def gena_h5_query(file_path: str, group_name: str, cols: List[str], cond: Callable) -> List:
    """
    query by cond function which takes a column array
    return :
        Lists of bool [rows] for every chunk
    """
    idx_list = []
    with h5py.File(file_path, 'r') as f:
        # find all chunks within the col
        chunks = f[group_name]
        num_chunks = len(chunks)
        #only load chunk data for cond fun
        chunk_id = []
        row_flag_in_chunk = []
        for c in range(num_chunks):
            chunk_id_str = 'c_'+str(c).zfill(4)
            dat = []
            for k in cols:
                dat.append(chunks[chunk_id_str][k][:])
            select = cond(dat)
            if np.any(select):
                chunk_id.append(c)
                row_flag_in_chunk.append(select)
    return chunk_id,row_flag_in_chunk

def query_time(t1,t2):
    #cond takes args as form ['Timestamp']
    def cond(Timestamp: List):
        return (Timestamp[0] > t1) & (Timestamp[0] < t2)
    return cond




'''Performance

In [5]: %time df = gio.gena_h5_load_from('gena_8.h5','L1/A',low,high)
CPU times: user 1min 13s, sys: 11.3 s, total: 1min 24s
Wall time: 1min 27s
with partial load. now reduce to 19 for 3
memory get reduce too

But only In [4]: %time df = gio.gena_h5_load('gena_8.h5','L1/A')
CPU times: user 14.8 s, sys: 11.8 s, total: 26.7 s
Wall time: 35.3 s

to load all
'''

#------------------------------------------------------------------ pqt 
from pathlib import Path
def init_pqt2(filepath):
    """
    在指定路径下创建多层级文件夹结构
    
    参数:
    base_path (str): 基础路径，如 '/data/project'
    """
    # 定义需要创建的文件夹结构
    subfolders = {
        'META': [],
        'P': [],
        'L1': ['A', 'B'],
        'L2': ['A', 'B']
    }
    
    # 使用Pathlib处理路径，确保兼容性
    base = Path(filepath)
    
    try:
        base.mkdir(parents=True, exist_ok=True)
        
        for parent, children in subfolders.items():
            parent_path = base / parent
            parent_path.mkdir(exist_ok=True)
            
            for child in children:
                child_path = parent_path / child
                child_path.mkdir(exist_ok=True)
                
        print(f"Create '{filepath}' successed")
        
    except Exception as e:
        print(f"Error when creating folders: {e}")
    
    return 0

#def dict_2_df(dict,df_type):
#    return pd.DataFrame(dict).astype(df_type)


def save_pqt2(df,filepath,subfolder,id):
    df.to_parquet(filepath+subfolder+str(id), engine='fastparquet',compression='gzip')
    return 0