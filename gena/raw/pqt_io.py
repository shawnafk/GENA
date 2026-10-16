import pandas as pd
import os
from . import structure as st
import sys
# -------------------------- Constants --------------------------
PARQUET_ENGINE = "fastparquet"
TIME_UNIT = "s"
COMPRESSION = "gzip"
SUPPORTED_SUFFIX = {
    "pdf": st.pdf_type,
    "l1_a": st.df1_type,
    "l1_b": st.df1_type,
    "l1a_a": st.df1a_type,
    "l1a_b": st.df1a_type,
    "l1b_a": st.df1b_type,
    "l1b_b": st.df1b_type,
    "l2_a": st.df2_type,
    "l2_b": st.df2_type,
    "l2b_a": st.df2_type,
    "l2b_b": st.df2_type,
    "l2s_a": st.df2_type,
    "l2s_b": st.df2_type,
    "l2g_a": st.df2_type,
    "l2g_b": st.df2_type,
}

# -------------------------- Helper Function --------------------------
def _get_empty_df(suffix):
    if suffix not in SUPPORTED_SUFFIX:
        print("Error: unsupported suffix")
        sys.exit(1)
    dtype_map = SUPPORTED_SUFFIX[suffix]
    return pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in dtype_map.items()})

# -------------------------- Initialize Files --------------------------
def init_pqt(filepath):
    try:
        os.makedirs(filepath, exist_ok=True)

        # Initialize pdf
        pdf = _get_empty_df("pdf")
        pdf_path = os.path.join(filepath, "pdf.pqt")
        pdf.to_parquet(pdf_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        # Initialize l1_a / l1_b
        l1_a = _get_empty_df("l1_a")
        l1_a_path = os.path.join(filepath, "l1_a.pqt")
        l1_a.to_parquet(l1_a_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        l1_b = _get_empty_df("l1_b")
        l1_b_path = os.path.join(filepath, "l1_b.pqt")
        l1_b.to_parquet(l1_b_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        # Initialize l1a_a / l1a_b
        l1a_a = _get_empty_df("l1a_a")
        l1a_a_path = os.path.join(filepath, "l1a_a.pqt")
        l1a_a.to_parquet(l1a_a_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        l1a_b = _get_empty_df("l1a_b")
        l1a_b_path = os.path.join(filepath, "l1a_b.pqt")
        l1a_b.to_parquet(l1a_b_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        # Initialize l1b_a / l1b_b
        l1b_a = _get_empty_df("l1b_a")
        l1b_a_path = os.path.join(filepath, "l1b_a.pqt")
        l1b_a.to_parquet(l1b_a_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        l1b_b = _get_empty_df("l1b_b")
        l1b_b_path = os.path.join(filepath, "l1b_b.pqt")
        l1b_b.to_parquet(l1b_b_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        # Initialize l2_a / l2_b
        l2_a = _get_empty_df("l2_a")
        l2_a_path = os.path.join(filepath, "l2_a.pqt")
        l2_a.to_parquet(l2_a_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        l2_b = _get_empty_df("l2_b")
        l2_b_path = os.path.join(filepath, "l2_b.pqt")
        l2_b.to_parquet(l2_b_path, engine=PARQUET_ENGINE, compression=COMPRESSION)


        # Initialize l2_a / l2_b
        l2b_a = _get_empty_df("l2b_a")
        l2b_a_path = os.path.join(filepath, "l2b_a.pqt")
        l2b_a.to_parquet(l2b_a_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        l2b_b = _get_empty_df("l2b_b")
        l2b_b_path = os.path.join(filepath, "l2b_b.pqt")
        l2b_b.to_parquet(l2b_b_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        # Initialize l2_a / l2_b
        l2s_a = _get_empty_df("l2s_a")
        l2s_a_path = os.path.join(filepath, "l2s_a.pqt")
        l2s_a.to_parquet(l2s_a_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        l2s_b = _get_empty_df("l2s_b")
        l2s_b_path = os.path.join(filepath, "l2s_b.pqt")
        l2s_b.to_parquet(l2s_b_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        # Initialize l2_a / l2_b
        l2g_a = _get_empty_df("l2g_a")
        l2g_a_path = os.path.join(filepath, "l2g_a.pqt")
        l2g_a.to_parquet(l2g_a_path, engine=PARQUET_ENGINE, compression=COMPRESSION)

        l2g_b = _get_empty_df("l2g_b")
        l2g_b_path = os.path.join(filepath, "l2g_b.pqt")
        l2g_b.to_parquet(l2g_b_path, engine=PARQUET_ENGINE, compression=COMPRESSION)


    except Exception as e:
        print("init_pqt failed:", str(e))
        sys.exit(1)

# -------------------------- Save (Append Mode) --------------------------
def save_pqt(df, filepath, file):
    try:
        if df.empty:
            print("Warning: empty dataframe, skip writing")
            return

        os.makedirs(filepath, exist_ok=True)
        file_path = os.path.join(filepath, f"{file}.pqt")
        df.to_parquet(
            file_path,
            engine=PARQUET_ENGINE,
            compression=COMPRESSION,
            append=True
        )

    except Exception as e:
        print("save_pqt failed:", str(e))
        sys.exit(1)

# -------------------------- Save Daily Files --------------------------
def save_pqt_daily(df, folder_path, suffix):
    try:
        if df.empty:
            print("Warning: empty dataframe, skip writing")
            return

        if "Timestamp" not in df.columns:
            print("Error: Timestamp column missing")
            sys.exit(1)

        # Convert timestamp to date
        temp_dates = pd.to_datetime(df["Timestamp"], unit=TIME_UNIT).dt.date
        temp_dates_str = temp_dates.apply(lambda x: x.strftime("%Y%m%d"))
        unique_dates = temp_dates.unique()

        empty_df = _get_empty_df(suffix)
        os.makedirs(folder_path, exist_ok=True)

        for date in unique_dates:
            fmt_date = date.strftime("%Y%m%d")
            file_name = os.path.join(folder_path, f"{fmt_date}_{suffix}.pqt")
            sub_df = df[temp_dates_str == fmt_date]

            if sub_df.empty:
                print(f"Warning: no data for {fmt_date}, skip")
                continue

            if not os.path.exists(file_name):
                empty_df.to_parquet(file_name, engine=PARQUET_ENGINE, compression=COMPRESSION)

            sub_df.to_parquet(
                file_name,
                engine=PARQUET_ENGINE,
                compression=COMPRESSION,
                append=True
            )

    except Exception as e:
        print("save_pqt_daily failed:", str(e))
        sys.exit(1)
