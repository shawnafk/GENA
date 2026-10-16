from gena import df as gdf
from gena import filters as flt
df = gdf.read_csv('./CE7ZJ_1_202407-202410.raw_L1A_a.csv')
subdf = flt.start_eff(flt.count_eff(df))
# True=缺失
missing_s = subdf.iloc[2000:3000][['X1_start (ns)', 'X2_start (ns)', 'Y1_start (ns)', 'Y2_start (ns)']] == 0 
missing_e = subdf.iloc[2000:3000][['X1_end (ns)', 'X2_end (ns)', 'Y1_end (ns)', 'Y2_end (ns)']] == 0

from scipy.stats import chi2_contingency
import pandas as pd
import numpy as np
# 检验 x1s 与 x2s 的独立性
ct = pd.crosstab(missing_e['Y1_end (ns)'], missing_e['Y2_end (ns)'])
chi2, p, dof, _ = chi2_contingency(ct)
cramers_v = np.sqrt(chi2 / (ct.sum().sum() * (min(ct.shape)-1)))
print(f"x1s vs y2s: p={p:.15e}, Cramér's V={cramers_v:.3f}")

# 同理检验 y1s vs y2s、x1s vs y1s 等
