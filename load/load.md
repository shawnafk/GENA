load.api
    负责访问数据源，得到 source variable

load.loaders
    负责解释 variable 的产品含义，并组装成具体 Frame

load.facade
    负责统一对外入口 load_data、load_pos、load_rot、load_sw，完成参数分发和常用组合加载
