output has 6 levels, 
1. L0: parameters directly converted from raw binary.
2. L1: measured values (tdc time and secondary electron amplitude) converted from L0
3. L1A: 8 anodes ((x1 x2 y1 y2) * (start end)) signal not all zero, splited A and B prob
4. L1B: 8 anodes ((x1 x2 y1 y2) * (start end)) signal all nonzero, splited A and B prob
5. L1C: 8 anodes ((x1 x2 y1 y2) * (start end)) four start singal (x1 x2 y1 y2) (start) all nonzero, splited A and B prob
6. L2: tdc measure time (ns) to coordinates (mm). obtained from L1_2_L2 function. (default is L1C to L2)

