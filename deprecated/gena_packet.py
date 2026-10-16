import numpy as np
def find_head_positions(f,len_range):
    headseek_num = []
    f.seek(len_range[0],0)
    remaining_bytes = len_range[1] - len_range[0];
    while (remaining_bytes > 3):
        if remaining_bytes < 4:
            read_size = remaining_bytes
        else:
            read_size = 4
         #read 4B every time
        headdata = f.read(read_size)
        # 寻找大包的包头 E2 25 0C 4E，以获取包头位置指针
        #if perfectly saved, all packet have 2064 length, packet would start with E2 (if missing 4)
        if (headdata[0] == 226) & (headdata[1] == 37) & (headdata[2] == 12) & (headdata[3] == 78):
            headseek_num.append(f.tell()-4)
        #if missing three
        elif (headdata[1] == 226) & (headdata[2] == 37) & (headdata[3] == 12):
            f.seek(-3,1)
        #same, if missing two
        elif (headdata[2] == 226) & (headdata[3] == 37):
            f.seek(-2,1)
        #however, if there is missing one Byte, then the last one of the 4 B would be E5, which belong to next 2064
        elif (headdata[3] == 226):
            #zhengjs, seek front -1 from current location
            f.seek(-1,1)
            #then re-read 4 to make the four always E2 25 0C 4E
            #put pointer just before E2
        remaining_bytes = len_range[1] - f.tell();
    return np.array(headseek_num)

def addr(f,n,head):
    f.seek(head[n])
    h = f.read(1)
    print('%#x'%h[0])
    f.seek(head[n])
    f.seek(73,1)
    p = f.read(1)
    print('%#x'%p[0])
    p4 = f.read(4)
    v = int.from_bytes(p4) & 0b111111111111
    print(v)



def get_packet_time(p):
    p.seek(8,1)
    ts = p.read(4)
    tms = p.read(2)
    Tim_s = int.from_bytes(ts)
    Tim_ms = int.from_bytes(tms)
    return Tim_s + Tim_ms/1000.

def get_all_A(f,hp):
    Time=[]
    A0 = []
    A1 = []
    A2 = []
    A3 = []
    for idx in hp:
        f.seek(idx)
        Time.append(get_packet_time(f))
        f.seek(idx+73)
        p = f.read(1)
        #print('%#x'%p[0])
        p4 = f.read(4)
        v = int.from_bytes(p4) & 0b111111111111
        #addr = 0x35
        if p[0] == 53:
            A0.append(v)
        #addr = 0x37
        if p[0] == 55:
            A1.append(v)
        #addr = 0x39
        if p[0] == 57:
            A2.append(v)
        #addr = 0x3B
        if p[0] == 59:
            A3.append(v)
        #!addr = 0x3D
        #if p[0] == 61:
        #    A3.append(v)
    return Time,A0,A1,A2,A3


def get_all_B(f,hp):
    Time=[]
    B0 = []
    B1 = []
    B2 = []
    B3 = []
    B= []
    for idx in hp:
        f.seek(idx)
        Time.append(get_packet_time(f))
        f.seek(idx+73)
        p = f.read(1)
        print('%#x'%p[0])
        v = f.read(1)
        if p[0] == 142:
            B0.append(v)
        #addr = 0x90
        if p[0] == 144:
            B1.append(v)
        if p[0] == 146:
            B2.append(v)
        if p[0] == 148:
            B3.append(v)
        if p[0] == 150:
            B.append(v)
    return Time,B2,B1,B2,B3,B

f = open('/home/zhengjs/GENA_RAW/ENA20240727','rb')
hp = find_head_positions(f,[0,f.seek(0,2)])
TA,A0,A1,A2,A3 = get_all_A(f,hp)
TB,B0,B1,B2,B3,B = get_all_B(f,hp)