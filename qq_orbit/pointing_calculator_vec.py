#basic libraries
import numpy as np

#astropy lib for astro locations
from astropy import time as AT
from astropy.coordinates import  get_body
from scipy.spatial.transform import Rotation as Rot

#calculate current sheet location
from geopack import geopack, t96

RE=6371
AU = 149597870700/1e3

#self defined functions
#va, vb, N*3 numpy array
#return angle in radian
def angle_bwteen_vec(va,vb):
    return np.arccos(np.einsum('ij,ij->i,',va,vb)/(np.linalg.norm(va,axis=1)*np.linalg.norm(vb,axis=1)))

#rot x,y,z, use euler angle, see new gio/coord lib
#always from earth looking mode, with long edge (mechanic X axis) parallel to Equator plane
#return new pointing, to Y=0 plane
#angle is N array
def gen_rot(u0, angle):
    q1 = u0[0]*np.sin(angle/2)
    q2 = u0[1]*np.sin(angle/2)
    q3 = u0[2]*np.sin(angle/2)
    w = np.cos(angle/2)
    quat_ary = np.column_stack([q1, q2, q3, w])
    return Rot.from_quat(quat_ary)

def rot_inertia(angle,minus_z_vec,axis):
    Rot = gen_rot(axis,angle)
    new_minus_z  =  Rot.as_matrix()@minus_z_vec
    return new_minus_z

def find_intersection_with_Y0(direction_vec, start_point):
    """
    计算由起点和方向矢量定义的直线与Y=0平面的交点
    
    参数:
        direction_vec: 3x1 numpy数组，表示直线的方向矢量
        start_point: 3x1 numpy数组，表示直线的起点（默认原点）
        
    返回:
        intersection: 3x1 numpy数组，表示交点坐标（x, 0, z）
        若直线与Y=0平面平行或在平面上，则返回None
    """
    # 确保输入是3x1的数组
    direction_vec = np.array(direction_vec).reshape(3, 1)
    start_point = np.array(start_point).reshape(3, 1)
    
    # 提取方向矢量的Y分量和起点的Y分量
    dir_y = direction_vec[1, 0]
    start_y = start_point[1, 0]
    
    # 检查直线是否与Y=0平面平行（方向矢量Y分量为0）
    if np.isclose(dir_y, 0):
        if np.isclose(start_y, 0):
            print("直线在Y=0平面上，有无穷多个交点")
        else:
            print("直线与Y=0平面平行，无交点")
        return None
    
    # 计算参数t：直线方程为 start_point + t*direction_vec，求Y=0时的t
    # 即 start_y + t*dir_y = 0 → t = -start_y / dir_y
    t = -start_y / dir_y
    
    # 计算交点（此时Y分量必为0）
    intersection = start_point + t * direction_vec
    return intersection[0,0],intersection[2,0]


from spacepy import coordinates as coord
from spacepy.time import Ticktock
def convert(isot,value,in_type='ECI2000',out_type='GSM'):
    temp = coord.Coords(value,in_type,'car',ticks=Ticktock(isot)).convert(out_type,'car')
    return np.array([temp.x,temp.y,temp.z]).flatten()

#Calculation Flow
#ICRS coord
#'2025-09-09 00:00:00'
#tailward is +, headward is -
#northward is +, southward is -
def print_info(isot_str,angle_tail,angle_updown,loc_sat=None):
    
    isot = AT.Time(isot_str)
    
    loc_sat=None
    #use moon location to approximate queqiao location
    if loc_sat == None:
        moon = get_body('moon',isot).represent_as('cartesian')
        loc_sat = moon.xyz.value*AU/RE
    loc_sat_gsm = convert(isot,loc_sat)
    
    print('----------------------------------------------')
    print('UTC Time is', isot_str)
    #decide the sign of rotation angle based on lunar date
    #to do
    
    
    #decide the sign of rotation angle based on GSM_Y
    moon_half = 0
    if loc_sat_gsm[1]<0:
        moon_half = 1
        print('Second half of the month, A in the North')
    else:
        print('First half of the month, A in the South')
    if moon_half == 0:
        new_angle_tail = -angle_tail
        new_angle_updown = -angle_updown
    else:
        new_angle_tail = angle_tail
        new_angle_updown = angle_updown

    
    print('QueQiao location is (GSM)', loc_sat_gsm)

    loc_sun = get_body('sun',isot).represent_as('cartesian').xyz.value*AU/RE
    sun_sat_vec = loc_sun - loc_sat

    #By default, looking at Earth
    loc_point=np.array([0,0,0])
    minus_z_vec  = loc_point - loc_sat
    A_center_vec  = 

    minus_z_vec_gsm = convert(isot,minus_z_vec)
    X0,Z0 = find_intersection_with_Y0(minus_z_vec_gsm,loc_sat_gsm)
    print("Initial minus z pointing at GSM Y0 plane at", X0, Z0, "RE")
    sunavoid = angle_bwteen_vec(sun_sat_vec,minus_z_vec)
    print("Initial Sun avoid angle is", np.degrees(sunavoid), "degs")
    print('----------------------------------------------')

    #---------- rotation along ICRS Z ---
    minus_z_vec = rot_inertia(np.radians(new_angle_tail),minus_z_vec,np.array([0,0,1]))
    minus_z_vec_gsm = convert(isot,minus_z_vec)
    X0,Z0 = find_intersection_with_Y0(minus_z_vec_gsm,loc_sat_gsm)

    print("After rotate tailward ", str(new_angle_tail), "degs, minus z pointing GSM Y0 plane at", X0, Z0, "RE")
    sunavoid = angle_bwteen_vec(sun_sat_vec,minus_z_vec)
    print("Sun avoid angle is", np.degrees(sunavoid), "degs")
    print('----------------------------------------------')

    #---------- rotation along sun-earth line ---
    sun_earth = loc_sun - np.array([0,0,0])
    norm_sun_earth = sun_earth/np.linalg.norm(sun_earth)
    minus_z_vec = rot_inertia(np.radians(new_angle_updown),minus_z_vec,norm_sun_earth)

    minus_z_vec_gsm = convert(isot,minus_z_vec)
    X0,Z0 = find_intersection_with_Y0(minus_z_vec_gsm,loc_sat_gsm)
    print("After rotate northward ",str(new_angle_updown), "degs, minus z pointing GSM Y0 plane at", X0, Z0, "RE")

    sunavoid = angle_bwteen_vec(sun_sat_vec,minus_z_vec)
    print("Sun avoid angle is", np.degrees(sunavoid), "degs")
    print('----------------------------------------------')




print_info('2025-09-05 00:00:00',35,0)
print_info('2025-09-09 00:00:00',35,0)
print_info('2025-09-10 00:00:00',35,0)
print_info('2025-09-11 00:00:00',30,-2)
print_info('2025-09-12 00:00:00',30,-2)