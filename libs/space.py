from geopack import geopack
import numpy as np
import pandas as pd
from scipy import optimize
#get key variables from omni at given time
#MP, solar wind dynamic pressure and IMP_Bz
#T96 field lines: solar wind dynamic pressure, dst index, IMF_By and IMF_Bz
#Bowshock: Plasma beta, Mag_Mach_Num, solar wind dynamic pressure, and IMF_Bz

#for 3 views, need  front: y-z at given x, side: x-z at x=0, top: x-y at z=0
#streamlines
#for 3d, need surfaces
class T96_fieldlines():
    def __init__(self,ut,swdp,dst,by,bz,Nlat=15,Nlong=5):
        self.ut  = ut
        self.swdp = swdp
        self.dst = dst
        self.by  = by
        self.bz  = bz
        self.Nlat = Nlat
        self.Nlong = Nlong
        self.xz = []
        self.xyz = []
        self.y0 = []
        #self.xz = self.trace_lines_xz()
        #self.xyz = self.trace_lines_xyz()
        #rbs = np.sqrt(bs[0,:]**2 + bs[1,:]**2)
        #idx_bs = rbs < rlim
        #self.lines = [bs[0,idx_bs],bs[1,idx_bs]]
        #self.xz = self.trace_lines_xz()
    #solar wind pressure pdyn (nanopascals)
    #dst (nanotesla)
    #byimf (nanotesla)
    #bzimf (nanotesla)
    def trace_lines_3d(self):
        ut = self.ut
        pm=[self.swdp,self.dst,self.by,self.bz]
        geopack.recalc(ut)
        #mp
        x=[]
        y=[]
        z=[]
        xyz=[]
        #high latitude, more
        #for lat in np.concatenate((np.linspace(-88,-60,10),np.linspace(-60,60,5),np.linspace(60,88,10))):
        for lat in np.linspace(-90,90,self.Nlat):
            for long in np.linspace(-180,180,self.Nlong):
                x.append(3.5*np.cos(np.radians(long))*np.sin(np.radians(lat)))
                y.append(3.5*np.sin(np.radians(long))*np.sin(np.radians(lat)))
                z.append(3.5*np.cos(np.radians(lat)))

        for dir in [-1,1]:
            for _x,_y,_z in zip(x,y,z):
                res =  geopack.trace(_x,_y,_z,dir,rlim=100,parmod=pm,exname='t96')
                xyz.append([res[3],res[4],res[5]])
        return xyz
    
    def trace_lines_xz(self):
        ut = self.ut
        pm=[self.swdp,self.dst,self.by,self.bz]
        geopack.recalc(ut)
        x=[]
        y=[]
        z=[]
        #meridional plane
        for lat in np.linspace(-90,90,self.Nlat):
            for long in [0, 180]:
                x.append(3.5*np.cos(np.radians(long))*np.sin(np.radians(90-lat)))
                y.append(3.5*np.sin(np.radians(long))*np.sin(np.radians(90-lat)))
                z.append(3.5*np.cos(np.radians(90-lat)))
        #side
        xz = []
        #given x plane
        for dir in [-1,1]:
                for _x,_y,_z in zip(x,y,z):
                    res =  geopack.trace(_x,_y,_z,dir,rlim=100,parmod=pm,exname='t96')
                    xz.append([res[3],res[5]])
        return xz

    def trace_lines_y0(
        self,
        *,
        seed_radius=1.05,
        footprint_latitude_range=(70.0, 82.0),
        rlim=80,
    ):
        """Trace T96 tail field lines from nightside magnetic footpoints.

        The returned list contains ``(3, N)`` arrays in GSM Earth radii. The
        seeds sample northern nightside SM footpoint latitudes densely across
        ``footprint_latitude_range``. They are converted to GSM after
        :func:`geopack.recalc`, then traced outward through the magnetosphere.
        This avoids duplicating the same physical line from both hemispheres
        and concentrates ``Nlat`` samples on field lines that approach or
        cross the magnetotail current sheet. Nonzero IMF ``By`` may carry the
        resulting three-dimensional lines away from GSM ``Y=0``.
        """
        parameters = np.asarray(
            [self.swdp, self.dst, self.by, self.bz],
            dtype=float,
        )
        if not np.isfinite(parameters).all():
            raise ValueError("T96 parameters must be finite")
        if not np.isfinite(seed_radius) or seed_radius <= 1:
            raise ValueError("seed_radius must be finite and greater than 1 R_E")
        if not np.isfinite(rlim) or rlim <= seed_radius:
            raise ValueError("rlim must be finite and greater than seed_radius")
        if not isinstance(self.Nlat, (int, np.integer)) or self.Nlat < 1:
            raise ValueError("Nlat must be a positive integer")
        latitude_range = np.asarray(footprint_latitude_range, dtype=float)
        if (
            latitude_range.shape != (2,)
            or not np.isfinite(latitude_range).all()
            or not 0 < latitude_range[0] < latitude_range[1] < 90
        ):
            raise ValueError(
                "footprint_latitude_range must be two increasing values "
                "between 0 and 90 degrees"
            )

        geopack.recalc(float(self.ut))
        lines = []
        for latitude in np.linspace(*latitude_range, self.Nlat):
            latitude_rad = np.radians(latitude)
            seed_sm = (
                -seed_radius * np.cos(latitude_rad),
                0.0,
                seed_radius * np.sin(latitude_rad),
            )
            seed_gsm = geopack.smgsm(*seed_sm, 1)
            result = geopack.trace(
                *seed_gsm,
                1,
                rlim=rlim,
                r0=1,
                parmod=parameters,
                exname="t96",
                maxloop=2000,
            )
            line = np.asarray(result[3:6], dtype=float)
            valid = np.isfinite(line).all(axis=0)
            line = line[:, valid]
            if line.shape[1] >= 2:
                lines.append(line)
        self.y0 = lines
        return lines

    def draw3d(self, ax,**kwargs):
        if self.xyz == []:
            self.xyz = self.trace_lines_3d()
        for l in self.xyz:
            ax.plot(l[0], l[1], l[2],**kwargs)

    def draw_xz(self, ax, **kwargs):
        #if self.xyz == []:
        #    self.xyz = self.trace_lines_3d()
        if self.xz == []:
            self.xz = self.trace_lines_xz()
        for l in self.xz:
            ax.plot(l[0], l[1],**kwargs)
        #for l in self.xyz:
        #    ax.scatter(l[0], l[2], color = c, s = 0.5, **kwargs)

    def draw_xy(self, ax, **kwargs):
        if self.xyz == []:
            self.xyz = self.trace_lines_3d()
        for l in self.xyz:
            ax.scatter(l[0], l[1], s = 0.5, **kwargs)

    def draw_yz(self, ax, **kwargs):
        if self.xyz == []:
            self.xyz = self.trace_lines_3d()
        for l in self.xyz:
            ax.scatter(l[1], l[2], s = 0.5, **kwargs)


class T96NeutralSheet:
    """Numerical T96+IGRF ``Bx=0`` neutral-sheet surface in GSM coordinates.

    The central surface is found independently at every ``(X, Y)`` grid point
    by bracketing the root of the total GSM ``Bx`` component along ``Z``. The
    optional thickness is a visual half-thickness in ``R_E``; it is not a
    density-derived plasma-sheet boundary.
    """

    def __init__(
        self,
        ut,
        swdp,
        dst,
        by,
        bz,
        *,
        x_range=(-70.0, -5.0),
        y_range=(-30.0, 30.0),
        z_range=(-20.0, 20.0),
        nx=33,
        ny=25,
        nz=17,
        half_thickness=3.0,
    ):
        self.ut = float(ut)
        self.swdp = float(swdp)
        self.dst = float(dst)
        self.by = float(by)
        self.bz = float(bz)
        self.x_range = self._validate_range('x_range', x_range)
        self.y_range = self._validate_range('y_range', y_range)
        self.z_range = self._validate_range('z_range', z_range)
        self.nx = self._validate_count('nx', nx)
        self.ny = self._validate_count('ny', ny)
        self.nz = self._validate_count('nz', nz)
        self.half_thickness = float(half_thickness)
        if not np.isfinite(self.half_thickness) or self.half_thickness <= 0:
            raise ValueError('half_thickness must be a positive finite value')
        parameters = np.asarray([self.swdp, self.dst, self.by, self.bz])
        if not np.isfinite(self.ut) or not np.isfinite(parameters).all():
            raise ValueError('T96 neutral-sheet parameters must be finite')
        self._surface = None
        self._psi = None

    @staticmethod
    def _validate_range(name, values):
        values = np.asarray(values, dtype=float)
        if (
            values.shape != (2,)
            or not np.isfinite(values).all()
            or values[0] >= values[1]
        ):
            raise ValueError(f'{name} must contain two increasing finite values')
        return tuple(values)

    @staticmethod
    def _validate_count(name, value):
        if not isinstance(value, (int, np.integer)) or value < 2:
            raise ValueError(f'{name} must be an integer of at least 2')
        return int(value)

    @property
    def parameters(self):
        return np.asarray([self.swdp, self.dst, self.by, self.bz], dtype=float)

    @property
    def dipole_tilt_radians(self):
        if self._psi is None:
            self._psi = float(geopack.recalc(self.ut))
        return self._psi

    @property
    def dipole_tilt_degrees(self):
        return float(np.degrees(self.dipole_tilt_radians))

    def magnetic_field(self, x, y, z):
        """Return total T96+IGRF GSM field components in nT at one point."""
        psi = self.dipole_tilt_radians
        external = geopack.call_external_model(
            't96',
            self.parameters,
            psi,
            float(x),
            float(y),
            float(z),
        )
        internal = geopack.igrf_gsm(float(x), float(y), float(z))
        return np.asarray(external, dtype=float) + np.asarray(internal, dtype=float)

    def bx(self, x, y, z):
        """Return the total T96+IGRF GSM ``Bx`` component in nT."""
        return float(self.magnetic_field(x, y, z)[0])

    def _z_root(self, x, y):
        samples = np.linspace(*self.z_range, self.nz)
        values = np.asarray([self.bx(x, y, z) for z in samples])
        roots = []
        for index, value in enumerate(values):
            if value == 0:
                roots.append(samples[index])
            if index == len(samples) - 1:
                continue
            next_value = values[index + 1]
            if value * next_value < 0:
                roots.append(optimize.brentq(
                    lambda z: self.bx(x, y, z),
                    samples[index],
                    samples[index + 1],
                ))
        if not roots:
            return np.nan
        return float(min(roots, key=abs))

    def surface_xyz(self):
        """Return ``(X, Y, Z)`` GSM grids for the central ``Bx=0`` surface."""
        if self._surface is None:
            x = np.linspace(*self.x_range, self.nx)
            y = np.linspace(*self.y_range, self.ny)
            X, Y = np.meshgrid(x, y)
            Z = np.empty_like(X)
            for row in range(self.ny):
                for column in range(self.nx):
                    Z[row, column] = self._z_root(X[row, column], Y[row, column])
            self._surface = (X, Y, Z)
        return self._surface

    def thickness_xyz(self):
        """Return lower and upper visual sheet bounds in GSM ``R_E``."""
        X, Y, Z = self.surface_xyz()
        return (X, Y, Z - self.half_thickness), (X, Y, Z + self.half_thickness)

    def z_at(self, x, y):
        """Return the central-sheet ``Z`` at a point within the rendered domain."""
        if not (
            self.x_range[0] <= x <= self.x_range[1]
            and self.y_range[0] <= y <= self.y_range[1]
        ):
            return np.nan
        return self._z_root(float(x), float(y))


#surfaces

#defined by r(theta, phi no) functions
#theta, determins r range
#phi, azimuthal angle, 0,2*pi
#common features

class R_surface:
    def __init__(self, center, latitude,longitude):
        self.center = center
        self.lat = latitude
        self.long = longitude
    
    def r_func(self,latitude):
        #to be defined in subclass
        pass
    
    def r2xyz(self,r,latitude,longitude):
        R = r * np.sin(latitude)
        x = self.center[0]+ r * np.cos(latitude)
        y = self.center[1]+ R * np.cos(longitude)
        z = self.center[2]+ R * np.sin(longitude)
        return x,y,z
    
    def surface_xyz(self):
        Lat, Long = np.meshgrid(self.lat, self.long)
        return self.r2xyz(self.r_func(Lat),Lat,Long)
    
    #longitude is perendicular to zy plane, zy plane is equtorial plane
    #cross section at y=0 plane
    #plane range given by lat range
    #long = pi/2 to get z
    def cut_ym0(self):
        r = self.r_func(self.lat)
        x,_,z = self.r2xyz(r,self.lat,np.pi/2)
        return x,z

    #long = pi/2 to get y
    def cut_zm0(self):
        r = self.r_func(self.lat)
        x,y,_ = self.r2xyz(r,self.lat,0)
        return x,y
    
    #cross section at x=d plane
    def cut_xd(self,d):
        def f(lat):
            return self.r_func(lat) * np.cos(lat) - d
        try:
            lat_cross = optimize.fsolve(f,0)
            if lat_cross == [0]:
                _,y,z = [[],[],[]]
            else:   
                _,y,z = self.r2xyz(self.r_func(lat_cross),lat_cross,self.long)
        except:
            _,y,z = [[],[],[]]
        return y,z
    
    #show method
    def draw3d(self, ax,alp = 0.5,**kwargs):
        xx,yy,zz = self.surface_xyz()
        ax.plot_surface(
            xx,yy,zz,
            #edgecolor=self.color,
            alpha=alp,
            **kwargs
            #cmap='viridis'  # 颜色映射（可选）
        )
    
    def draw_xz(self,ax,**kwargs):
        x,z = self.cut_ym0()
        ax.plot(x,z, **kwargs)
    
    def draw_xy(self,ax,**kwargs):
        x,y = self.cut_zm0()
        ax.plot(x,y, **kwargs)
        
    def draw_yz(self,ax,d,**kwargs):
        y,z = self.cut_xd(d)
        print(len(y))
        if len(y)!=0 or len(z) != 0:
            ax.plot(y,z, **kwargs)
        else: 
            print("no cross point at x=",d)
            _,yy,zz = self.surface_xyz()
            ax.scatter(yy.flatten(),zz.flatten(),s=.5)
            
            
class Star(R_surface):
    def __init__(self, center = [0,0,0], R = 1):
        self.center = center
        self.R = R
        lat = np.radians(np.linspace(-180, 180, 128))
        long = np.radians(np.linspace(-180, 180, 128))
        super().__init__(center,lat, long)

    def r_func(self,lat):
        r = self.R * np.ones_like(lat)
        return r

class Bowshock(R_surface):
    def __init__(self,B_z,beta,M_ms,D_p,scale=0.8):
        self.B_z  = B_z
        self.beta = beta
        self.M_ms = M_ms
        self.D_p  = D_p
        lat = np.radians(np.linspace(-180,180,32))*scale
        long = np.linspace(0, 2*np.pi, 15)
        super().__init__([0,0,0], lat, long)

    def r_func(self,theta):
        B_z  = self.B_z 
        beta = self.beta
        M_ms = self.M_ms
        D_p  = self.D_p 
 
        a_1=11.1266
        a_2=0.0010
        a_3=-0.0005
        a_4=2.5966
        a_5=0.8182
        a_6=-0.0170
        a_7=-0.0122
        a_8=1.3007
        a_9=-0.0049
        a_10=-0.0328
        a_11=6.047
        a_12=1.029
        a_13=0.0231
        a_14=-0.002
        if B_z >0:
            r_0=a_1 * (1+a_2*B_z)*(1+a_9 * beta ) * (1+a_4 * ((a_8-1)*M_ms**2+2)/ (a_8+1)/M_ms**2) * D_p**(-1 / a_11)
            alpha = a_5 * (1+a_13 * B_z)*(1+ a_7 * D_p)*(1+a_10*np.log(1+beta)) * (1+a_14 * M_ms)
        else:
            r_0=a_1 * (1+a_3*B_z)*(1+a_9 * beta ) * (1+a_4 * ((a_8-1)*M_ms**2+2)/ (a_8+1)/M_ms**2) * D_p**(-1/a_11)
            alpha=a_5*(1+a_6 * B_z)*(1+ a_7 * D_p)*(1+a_10*np.log(1+beta)) * (1+a_14 * M_ms)
        epsi = a_12
        r = r_0 * ( (1+epsi)/(1+epsi * np.cos(theta) ) ) ** alpha
        return r

class Magnetopause(R_surface):
    def __init__(self,B_z,D_p,scale=0.9):
        self.B_z  = B_z
        self.D_p  = D_p
        lat = np.radians(np.linspace(-180,180,32))*scale
        long = np.linspace(0, 2*np.pi, 16) 
        super().__init__([0,0,0],lat, long)
    
    def r_func(self,lat):
        B_z  = self.B_z 
        D_p  = self.D_p 
        a_1 = 11.646
        a_2 = 0.216
        a_3 = 0.122
        a_4 = 6.215
        a_5 = 0.578
        a_6 = -0.009
        a_7 = 0.012
        if B_z > 0:
            r0  = a_1 * D_p**(-1/a_4)
        elif B_z < 0 and B_z > -8:
            r0 = (a_1 + a_2* B_z) * D_p**(-1/a_4)
        else:
            r0 =( a_1 + 8*a_3 - 8 * a_2 + a_3 * B_z ) * D_p**(-1/a_4)
        alpha = (a_5 + a_6 * B_z) * (1+ a_7 * D_p)

        r = r0 * ( 2/(1+np.cos(lat) ) ) ** alpha
        return r

def get_obj(ISOT,DST=0,D_p=1,beta=5,M_ms=5,B_y=0,B_z=0):
    UT = pd.to_datetime(ISOT).timestamp()
    obj={}
    obj['earth'] = Star()
    obj['bs'] = Bowshock(B_z,beta,M_ms,D_p)
    obj['mp'] = Magnetopause(B_z,D_p)
    obj['fl'] = T96_fieldlines(UT, D_p, DST, B_y, B_z)
    return obj

from load import load_sw
# ---------------------------------------------------------------------------
# Solar wind parameter interpolation  (heavily used by cross-mission analysis)
# ---------------------------------------------------------------------------

def get_obj_ave(st_time,use_average_sw_params=True):
    if use_average_sw_params:
        #add verify and backup
        sw_paras = load_sw(st_time)
        ave_sw_paras = {k: np.mean(v) for k, v in sw_paras.items()}
        by = 0
        bz = ave_sw_paras['BZ_GSM']
        mgs = ave_sw_paras['Mgs_mach_num']
        beta = ave_sw_paras['Beta']
        dst = ave_sw_paras['SYM_H']
        dp = ave_sw_paras['Pressure'] 
        bs_ave = Bowshock(bz, beta, mgs, dp)
        mp_ave = Magnetopause(bz, dp)
        fl = T96_fieldlines(np.average(st_time),dp,dst,by,bz)
    else:
        by = 0
        bz = 0
        mgs = 8
        beta = 5
        dst = -40
        dp = 5 
        bs_ave = Bowshock(bz, beta, mgs, dp)
        mp_ave = Magnetopause(bz, dp)
        fl = T96_fieldlines(np.average(st_time),dp,dst,by,bz)
    obj={}
    obj['earth'] = Star()
    obj['bs'] = bs_ave
    obj['mp'] = mp_ave
    obj['fl'] = fl
    return obj

def objects(ax, obj, dir='xz', d=-60):
    if 'xz' == dir:
        obj['earth'].draw_xz(ax, color='w')
        obj['bs'].draw_xz(ax, color='r')
        obj['mp'].draw_xz(ax, color='k')
        obj['fl'].draw_xz(ax, color='w')
    if 'xy' == dir:
        obj['earth'].draw_xy(ax, color='w')
        obj['bs'].draw_xy(ax, color='r')
        obj['mp'].draw_xy(ax, color='k')
        obj['fl'].draw_xy(ax, color='w')
    if 'yz' == dir:
        obj['bs'].draw_yz(ax, d=d, color='r')
        obj['mp'].draw_yz(ax, d=d, color='k')
        obj['fl'].draw_yz(ax, color='w')



if __name__ == "__main__":
    print("unit test here")
