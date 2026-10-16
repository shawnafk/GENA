# -------------------------------------------- astro object -----------------------------
#astroobj
from jplephem.spk import SPK
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
import json
import numpy as np
import os
import pandas as pd
from pathlib import Path
import shutil
import subprocess
from urllib.parse import urlencode
from xml.etree import ElementTree as ET

from common.constants import RE

current_dir = os.path.dirname(os.path.abspath(__file__))
kernel = SPK.open(current_dir+'/de432s.bsp')
from astropy.coordinates import  get_body
from astropy.time import Time
from datetime import datetime
import pytz
format = '%Y-%m-%d %H:%M:%S%z'
#dateime object has timezone information
#datetime.strptime(time_str, FORMAT)
#datetime.fromtimestamp(timestamp)
#dtobj.timestamp()
#dtobj.strftime('%Y-%m-%d %H:%M:%S')
#create Datetime object from str, timestamp to utc 1970-01-01 00:00:00
def Str_2_Datetime(str,tz='+0000'):
    return datetime.strptime(str+tz, format)

def Timestamp_2_Datetime(ts,tz='UTC'):
    return datetime.fromtimestamp(ts,pytz.timezone(tz))

def Datetime_2_Str(dt):
    return dt.strftime(format)

#convert timestamp to iso time
def Datetime_2_ISOT(dt):
    isot = Time(dt, scale='utc')
    return isot


#def def_timeinput(func):
#    def change_t(t):
#        if t is timestamp:
#
#        elif t is isot_str:
#        
#        elif t is ...:
#        
#        return func(t)
#    return change_t


#given timestamp get moon localtion in J2000
def moon(t):
    if np.isscalar(t) or isinstance(t, pd.Timestamp):
        t_iso = Datetime_2_ISOT(t)
        t_jd = t_iso.jd
    else:
        #given any t, to datetime
        #take array or scalar
        t_iso = [Datetime_2_ISOT(_) for _ in t]
        t_jd =  [_.jd for _ in t_iso]
    moon_bary = kernel[3,301].compute(t_jd).transpose()
    return moon_bary

def sun(t):
    if np.isscalar(t) or isinstance(t, pd.Timestamp):
        t_iso = Datetime_2_ISOT(t)
        t_jd = t_iso.jd
    else:
        t_iso = [Datetime_2_ISOT(_) for _ in t]
        t_jd =  [_.jd for _ in t_iso]
    earth_bary = kernel[0,3].compute(t_jd).transpose()
    sun_bary = -earth_bary
    return sun_bary

def myget_body(t,body='moon'):
    #m to km
    isot = [Datetime_2_ISOT(_) for _ in t]
    AU = 149597870700/1e3
    N = isot.shape[0]
    body_ej2000_xyz = np.zeros((N,3))
    #t = Time(iso_format_string, scale='utc')
    #t = Time(iso_format_string, format='isot', scale='utc')
    for i in range(N):
        body_ej2000 = get_body(body,isot[i]).represent_as('cartesian')
        body_ej2000_xyz[i,0] = np.double(body_ej2000.x)*AU
        body_ej2000_xyz[i,1] = np.double(body_ej2000.y)*AU
        body_ej2000_xyz[i,2] = np.double(body_ej2000.z)*AU
    #scj2000 = sc.transform_to(GCRS(obstime=Time("J2000"))).represent_as('cartesian')
    return body_ej2000_xyz

from scipy.interpolate import interp1d

_ORBIT_PATH = Path(current_dir) / 'orbits_all.npz'
_ORBIT_METADATA_PATH = Path(current_dir) / 'orbits_all.meta.json'
_SSC_ENDPOINT = 'https://sscweb.gsfc.nasa.gov/WS/sscr/2'
_SSC_NAMESPACE = {'ssc': 'http://sscweb.gsfc.nasa.gov/schema'}
_SSC_SOURCES = {
    'THA': (('themisa', 2), ('themisapred', 2)),
    'ATMS1': (('artemisp1', 2), ('artemisp1pred', 2)),
    'ATMS2': (('artemisp2', 2), ('artemisp2pred', 2)),
    'THD': (('themisd', 2), ('themisdpred', 2)),
    'THE': (('themise', 2), ('themisepred', 2)),
    'MMS1': (('mms1', 16),),
}


def _ssc_xml(path, *, timeout=120, attempts=4):
    """Return one SSCWeb XML response using its documented curl transport."""
    curl = shutil.which('curl')
    if curl is None:
        raise RuntimeError('SSCWeb updates require the curl executable')
    url = f'{_SSC_ENDPOINT}/{path}'
    command = [
        curl,
        '--fail',
        '--silent',
        '--show-error',
        '--location',
        '--retry',
        str(max(attempts - 1, 0)),
        '--retry-all-errors',
        '--retry-delay',
        '2',
        '--retry-max-time',
        str(timeout),
        '--max-time',
        str(timeout),
        '--header',
        'Accept: application/xml',
        '--header',
        'User-Agent: GENA-CODE-Dev ephemeris updater',
        url,
    ]
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            timeout=timeout * max(attempts, 1),
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        detail = getattr(exc, 'stderr', b'').decode(errors='replace').strip()
        raise RuntimeError(f'SSCWeb request failed: {url} ({detail})') from exc
    try:
        return ET.fromstring(result.stdout)
    except ET.ParseError as exc:
        raise RuntimeError(f'SSCWeb returned invalid XML: {url}') from exc


def _ssc_observatory_coverage(*, timeout=120):
    """Return SSC observatory UTC availability intervals."""
    root = _ssc_xml('observatories', timeout=timeout)
    coverage = {}
    for item in root.findall('.//ssc:Observatory', _SSC_NAMESPACE):
        identifier = item.findtext('ssc:Id', namespaces=_SSC_NAMESPACE)
        start = item.findtext('ssc:StartTime', namespaces=_SSC_NAMESPACE)
        stop = item.findtext('ssc:EndTime', namespaces=_SSC_NAMESPACE)
        if identifier and start and stop:
            coverage[identifier] = (
                pd.Timestamp(start).tz_convert('UTC'),
                pd.Timestamp(stop).tz_convert('UTC'),
            )
    return coverage


def _format_ssc_time(value):
    return pd.Timestamp(value).tz_convert('UTC').strftime('%Y%m%dT%H%M%SZ')


def _download_ssc_locations(
    observatory,
    start,
    stop,
    *,
    resolution_factor,
    timeout=120,
):
    """Download one SSC GSE trajectory and return Unix seconds plus ``R_E``."""
    interval = f'{_format_ssc_time(start)},{_format_ssc_time(stop)}'
    query = urlencode({'resolutionFactor': int(resolution_factor)})
    root = _ssc_xml(
        f'locations/{observatory}/{interval}/gse/?{query}',
        timeout=timeout,
    )
    status = root.findtext('.//ssc:StatusCode', namespaces=_SSC_NAMESPACE)
    if status and status.lower() != 'success':
        detail = root.findtext('.//ssc:StatusSubCode', namespaces=_SSC_NAMESPACE)
        raise RuntimeError(f'SSCWeb returned {status}/{detail} for {observatory}')
    data = root.find('.//ssc:Data', _SSC_NAMESPACE)
    if data is None:
        return np.empty((0, 4), dtype=float)
    returned_id = data.findtext('ssc:Id', namespaces=_SSC_NAMESPACE)
    if returned_id != observatory:
        raise RuntimeError(f'SSCWeb returned {returned_id!r} for {observatory!r}')
    coordinates = data.find('ssc:Coordinates', _SSC_NAMESPACE)
    if coordinates is None:
        return np.empty((0, 4), dtype=float)
    coordinate_system = coordinates.findtext(
        'ssc:CoordinateSystem', namespaces=_SSC_NAMESPACE
    )
    if coordinate_system is None or coordinate_system.lower() != 'gse':
        raise RuntimeError(
            f'SSCWeb returned {coordinate_system!r} coordinates for {observatory}'
        )
    times = [node.text for node in data.findall('ssc:Time', _SSC_NAMESPACE)]
    axes = [
        [float(node.text) for node in coordinates.findall(f'ssc:{axis}', _SSC_NAMESPACE)]
        for axis in 'XYZ'
    ]
    if not times:
        return np.empty((0, 4), dtype=float)
    if any(len(axis) != len(times) for axis in axes):
        raise RuntimeError(f'SSCWeb returned inconsistent coordinate lengths for {observatory}')
    dates = pd.DatetimeIndex(pd.to_datetime(times, utc=True)).as_unit('ns')
    timestamps = dates.asi8.astype(float) / 1e9
    positions_re = np.column_stack(axes) / (RE / 1000)
    result = np.column_stack((timestamps, positions_re))
    valid = np.isfinite(result).all(axis=1)
    result = result[valid]
    if len(result) and np.max(np.linalg.norm(result[:, 1:], axis=1)) > 100:
        raise RuntimeError(f'SSCWeb returned nonphysical {observatory} coordinates')
    return result


def _download_ssc_interval(
    observatory,
    start,
    stop,
    *,
    resolution_factor,
    chunk_days,
    timeout,
    verbose,
):
    chunks = []
    cursor = pd.Timestamp(start).tz_convert('UTC')
    stop = pd.Timestamp(stop).tz_convert('UTC')
    delta = pd.Timedelta(days=chunk_days)
    while cursor <= stop:
        chunk_stop = min(cursor + delta, stop)
        if verbose:
            print(
                f'download SSCWeb {observatory}: {cursor.isoformat()} to '
                f'{chunk_stop.isoformat()}',
                flush=True,
            )
        chunks.append(
            _download_ssc_locations(
                observatory,
                cursor,
                chunk_stop,
                resolution_factor=resolution_factor,
                timeout=timeout,
            )
        )
        if chunk_stop >= stop:
            break
        cursor = chunk_stop
    chunks = [chunk for chunk in chunks if len(chunk)]
    return np.vstack(chunks) if chunks else np.empty((0, 4), dtype=float)


def _merge_orbit_rows(*arrays):
    arrays = [np.asarray(array, dtype=float) for array in arrays if len(array)]
    if not arrays:
        return np.empty((0, 4), dtype=float)
    merged = np.vstack(arrays)
    merged = merged[np.argsort(merged[:, 0], kind='stable')]
    _, reverse_indices = np.unique(merged[::-1, 0], return_index=True)
    keep = np.sort(len(merged) - 1 - reverse_indices)
    return merged[keep]


def _load_orbit_metadata(path):
    if not path.is_file():
        return {'keys': {}}
    try:
        return json.loads(path.read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError):
        return {'keys': {}}


def _update_one_orbit_key(
    key,
    existing,
    target,
    availability,
    key_metadata,
    *,
    chunk_days,
    timeout,
    verbose,
):
    cadence_seconds = float(np.median(np.diff(existing[:, 0])))
    cadence = pd.Timedelta(seconds=cadence_seconds)
    key_metadata = dict(key_metadata or {})
    key_metadata['sources'] = dict(key_metadata.get('sources', {}))
    primary_end = pd.to_datetime(existing[-1, 0], unit='s', utc=True)
    additions = []
    for source_index, (observatory, resolution_factor) in enumerate(_SSC_SOURCES[key]):
        if observatory not in availability:
            raise RuntimeError(f'SSCWeb observatory is unavailable: {observatory}')
        available_start, available_stop = availability[observatory]
        source_record = key_metadata['sources'].get(observatory, {})
        recorded_stop = source_record.get('downloaded_through')
        if source_index == 0:
            source_cursor = (
                pd.Timestamp(recorded_stop).tz_convert('UTC') + cadence
                if recorded_stop
                else min(primary_end, available_stop) + cadence
            )
        else:
            source_cursor = primary_end + cadence
        source_start = max(source_cursor, available_start)
        source_stop = min(target + cadence, available_stop)
        if source_start > source_stop:
            continue
        downloaded = _download_ssc_interval(
            observatory,
            source_start,
            source_stop,
            resolution_factor=resolution_factor,
            chunk_days=chunk_days,
            timeout=timeout,
            verbose=verbose,
        )
        if not len(downloaded):
            continue
        additions.append(downloaded)
        downloaded_stop = pd.to_datetime(downloaded[-1, 0], unit='s', utc=True)
        key_metadata['sources'][observatory] = {
            'downloaded_through': downloaded_stop.isoformat(),
            'resolution_factor': resolution_factor,
        }
        if source_index == 0:
            primary_end = downloaded_stop
    updated = _merge_orbit_rows(existing, *additions)
    key_metadata['coverage_start'] = pd.to_datetime(
        updated[0, 0], unit='s', utc=True
    ).isoformat()
    key_metadata['coverage_end'] = pd.to_datetime(
        updated[-1, 0], unit='s', utc=True
    ).isoformat()
    if updated[-1, 0] < target.timestamp():
        raise RuntimeError(
            f'SSCWeb could only update {key} through '
            f'{key_metadata["coverage_end"]}, short of {target.isoformat()}'
        )
    return updated, key_metadata


def update_satellite_orbits(
    keys=None,
    *,
    end=None,
    orbit_path=None,
    metadata_path=None,
    chunk_days=31,
    timeout=120,
    workers=1,
    verbose=False,
):
    """Extend the local SSCWeb GSE orbit archive atomically.

    Definitive THEMIS/ARTEMIS data is preferred. If it does not yet reach
    ``end``, the corresponding SSC predictive observatory fills the remaining
    interval. Coordinates downloaded in kilometres are stored in ``R_E``.
    """
    orbit_path = Path(orbit_path or _ORBIT_PATH)
    metadata_path = Path(metadata_path or _ORBIT_METADATA_PATH)
    if keys is None:
        keys = list(_SSC_SOURCES)
    elif isinstance(keys, str):
        keys = [keys]
    keys = list(dict.fromkeys(str(key).upper() for key in keys))
    unknown = [key for key in keys if key not in _SSC_SOURCES]
    if unknown:
        raise KeyError(f'unsupported SSC orbit key(s): {unknown}')
    if not isinstance(chunk_days, int) or chunk_days < 1:
        raise ValueError('chunk_days must be a positive integer')
    if not isinstance(workers, int) or workers < 1:
        raise ValueError('workers must be a positive integer')
    target = pd.Timestamp.now(tz='UTC') if end is None else pd.Timestamp(end)
    target = target.tz_localize('UTC') if target.tzinfo is None else target.tz_convert('UTC')

    with np.load(orbit_path) as archive:
        datasets = {name: np.asarray(archive[name], dtype=float) for name in archive.files}
    missing = [key for key in keys if key not in datasets]
    if missing:
        raise KeyError(f'orbit archive is missing key(s): {missing}')
    keys_to_update = [key for key in keys if target.timestamp() > datasets[key][-1, 0]]
    if not keys_to_update:
        return {
            key: (
                pd.to_datetime(datasets[key][0, 0], unit='s', utc=True),
                pd.to_datetime(datasets[key][-1, 0], unit='s', utc=True),
            )
            for key in keys
        }

    availability = _ssc_observatory_coverage(timeout=timeout)
    metadata = _load_orbit_metadata(metadata_path)
    metadata.setdefault('keys', {})
    update_kwargs = {
        'target': target,
        'availability': availability,
        'chunk_days': chunk_days,
        'timeout': timeout,
        'verbose': verbose,
    }
    if workers == 1 or len(keys_to_update) == 1:
        for key in keys_to_update:
            datasets[key], metadata['keys'][key] = _update_one_orbit_key(
                key,
                datasets[key],
                key_metadata=metadata['keys'].get(key),
                **update_kwargs,
            )
    else:
        with ThreadPoolExecutor(max_workers=min(workers, len(keys_to_update))) as executor:
            futures = {
                executor.submit(
                    _update_one_orbit_key,
                    key,
                    datasets[key],
                    key_metadata=metadata['keys'].get(key),
                    **update_kwargs,
                ): key
                for key in keys_to_update
            }
            for future in as_completed(futures):
                key = futures[future]
                datasets[key], metadata['keys'][key] = future.result()

    orbit_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = orbit_path.with_name(f'.{orbit_path.name}.tmp')
    try:
        with temporary.open('wb') as stream:
            np.savez(stream, **datasets)
        temporary.replace(orbit_path)
    finally:
        temporary.unlink(missing_ok=True)
    metadata.update({
        'updated_at': pd.Timestamp.now(tz='UTC').isoformat(),
        'source': _SSC_ENDPOINT,
        'coordinate_system': 'GSE',
        'unit': 'R_E',
    })
    metadata_temporary = metadata_path.with_name(f'.{metadata_path.name}.tmp')
    try:
        metadata_temporary.write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + '\n',
            encoding='utf-8',
        )
        metadata_temporary.replace(metadata_path)
    finally:
        metadata_temporary.unlink(missing_ok=True)
    _satellite_interpolators.cache_clear()
    return {
        key: (
            pd.to_datetime(datasets[key][0, 0], unit='s', utc=True),
            pd.to_datetime(datasets[key][-1, 0], unit='s', utc=True),
        )
        for key in keys
    }


@lru_cache(maxsize=None)
def _satellite_interpolators(key):
    """Load one satellite's GSE ``R_E`` orbit and build cached interpolators."""
    with np.load(_ORBIT_PATH, 'rb') as orbits:
        o_dat = orbits[key]
    interpolators = tuple(
        interp1d(
            o_dat[:, 0],
            o_dat[:, component],
            kind='cubic',
            fill_value="extrapolate",
        )
        for component in (1, 2, 3)
    )
    return (*interpolators, float(o_dat[0, 0]), float(o_dat[-1, 0]))


def satellites(t_eval, key='MMS1', *, extrapolate=False, update=False):
    """Interpolate GSE satellite positions in Earth radii.

    Parameters
    ----------
    t_eval : array-like
        Unix timestamps in seconds.
    key : str, optional
        Orbit key stored in ``orbits_all.npz``.
    extrapolate : bool, optional
        Permit interpolation outside the stored orbit coverage. Disabled by
        default because cubic extrapolation can produce nonphysical positions.
    update : bool, optional
        If requested timestamps extend beyond local coverage, download and
        atomically merge SSCWeb data through the latest requested time before
        rebuilding the interpolation cache. Network access is opt-in.

    Returns
    -------
    numpy.ndarray
        Array with shape ``(N, 3)`` in GSE coordinates and ``R_E`` units.
    """
    key = str(key).upper()
    fx, fy, fz, start, stop = _satellite_interpolators(key)
    values = np.atleast_1d(np.asarray(t_eval, dtype=float))
    if not np.isfinite(values).all():
        raise ValueError('satellite ephemeris timestamps must be finite')
    outside = (values < start).any() or (values > stop).any()
    if outside and update:
        if (values < start).any():
            raise ValueError('automatic SSCWeb update only extends future coverage')
        update_satellite_orbits(key, end=pd.to_datetime(values.max(), unit='s', utc=True))
        fx, fy, fz, start, stop = _satellite_interpolators(key)
        outside = (values < start).any() or (values > stop).any()
    if not extrapolate and outside:
        requested_start = pd.to_datetime(values.min(), unit='s', utc=True)
        requested_stop = pd.to_datetime(values.max(), unit='s', utc=True)
        coverage_start = pd.to_datetime(start, unit='s', utc=True)
        coverage_stop = pd.to_datetime(stop, unit='s', utc=True)
        raise ValueError(
            f'{key} ephemeris covers {coverage_start} to {coverage_stop}; '
            f'requested {requested_start} to {requested_stop}'
        )
    return np.column_stack((fx(t_eval), fy(t_eval), fz(t_eval)))
