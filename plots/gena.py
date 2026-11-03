import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde

from gena import filter as flt
from gena.constants import init_param
from gena.post import statistic as sts
from plots import tools
def _interpolate_nonzero_parameter(pdf, name):
    """Fill zero-valued engineering samples for plotting only."""
    timestamp = pdf["Timestamp"].to_numpy(dtype=float)
    values = pdf[name].to_numpy(dtype=float)
    valid = np.isfinite(timestamp) & np.isfinite(values) & (values != 0)
    if valid.sum() == 0:
        return np.full(len(pdf), np.nan)
    if valid.sum() == 1:
        return np.full(len(pdf), values[valid][0])
    return interp1d(
        timestamp[valid],
        values[valid],
        bounds_error=False,
        fill_value="extrapolate",
    )(timestamp)

@tools.pass_axes()
def energy(l2, scale='log', E_bin=None, f=None, ax=None, **kwargs):
    if scale not in ('log', 'linear'):
        raise ValueError(f"Unsupported scale: {scale!r}")
    if E_bin is None:
        E_bin = np.logspace(1, 5, 128) / 1000 if scale == 'log' else np.linspace(0, 100, 128)
    ax.hist(l2['Energy'], E_bin, **kwargs)
    ax.set_xlabel('E (keV)')
    ax.set_ylabel('N')
    if scale == 'log':
        ax.set_xscale('log')
        ax.set_xlim(0.5, 100)
    else:
        ax.set_xlim(0, 100)
    return f, ax


@tools.pass_axes()
def layer(l2_df, layer=1, xyrange=None, f=None, ax=None, **kwargs):
    if xyrange is None:
        xyrange = [np.arange(-100, 100, 1), np.arange(-50, 50, 1)]
    if 1 == layer:
        im = ax.hist2d(l2_df['X1'], l2_df['Y1'], xyrange, **kwargs)
    if 2 == layer:
        im = ax.hist2d(l2_df['X2'], l2_df['Y2'], xyrange, **kwargs)
    tools.draw_rectangle(ax, 90, 30)
    tools.draw_rectangle(ax, 150, 60)
    f.colorbar(im[3], ax=ax)
    ax.set_xlabel('X mm')
    ax.set_ylabel('Y mm')
    ax.set_aspect('equal')
    f.autofmt_xdate()
    return f, ax


@tools.pass_axes()
def fov(l2, xrange=None, yrange=None, fxy=(59.86, 19.5), sxy=(45.06, 9.85), f=None, ax=None, **kwargs):
    if xrange is None:
        xrange = np.arange(-35, 35, 2.5)
    if yrange is None:
        yrange = np.arange(-20, 20, 2.5)
    h = ax.hist2d(l2['PX'], l2['PY'], [xrange, yrange], **kwargs)
    tools.draw_rectangle(ax, *fxy)
    tools.draw_rectangle(ax, *sxy)
    ax.set_aspect('equal')
    f.colorbar(h[3], ax=ax)
    ax.set_xlabel('AX deg')
    ax.set_ylabel('AY deg')
    f.autofmt_xdate()
    return f, ax


@tools.pass_axes(subplot_spec=(2, 1, True, False))
def ho_ratio(l2, window='5min', f=None, ax=None,**kwargs):
    l2_H = l2[l2['Spec'] == 1]
    l2_O = l2[l2['Spec'] == 2]
    tH_n, nH_flux, _ = l2_H.window_hist0d(window=window, **kwargs)
    tO_n, nO_flux, _ = l2_O.window_hist0d(window=window, **kwargs)

    ax[0].plot(tH_n, nH_flux, label='H')
    ax[0].plot(tO_n, nO_flux, label='O')
    ax[0].legend()
    species = pd.concat([pd.Series(nH_flux, index=tH_n), pd.Series(nO_flux, index=tO_n)],
                        axis=1).sort_index().fillna(0)
    total = species.sum(axis=1)
    ax[1].plot(species.index, species[1] / total.where(total > 0))
    ax[1].set_ylabel('O/Tot')
    f.autofmt_xdate()
    return f, ax


@tools.pass_axes(subplot_spec=(3, 1, True, False))
def ep_tof(l2, f=None, ax=None, **kwargs):
    idx_H = l2['Spec'] == 1
    idx_O = l2['Spec'] == 2
    print(f"num of H atom = {sum(idx_H)}")
    print(f"num of O atom = {sum(idx_O)}")
    EP_range = np.arange(200, 900, 10)
    Tof_range = np.logspace(-3, 0, 50)
    ax[0].hist2d(l2['EP'], 1 / l2['Tof'], [EP_range, Tof_range], **kwargs)
    ax[1].hist2d(l2['EP'][idx_H], 1 / l2['Tof'][idx_H], [EP_range, Tof_range], **kwargs)
    ax[2].hist2d(l2['EP'][idx_O], 1 / l2['Tof'][idx_O], [EP_range, Tof_range], **kwargs)
    f.autofmt_xdate()
    return f, ax


@tools.pass_axes(subplot_spec=(3, 1, True, False))
def ep_energy_per_nucleon(l2, f=None, ax=None, **kwargs):
    idx_H = l2['Spec'] == 1
    idx_O = l2['Spec'] == 2
    print(f"num of H atom = {sum(idx_H)}")
    print(f"num of O atom = {sum(idx_O)}")
    A = np.ones(len(l2))
    A[idx_O] = 16
    EP_range = np.arange(200, 900, 10)
    E_range = np.arange(1, 100, 2)
    ax[0].hist2d(l2['EP'], l2['Energy'] / A, [EP_range, E_range], **kwargs)
    ax[1].hist2d(l2['EP'][idx_H], l2['Energy'][idx_H], [EP_range, E_range], **kwargs)
    ax[2].hist2d(l2['EP'][idx_O], l2['Energy'][idx_O] / 16, [EP_range, E_range], **kwargs)
    for _ in ax:
        _.set_xscale('log')
        _.set_yscale('log')
        _.set_xlabel('Pulse (mV)')
        _.set_ylabel('E (keV)')
    f.autofmt_xdate()
    return f, ax

def _meshvalues(df, vnames, vb, win, f=None, ax=None):
    all_hist_matrix = []
    for vname in vnames:
        times, bin_centers, hist, _ = df.window_hist1d(vname, vb, window=win)
        all_hist_matrix.append(hist)
    vmax = max(2.0, float(np.max([np.max(hist) for hist in all_hist_matrix])))
    vmin = min(10.0, max(1.0, vmax / 10.0))
    for _ax, hist_matrix in zip(ax, all_hist_matrix):
        mesh = _ax.pcolormesh(times, bin_centers, hist_matrix.T, cmap='viridis', shading='nearest', norm=LogNorm(vmin=vmin, vmax=vmax))
    tools.add_bar_top(ax[0], mesh)


@tools.pass_axes()
def meshvalues(df, vnames, vb, win, f=None, ax=None):
    _meshvalues(df, vnames, vb, win, f=f, ax=ax)
    return f, ax


@tools.pass_axes(subplot_spec=(4, 1, True, False))
def start_signal(l1_df, *, f=None, ax=None, win='5min'):
    vnames = ['TX_1s', 'TX_2s', 'TY_1s', 'TY_2s']
    _meshvalues(l1_df, vnames, np.arange(1, 100, 2), win, f=f, ax=ax)
    for i, _ in enumerate(ax):
        _.set_ylabel(vnames[i])
    for i in range(2):
        ax[i].axhline(50, color='r')
        ax[i + 2].axhline(40, color='r')
    f.autofmt_xdate()
    return f, ax


@tools.pass_axes(subplot_spec=(4, 1, True, False))
def end_signal(l1_df, f=None, ax=None, win='5min'):
    vnames = ['TX_1e', 'TX_2e', 'TY_1e', 'TY_2e']
    _meshvalues(l1_df, vnames, np.arange(1, 700, 2), win, f=f, ax=ax)
    for i, _ in enumerate(ax):
        _.set_ylabel(vnames[i])
    f.autofmt_xdate()
    return f, ax


@tools.pass_axes(subplot_spec=(2, 1, True, False))
def channels(l1a, f=None, ax=None, *, window="1min"):
    for axis, edge in zip(ax, "se"):
        for channel in ("TX_1", "TX_2", "TY_1", "TY_2"):
            name = channel + edge
            valid = l1a[name].between(5, 633, inclusive="neither")
            t, h, _ = sts.window_hist0d(l1a[valid], window=window)
            axis.plot(t, h, label=name)
        axis.set_yscale('log')
        axis.legend()
    f.autofmt_xdate()
    return f, ax


@tools.pass_axes(subplot_spec=(4, 1, True, False))
def parameters(pdf, flag, show_temp=True, f=None, ax=None, **kwargs):
    t_date = pdf['Datetime']
    if flag == 'a':
        ax[0].plot(t_date, abs(_interpolate_nonzero_parameter(pdf, 'HV_A')), **kwargs)
        ax[1].plot(t_date, _interpolate_nonzero_parameter(pdf, 'Th_A_MCP'), **kwargs)
        ax[2].plot(t_date, _interpolate_nonzero_parameter(pdf, 'Th_A_X'), label='AX', **kwargs)
        ax[2].plot(t_date, _interpolate_nonzero_parameter(pdf, 'Th_A_Y'), label='AY', **kwargs)
        if show_temp:
            ax[3].plot(t_date, _interpolate_nonzero_parameter(pdf, 'Temp_A'), **kwargs)
        ax[2].legend(frameon=False, ncol=3)
    elif flag == 'b':
        ax[0].plot(t_date, abs(_interpolate_nonzero_parameter(pdf, 'HV_B')), **kwargs)
        ax[1].plot(t_date, _interpolate_nonzero_parameter(pdf, 'Th_B_MCP'), **kwargs)
        b_gain = _interpolate_nonzero_parameter(pdf, 'Th_B_G')
        ax[2].plot(t_date, _interpolate_nonzero_parameter(pdf, 'Th_B_X1') * 2 + b_gain, label='BX1', **kwargs)
        ax[2].plot(t_date, _interpolate_nonzero_parameter(pdf, 'Th_B_Y1') * 2 + b_gain, label='BY1', **kwargs)
        ax[2].plot(t_date, _interpolate_nonzero_parameter(pdf, 'Th_B_X2') * 2 + b_gain, label='BX2', **kwargs)
        ax[2].plot(t_date, _interpolate_nonzero_parameter(pdf, 'Th_B_Y2') * 2 + b_gain, label='BY2', **kwargs)
        if show_temp:
            ax[3].plot(t_date, _interpolate_nonzero_parameter(pdf, 'Temp_B'), **kwargs)
        ax[2].legend(frameon=False, ncol=3)
    else:
        raise ValueError("flag must be 'a' or 'b'")
    ax[0].set_ylabel('Volt')
    ax[0].set_ylim(2600, 3000)
    ax[1].set_ylabel('MCP')
    ax[1].set_ylim(200, 1250)
    ax[2].set_ylabel('Anode')
    ax[2].set_ylim(20, 700)
    if show_temp:
        ax[3].set_ylabel('Temp')
        ax[3].set_ylim(-10, 40)
    ax[-1].yaxis.set_minor_locator(ticker.MultipleLocator(1))
    f.autofmt_xdate()
    return f, ax


@tools.pass_axes()
def counts( dfs, vmax=16 * 62 * 60 * 5, *, window='5min', labels=None, f=None, ax=None,):
    """Plot time-binned counts for one or more event products."""
    dfs = list(dfs)
    if labels is None:
        labels = [None] * len(dfs)
    else:
        labels = list(labels)
        if len(labels) != len(dfs):
            raise ValueError('labels must have the same length as dfs')

    plotted = False
    for df, label in zip(dfs, labels):
        if df.empty:
            continue
        t, num, _ = df.window_hist0d(window=window)
        ax.plot(t, num, label=label)
        plotted = True

    if plotted:
        if any(label is not None for label in labels):
            ax.legend(ncol=2)
        if vmax is None:
            ax.set_yscale('symlog', linthresh=1)
        else:
            ax.axhline(vmax, color='r')
            ax.set_ylim(1, vmax)
            ax.set_yscale('log')
    else:
        ax.text( 0.5, 0.5, 'No event products available', ha='center', va='center', transform=ax.transAxes,)
    ax.set_ylabel(f'Counts / {window}')
    f.autofmt_xdate()
    return f, ax


@tools.pass_axes(subplot_spec=(4, 3, False, False))
def histogram(l1_df, df_tot=None, uuid=None, pid=None, params=init_param(1), prob='a', *, f=None, ax=None):
    p = prob.upper()
    centers = [params[f'{p}_X_SUM'], params[f'{p}_Y_SUM']]
    rights = [params[f'{p}_XR'], params[f'{p}_YR']]
    lefts = [params[f'{p}_XL'], params[f'{p}_YL']]
    ax = np.array(ax).reshape(4, 3)
    for j, e in enumerate(['X', 'Y']):
        for i, c in enumerate(['s', 'e']):
            i = i + j * 2
            varname1 = f'T{e}_1{c}'
            varname2 = f'T{e}_2{c}'
            ax[i, 0].hist(l1_df[varname1], np.arange(0, int(8192 * 0.081) + 1), label=varname1, alpha=0.5)
            ax[i, 0].hist(l1_df[varname2], np.arange(0, int(8192 * 0.081) + 1), label=varname2, alpha=0.5)
            ax[i, 0].axvline(centers[j], color='k')
            ax[i, 1].hist(l1_df[varname1] + l1_df[varname2], np.arange(0, int(8192 * 0.081 * 2) + 1), label=varname1 + '+' + varname2, alpha=0.5)
            ax[i, 1].axvline(centers[j], color='k')
            ax[i, 2].hist(l1_df[varname1] - l1_df[varname2], np.arange(-int(8192 * 0.081), int(8192 * 0.081) + 1), label=varname1 + '-' + varname2, alpha=0.5)
            ax[i, 2].axvline(lefts[j], color='k')
            ax[i, 2].axvline(rights[j], color='k')
            for k in range(3):
                ax[i, k].legend()
            if df_tot is not None and uuid is not None and pid is not None:
                selected = df_tot.query('UUID == @uuid and PID == @pid')
                if not selected.empty:
                    value1 = selected[varname1].values[0]
                    value2 = selected[varname2].values[0]
                    print(f'{varname1}:{value1}')
                    print(f'{varname2}:{value2}')
                    ax[i, 0].axvline(value1, color='r', ls='--')
                    ax[i, 0].axvline(value2, color='g', ls='--')
                    ax[i, 1].axvline(value1 + value2, color='k', ls='--')
                    ax[i, 2].axvline(value1 - value2, color='k', ls='--')
    return f, ax


@tools.pass_axes(subplot_spec=(2, 2, True, True))
def joint(df, *, f=None, ax=None):
    idx = 0
    for d in ['X', 'Y']:
        for c in ['s', 'e']:
            xvarname = f'T{d}_1{c}'
            yvarname = f'T{d}_2{c}'
            H = ax[idx].hist2d(df[xvarname], df[yvarname], np.arange(0, int(8192 * 0.081) + 1), norm='log')
            ax[idx].set_xlabel(xvarname)
            ax[idx].set_ylabel(yvarname)
            colorbar = f.colorbar(H[-1], ax=ax[idx], orientation='horizontal', pad=0.1, location='top')
            maximum = int(np.nanmax(H[0])) if H[0].size else 0
            if maximum > 0:
                tick_count = min(5, maximum)
                ticks = np.unique(
                    np.rint(np.geomspace(1, maximum, num=tick_count)).astype(int)
                )
                colorbar.set_ticks(ticks)
            colorbar.ax.xaxis.set_minor_locator(ticker.NullLocator())
            colorbar.ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:g}'))
            idx += 1
    f.autofmt_xdate()
    f.subplots_adjust(wspace=0.35, hspace=0.55)
    return f, ax


@tools.pass_axes(subplot_spec=(2, 2, True, True))
def scatter_joint(df, row_id, f=None, ax=None):
    idx = 0
    for d in ['X', 'Y']:
        for c in ['s', 'e']:
            xvarname = f'T{d}_1{c}'
            yvarname = f'T{d}_2{c}'
            ax[idx].axvline(df.iloc[row_id][xvarname], color='r', ls='--')
            ax[idx].axhline(df.iloc[row_id][yvarname], color='r', ls='--')
            ax[idx].set_xlabel(xvarname)
            ax[idx].set_ylabel(yvarname)
            idx += 1
    f.autofmt_xdate()
    f.subplots_adjust(wspace=0.35, hspace=0.35)
    return f, ax


def _add_value_labels(ax, labels, values, bottoms=None):
    for label, value, bottom in zip(labels, values, bottoms or [0] * len(values)):
        bar = ax.bar(label, value, bottom=bottom)[0]
        if value > 0:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., bottom + height, f'{value}', ha='center', va='center', color='black', fontsize=8)

@tools.pass_axes()
def filter(l1, f=None, ax=None, **kwargs):
    xlabel = ['L0', 'triggered', 'start_stop', 'start_ns', 'golden_match']

    tg = flt.triggered_eff(l1)
    ss = flt.start_stop_eff(l1)
    so = flt.start_only_eff(l1)

    # TDC start values within ns bounds (replaces old start_ns_eff)
    ns_check = flt.check('length', dir='all', edge='all', channel='s')(l1)
    ns_within = ss & ns_check

    # Golden standard — all checks pass (replaces old match_eff)
    golden = flt.time_tof_eff(dir='all')(l1)
    golden_match = ss & golden

    _add_value_labels(ax, [xlabel[0]], [len(l1)])
    _add_value_labels(ax, xlabel[1:], [sum(m) for m in [tg, ss, ns_within, golden_match]])

    # Breakdown of "start_stop" bar: start_only on bottom, rest on top
    _add_value_labels(ax, [xlabel[2]], [sum(so)], [0])
    _add_value_labels(ax, [xlabel[2]], [sum(ss) - sum(so)], [sum(so)])

    ax.set_yscale('log')
    return f, ax


def _distribution_bins(prob, ssid, kind):
    params = init_param(ssid)
    if prob == 'a':
        xlim = abs(params['A_XR']) + abs(params['A_XL'])
        ylim = abs(params['A_YR']) + abs(params['A_YL'])
    elif prob == 'b':
        xlim = abs(params['B_XR']) + abs(params['B_XL'])
        ylim = abs(params['B_YR']) + abs(params['B_YL'])
    else:
        print('wrong id')
        return None
    if kind == '8':
        ybins = [1, ylim, ylim + 40, ylim + 80, ylim + 120, ylim + 160, ylim + 200, ylim + 240, ylim + 280, ylim + 320, 9999]
        xbins = [1, xlim, xlim + 50, xlim + 100, xlim + 150, xlim + 200, xlim + 250, xlim + 300, xlim + 350, xlim + 400, 9999]
    else:
        ybins = [1, ylim, ylim + 30, ylim + 60, ylim + 90, ylim + 120, 300]
        xbins = [1, xlim, xlim + 50, xlim + 100, xlim + 150, xlim + 200, xlim + 250, xlim + 300, 700]
    xlabels = [str(int(x1)) + '-' + str(int(x2)) for x1, x2 in zip(xbins[:-1], xbins[1:])]
    ylabels = [str(int(y1)) + '-' + str(int(y2)) for y1, y2 in zip(ybins[:-1], ybins[1:])]
    return xbins, xlabels, ybins, ylabels


@tools.pass_axes()
def distribution(df, prob, ssid=1, channels=None, include_zero=True, kind='4', f=None, ax=None):
    bins = _distribution_bins(prob, ssid, kind)
    if bins is None:
        return f, ax
    xbins, xlabels, ybins, ylabels = bins
    if channels is None:
        channels = ['TX_1s', 'TX_2s', 'TX_1e', 'TX_2e', 'TY_1s', 'TY_2s', 'TY_1e', 'TY_2e'] if kind == '8' else ['TX_1s', 'TX_2s', 'TY_1s', 'TY_2s']
    axes = ax if isinstance(ax, list) else [ax]
    for idx, cnl in enumerate(channels):
        target = axes[idx] if idx < len(axes) else axes[-1]
        series = df[cnl]
        label_bins = xlabels if cnl.startswith('TX') else ylabels
        value_bins = xbins if cnl.startswith('TX') else ybins
        values = series[series > 0] if include_zero else series
        interval_counts = pd.cut(values, bins=value_bins, labels=label_bins, right=False).value_counts().reindex(label_bins, fill_value=0)
        all_counts = pd.concat([pd.Series({'0': (series == 0).sum()}), interval_counts]) if include_zero else interval_counts
        total_count = all_counts.sum()
        all_percentages = all_counts / total_count * 100 if total_count else all_counts * 0
        bars = target.barh(all_counts.index, all_counts.values)
        ymin, ymax = target.get_ylim()
        y_range = ymax - ymin
        for bar, (index, value) in zip(bars, all_counts.items()):
            y_coord = bar.get_y() + bar.get_height() / 2
            normalized_y = (y_coord - ymin) / y_range if y_range else 0
            target.text(0.3, normalized_y, f'{value}', ha='left', va='center', transform=target.transAxes)
            target.text(0.5, normalized_y, f'{all_percentages[index]:.2f}%', ha='left', va='center', transform=target.transAxes)
        target.set_ylabel('Value Intervals')
        target.set_xlabel('Count')
        target.set_title(cnl)
    fig_w = 3.375
    f.set_size_inches(fig_w * 4, fig_w if kind == '8' else fig_w / 4 * 3)
    return f, ax


@tools.pass_axes(subplot_spec=(2, 1, False, False))
def match_distribution(df, folder='./', gap=10, f=None, ax=None):
    sumx = df['TX_1s'] + df['TX_2s']
    sumy = df['TY_1s'] + df['TY_2s']
    bins = np.arange(0, 600 + gap, gap)
    ratiox, barx = np.histogram(sumx, bins, density=True)
    ax[0].bar(barx[:-1], ratiox, gap)
    ax[0].axvline(45, color='r')
    ax[0].axvline(55, color='r')
    ax[0].set_ylabel('Sum X dist')
    ax[0].set_xlim(0, 600)
    ratioy, bary = np.histogram(sumy, bins, density=True)
    ax[1].axvline(35, color='r')
    ax[1].axvline(45, color='r')
    ax[1].bar(bary[:-1], ratioy, gap)
    ax[1].set_xlabel('Sum (ns)')
    ax[1].set_ylabel('Sum Y dist')
    ax[1].set_xlim(0, 600)
    fig_w = 3.375
    f.set_size_inches(fig_w * 4, fig_w * 2)
    return f, ax


@tools.pass_axes(subplot_spec=(1, 2, True, True))
def sum_xy(df, params=init_param(1), prob='a', f=None, ax=None):
    CX = params['A_X_SUM'] if prob == 'a' else params['B_X_SUM']
    CY = params['A_Y_SUM'] if prob == 'a' else params['B_Y_SUM']
    for idx, c in enumerate(['s', 'e']):
        sum_xs = (df[f'TX_1{c}'] + df[f'TX_2{c}'] - CX) / 2
        sum_ys = (df[f'TY_1{c}'] + df[f'TY_2{c}'] - CY) / 2
        ax[idx].hist2d(sum_xs, sum_ys, np.arange(0, int(8192 * 0.081) + 1), norm='log')
        ax[idx].set_xlabel(f'Sum X_{c}')
        ax[idx].set_ylabel(f'Sum Y_{c}')
        ax[idx].set_aspect('equal')
    f.autofmt_xdate()
    f.subplots_adjust(wspace=0.35, hspace=0.35)
    return f, ax


@tools.pass_axes(subplot_spec=(1, 2, True, True))
def diff_xy(df, params=init_param(1), prob='a', f=None, ax=None):
    XL, XR = (params['A_XL'], params['A_XR']) if prob == 'a' else (params['B_XL'], params['B_XR'])
    YL, YR = (params['A_YL'], params['A_YR']) if prob == 'a' else (params['B_YL'], params['B_YR'])
    for idx, c in enumerate(['s', 'e']):
        diff_xs = df[f'TX_1{c}'] - df[f'TX_2{c}']
        diff_ys = df[f'TY_1{c}'] - df[f'TY_2{c}']
        Xbound = np.round(np.max([abs(XR), abs(XL)])).astype(int)
        Ybound = np.round(np.max([abs(YR), abs(YL)])).astype(int)
        ax[idx].hist2d(diff_xs, diff_ys, [np.arange(-Xbound - 10, Xbound + 10), np.arange(-Ybound - 10, Ybound + 10)], norm='log')
        ax[idx].axvline(XL, color='r')
        ax[idx].axvline(XR, color='r')
        ax[idx].axhline(YL, color='r')
        ax[idx].axhline(YR, color='r')
        ax[idx].set_xlabel(f'Diff X_{c}')
        ax[idx].set_ylabel(f'Diff Y_{c}')
        ax[idx].set_aspect('equal')
    f.autofmt_xdate()
    f.subplots_adjust(wspace=0.35, hspace=0.35)
    return f, ax


@tools.pass_axes(subplot_spec=(2, 1, False, False))
def start_end_level2(df, msize=1, s_or_hist=1, f=None, ax=None):
    ax1, ax2 = ax
    if s_or_hist == 0:
        ax1.scatter(df['X1'], df['Y1'], label='start', s=msize)
        ax2.scatter(df['X2'], df['Y2'], label='end', s=msize)
        ax1.legend()
        ax2.legend()
    elif s_or_hist == 1:
        hh = ax1.hist2d(df['X1'], df['Y1'], np.arange(-100, 100, 10))
        f.colorbar(hh[3], ax=ax1)
        hh = ax2.hist2d(df['X2'], df['Y2'], np.arange(-100, 100, 10))
        f.colorbar(hh[3], ax=ax2)
    elif s_or_hist == 2:
        xy = np.vstack([df['X1'], df['Y1']])
        z = 1 - gaussian_kde(xy)(xy)
        ax1.scatter(df['X1'], df['Y1'], c=z, s=1, cmap='Spectral')
        xy = np.vstack([df['X2'], df['Y2']])
        z = 1 - gaussian_kde(xy)(xy)
        ax2.scatter(df['X2'], df['Y2'], c=z, s=1, cmap='Spectral')
    tools.draw_rectangle(ax1, 90, 30)
    tools.draw_rectangle(ax2, 90, 30)
    tools.draw_rectangle(ax1, 150, 60)
    tools.draw_rectangle(ax2, 150, 60)
    ax1.set_xlabel('X (mm)')
    ax1.set_ylabel('Y (mm)')
    ax2.set_xlabel('X (mm)')
    ax2.set_ylabel('Y (mm)')
    return f, ax


@tools.pass_axes()
def tdc_bar(df, cnl, channel=None, f=None, ax=None):
    series = df[cnl]
    zero_count = (series == -94 * 0.081).sum()
    if channel == 'x':
        bins = [1, 100, 200, 300, 400, 500, 600]
        labels = ['1 - 100', '100 - 200', '200 - 300', '300 - 400', '400 - 500', '500 - 600']
    elif channel == 'y':
        bins = [1, 50, 100, 150, 200, 250, 300]
        labels = ['1 - 50', '5₀ - 1₀₀', '1₀₀ - 1₅₀', '1₅₀ - 2₀₀', '2₀₀ - 2₅₀', '2₅₀ - 3₀₀']
    else:
        bins = [1, 100, 200, 300, 400, 500, 600, 1000]
        labels = ['1 - 100', '100 - 200', '200 - 300', '300 - 400', '400 - 500', '500 - 600', '>600']
    interval_counts = pd.cut(series[series > 0], bins=bins, labels=labels, right=False).value_counts().reindex(labels, fill_value=0)
    all_counts = pd.concat([pd.Series({'0': zero_count}), interval_counts])
    if channel is None:
        all_percentages = all_counts / all_counts.sum() * 100
        ax.pie(all_percentages, labels=all_percentages.index, autopct='%1.1f%%')
    else:
        total_count = all_counts.sum()
        all_percentages = all_counts / total_count * 100 if total_count else all_counts * 0
        bars = ax.bar(all_counts.index, all_counts.values)
        for bar, (index, value), locx in zip(bars, all_counts.items(), np.linspace(0, 1, len(all_counts))):
            ax.text(locx, 0.2, f'{value}', transform=ax.transAxes, ha='center')
            ax.text(locx, 0.1, f'{all_percentages[index]:.5f}%', transform=ax.transAxes, ha='center')
        ax.set_xlabel('Value Intervals')
        ax.set_ylabel('Count')
        ax.set_yticks([])
        ax.set_title(cnl)
        fig_w = 3.375
        f.set_size_inches(fig_w * 4, fig_w / 4 * 3)
    return f, ax
