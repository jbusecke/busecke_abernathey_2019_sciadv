# Functions used in busecke and abernathey 2018
import xarray as xr
import numpy as np
from astropy.convolution import (convolve_fft, Gaussian1DKernel)
import dask.array as dsa
from scipy import stats
import matplotlib.pyplot as plt

def calculate_Kmix(eke, u, v, u_zt, v_zt, L_D, L_obs, Gamma=0.35):
    """Calculate Diffusivity following the Supressed Mixing Length Theory (SMLT, Ferrari and Nikurashin 2010)
    Parameters
    ----------
    eke : xr.DataArray
        Eddy Kinetic Energy in m^2/s^2.
        
    u : xr.DataArray
        Large scale zonal velocity (e.g. low pass filtered) in m/s
        
    v : xr.DataArray
        Large scale meridional velocity (e.g. low pass filtered) in m/s
        
    u_zt : xr.DataArray
        Depth and time integrated zonal velocity in m/s
    
    v_zt : xr.DataArray
        Depth and time integrated meridional velocity in m/s
    
    L_D : xr.DataArray
        First Baroclinic Rossby radius of deformation in m
        
    L_obs : xr.DataArray
        Observed Eddy Lengthscale (diameter) in m
    

    Returns
    -------
    ds : xr.Dataset
        Dataset with all calculated parameters. For detailed description see (paper link)
    
    """
    
    ds = xr.Dataset()
    
    ds['u_rms'] = xr.ufuncs.sqrt(2 * eke)
    ds['u_bar'] = u
    ds['v_bar'] = v
    
    ds['delta_u_rms'] = ds['u_rms'] - ds['u_rms'].mean('time')
    ds['delta_u_bar'] =  ds['u_bar'] -  ds['u_bar'].mean('time')
    ds['delta_v_bar'] = ds['v_bar'] - ds['v_bar'].mean('time')
    
    ds['L_D'] = L_D
    
    ds['L_obs'] = L_obs
    

    ds['k'] = 2*np.pi/ds['L_obs']
    
    ds['gamma'] = 1/(4*24*60*60)
    ds['alpha'] = (ds['k']**2) / (ds['gamma']**2)

    omega=2*np.pi/(24*60*60) # Earths rotational frequency
    R = 6.371*1e6 # Earths Radius
    ref_da = u.mean(['time', 'XC'])
    ds['beta'] = xr.DataArray(2.0*omega*np.cos(np.deg2rad(ds.YC))/R,
                            dims=ref_da.dims,
                            coords=ref_da.coords,
                            name='beta')
    
    # Mixing Efficienct=y 
    ds['Gamma'] = Gamma
    
    ds['K_0'] = ds['Gamma'] * ds['u_rms'] * ds['L_obs']
    ds['delta_K_0'] = ds['Gamma'] * ds['L_obs'] * ds['delta_u_rms']
    
    # phase speed
    ds['theoretical_c'] = - ds['beta'] * (ds['L_D']**2)

    ds['c_x'] = u_zt + ds['theoretical_c']
    ds['c_y'] = v_zt

    ds['speed_diff_x'] = (ds['c_x'] - ds['u_bar'])
    ds['speed_diff_y'] = (ds['c_y'] - ds['v_bar'])
    
    ds['S_x'] = 1 / (1 + (ds['alpha'] * ds['speed_diff_x']**2))
    ds['S_y'] = 1 / (1 + (ds['alpha'] * ds['speed_diff_y']**2))
    
    
    
    ds['dS_x_dU'] = 2 * ds['alpha'] * ds['speed_diff_x'] * (ds['S_x']**2)
    ds['dS_y_dV'] = 2 * ds['alpha'] * ds['speed_diff_y'] * (ds['S_y']**2)
    
    ds['delta_S_x'] = ds['dS_x_dU'] * ds['delta_u_bar']
    ds['delta_S_y'] = ds['dS_y_dV'] * ds['delta_v_bar']
    
    # just the zonal components
    ds['S'] = ds['S_x']
    ds['delta_S'] = ds['delta_S_x']
    
    ds['eke_variation'] = ds['S'].mean('time') * ds['delta_K_0']
    ds['ml_variation'] = ds['K_0'].mean('time') * ds['delta_S']
    ds['delta_K_mix'] = ds['eke_variation'] + ds['ml_variation']
    ds['mean_K_mix'] = ds['Gamma'] * ds['L_obs'] * ds['S'].mean('time') * ds['u_rms'].mean('time') 
    ds['K_mix'] = ds['mean_K_mix'] + ds['delta_K_mix']
    return ds


def annual_range(da, time_spec = '1A', range_cut = [0.15, 0.85]):
    """Calculate the annual range of a dataarray"""
    annual = da.resample(time='1A').mean('time')
    annual_range = annual.load().quantile(range_cut[1], 'time') - annual.load().quantile(range_cut[0], 'time')
    return annual_range

def xarray_norm(in_a, norm_a):
    """normalize and scale `in_a` with mean and std of `norm_a`"""
    out_norm = ((in_a - np.nanmean(in_a.data)) / np.nanstd(in_a.data) *
                np.nanstd(norm_a)) + np.nanmedian(norm_a)
    return out_norm


##################################################################################
# The following functions originate from https://github.com/jbusecke/xarrayutils,# 
# but to ensure reproducibility and reduce dependencies, they are copied below   #
##################################################################################
def dict2box(di, xdim='lon', ydim='lat'):
    return np.array([di[xdim].start, di[xdim].stop,
                     di[ydim].start, di[ydim].stop])


def box_plot(box, ax=None, split_detection='True', **kwargs):
    """plots box despite coordinate discontinuities.
    INPUT
    -----
    box: np.array
        Defines the box in the coordinates of the current axis.
        Describing the box corners [x1, x2, y1, y2]
    ax: matplotlib.axis
        axis for plotting. Defaults to plt.gca()
    kwargs: optional
        anything that can be passed to plot can be put as kwarg
    """

    if len(box) != 4:
        raise RuntimeError("'box' must be a 4 element np.array, \
            describing the box corners [x1, x2, y1, y2]")
    xlim = plt.gca().get_xlim()
    ylim = plt.gca().get_ylim()
    x_split = False
    y_split = False

    if ax is None:
        ax = plt.gca()

    if split_detection:
        if np.diff([box[0], box[1]]) < 0:
            x_split = True

        if np.diff([box[2], box[3]]) < 0:
            y_split = True

    if y_split and not x_split:
        ax.plot([box[0], box[0], box[1], box[1], box[0]],
                 [ylim[1], box[2], box[2], ylim[1], ylim[1]], **kwargs)

        ax.plot([box[0], box[0], box[1], box[1], box[0]],
                 [ylim[0], box[3], box[3], ylim[0], ylim[0]], **kwargs)

    elif x_split and not y_split:
        ax.plot([xlim[1], box[0], box[0], xlim[1], xlim[1]],
                 [box[2], box[2], box[3], box[3], box[2]], **kwargs)

        ax.plot([xlim[0], box[1], box[1], xlim[0], xlim[0]],
                 [box[2], box[2], box[3], box[3], box[2]], **kwargs)

    elif x_split and y_split:
        ax.plot([xlim[1], box[0], box[0]], [box[2], box[2], ylim[1]],
                 **kwargs)

        ax.plot([xlim[0], box[1], box[1]], [box[2], box[2], ylim[1]],
                 **kwargs)

        ax.plot([xlim[1], box[0], box[0]], [box[3], box[3], ylim[0]],
                 **kwargs)

        ax.plot([xlim[0], box[1], box[1]], [box[3], box[3], ylim[0]],
                 **kwargs)

    elif not x_split and not y_split:
        ax.plot([box[0], box[0], box[1], box[1], box[0]],
                 [box[2], box[3], box[3], box[2], box[2]], **kwargs)


        
def box_plot_dict(di, xdim='lon', ydim='lat', **kwargs):
    """plot box from xarray selection dict e.g.
    `{'xdim':slice(a, b), 'ydim':slice(c,d), ...}`"""

    # extract box from dict
    box  = dict2box(di, xdim=xdim, ydim=ydim)
    # plot
    box_plot(box, **kwargs)
    
        
def filter_1D(data, std, dim='time', dtype=None):
    if dtype is None:
        dtype = set_dtype(data)

    kernel = Gaussian1DKernel(std)

    def smooth_raw(data):
        raw_data = getattr(data, 'values', data)
        result = convolve_fft(raw_data, kernel, boundary='wrap')
        result[np.isnan(raw_data)] = np.nan
        return result

    def temporal_smoother(data):
        dims = ([dim])
        return xr.apply_ufunc(smooth_raw, data,
                              vectorize=True,
                              dask='parallelized',
                              input_core_dims=[dims],
                              output_core_dims=[dims],
                              output_dtypes=[dtype])
    return temporal_smoother(data)

def set_dtype(aa):
    if isinstance(aa, xr.Dataset):
        dtype = aa[list(aa.data_vars)[0]].dtype
        print('No `dtype` chosen. Input is Dataset. \
        Defaults to %s' % dtype)
    elif isinstance(aa, xr.DataArray):
        dtype = aa.dtype
    return dtype


def weighted_mean(da_data, da_weight, **kwargs):
    """calculate average of da_data weighted by da_weight
    Parameters
    ----------
    da_data : xarray.DataArray
        Data to be averaged
    da_weight : xarray.DataArray
        weights to be used during averaging. Dimensions have to be
        matching with 'da_data'
    dim : {None, str, list}, optional
        Dimensions to average over
    preweighted: Bool, optional
        Specifies whether weights will be applied (False, default) or
        have already been
        applied to da_data (True).
    dim_check: Bool, optional
        Activates check for dimension consistency. If dimensions of 'da_weight'
        do not include all elements of 'dim' error is raised
    """
    data, weight_expanded = weighted_sum_raw(da_data, da_weight, **kwargs)
    out = data / weight_expanded
    out.attrs = data.attrs
    return out


def weighted_sum_raw(da_data, da_weight, dim=None,
                     preweighted=False, dimcheck=True, **kwargs):
    """calculate sum of da_data weighted by da_weight and the weights themselves
    Parameters
    ----------
    da_data : xarray.DataArray
        Data to be averaged
    da_weight : xarray.DataArray
        weights to be used during averaging. Dimensions have to be matching
        with 'da_data'
    dim : {None, str, list}, optional
        Dimensions to average over
    preweighted: Bool, optional
        Specifies whether weights will be applied (False, default) or have
        already been
        applied to da_data (True).
    dim_check: Bool, optional
        Activates check for dimension consistency. If dimensions of
        'da_weight' do not include all elements of 'dim' error is raised
    """
    if isinstance(dim, str):
        dim = [dim]

    # Check dimension consistency
    if dim:
        if dimcheck:
            if not set(dim).issubset(da_weight.dims):
                raise RuntimeError("Dimensions of 'da_weight' do not include all averaging dimensions.\
                Broadcast da_weight or deactivate 'dim_check'.")
    if 'keep_attrs' in kwargs.keys():
        keep_attrs = kwargs['keep_attrs']
    else:
        keep_attrs = False

    weight_expanded = _broadcast_weights(da_data, da_weight,
                                         keep_attrs=keep_attrs)

    if preweighted:
        data = da_data
    else:
        data = da_data * weight_expanded
        data.attrs = da_data.attrs

    return data.sum(dim, **kwargs), weight_expanded.sum(dim, **kwargs)


def _broadcast_weights(da_data, da_weight, keep_attrs=False):
    """broadcasts da_weights to the same shape as da_data and \
    masks the same missing values"""
    da_data = da_data.copy()
    ones = (da_data * 0) + 1
    weights_expanded = ones * da_weight
    # add attrs back in
    if keep_attrs:
        weights_expanded.attrs = da_data.attrs
    return weights_expanded



def _linregress_ufunc(a, b):
    '''ufunc to wrap scipy.stats.linregress for xr_linregress'''
    slope, intercept, r_value, p_value, std_err = stats.linregress(a, b)
    return np.array([slope, intercept, r_value, p_value, std_err])


def xr_linregress(a, b, dim='time', convert_to_dataset=True, dtype=None):
    """Applies scipy.stats.linregress over two xr.DataArrays or xr.Datasets.
    Parameters
    ----------
    a : {xr.DataArray}
        Independent variable for linear regression. E.g. time.
    b : {xr.DataArray, xr.Dataset}
        Dependent variable.
    dim : str
        Dimension over which to perform linear regression.
        Must be present in both `a` and `b` (the default is 'time').
    convert_to_dataset : bool
        Converts the output parameter to data_variables instead of an
        additional dimension (the default is True).
    dtype : dtype
         Dtype for the output. If None, defaults to dtype of `b`, or `b`s
         first data variable(the default is None).
    Returns
    -------
    type(b)
        Returns a dataarray containing the parameter values of
        scipy.stats.linregress for each data_variable in `b`.
    """

    if dtype is None:
        dtype = set_dtype(b)

    stats = xr.apply_ufunc(_linregress_ufunc, a, b,
                           input_core_dims=[[dim], [dim]],
                           output_core_dims=[['parameter']],
                           vectorize=True,
                           dask='parallelized',
                           output_dtypes=[dtype],
                           output_sizes={'parameter': 5}
                           )
    stats['parameter'] = xr.DataArray(['slope', 'intercept',
                                       'r_value', 'p_value',
                                       'std_err'], dims=['parameter'])
    if convert_to_dataset:
        # chug them all into a dataset
        ds_stats = xr.Dataset()
        for ff in stats.parameter.data:
            ds_stats[ff] = stats.sel(parameter=ff)
        out = ds_stats
    else:
        out = stats
    return out


def linear_trend(obj, dim):
    """Convenience wrapper for 'xr_linregress'. Calculates the trend per
    given timestep. E.g. if the data is passed as yearly values, the
    trend is in units/yr.
    """
    x = xr.DataArray(np.arange(len(obj[dim])).astype(np.float), dims=dim,
                     coords={dim: obj[dim]})
    trend = xr_linregress(x, obj, dim=dim, convert_to_dataset=False)
    return trend




        