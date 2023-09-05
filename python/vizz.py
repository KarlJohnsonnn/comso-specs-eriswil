import numpy as np
import xarray as xr
import datetime
from typing import Dict
from tqdm.auto import tqdm
import copy
import multiprocessing
import concurrent.futures
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patheffects as PathEffects
import matplotlib.cm as cm
import matplotlib.dates as md
import seaborn as sns
import re
from ipywidgets import interact, widgets, fixed

HHLd = np.array(
    [21.500002  , 20.514286  , 19.556965  , 18.62768   , 17.72607   ,
    16.851786  , 16.004465  , 15.183751  , 14.389286  , 13.620715  ,
    12.87768   , 12.159819  , 11.474235  , 10.83433   , 10.231516  ,
    9.650918  ,  9.092202  ,  8.555033  ,  8.039072  ,  7.5439944 ,
    7.069464  ,  6.615141  ,  6.180692  ,  5.76579   ,  5.370097  ,
    4.993273  ,  4.634993  ,  4.294922  ,  3.9727197 ,  3.6680555 ,
    3.3805976 ,  3.110009  ,  2.855953  ,  2.6181023 ,  2.3961203 ,
    2.18967   ,  1.9984195 ,  1.8220367 ,  1.6601849 ,  1.5125275 ,
    1.3787366 ,  1.2584761 ,  1.151409  ,  1.0572032 ,  0.9755266 ,
    0.9060427 ,  0.84841454,  0.8023148 ,  0.76740676,  0.7433537 ]
    )

RGRENZ = np.array(
    [1.25992106e-09, 1.58740110e-09, 1.99999994e-09, 2.51984211e-09,
    3.17480220e-09, 3.99999989e-09, 5.03968423e-09, 6.34960440e-09,
    7.99999977e-09, 1.00793685e-08, 1.26992088e-08, 1.59999995e-08,
    2.01587369e-08, 2.53984176e-08, 3.19999991e-08, 4.03174738e-08,
    5.07968352e-08, 6.39999982e-08, 8.06349476e-08, 1.01593670e-07,
    1.27999996e-07, 1.61269895e-07, 2.03187341e-07, 2.55999993e-07,
    3.22539790e-07, 4.06374681e-07, 5.11999986e-07, 6.45079581e-07,
    8.12749363e-07, 1.02399997e-06, 1.29015916e-06, 1.62549873e-06,
    2.04799994e-06, 2.58031832e-06, 3.25099745e-06, 4.09599988e-06,
    5.16063665e-06, 6.50199490e-06, 8.19199977e-06, 1.03212733e-05,
    1.30039898e-05, 1.63839995e-05, 2.06425466e-05, 2.60079796e-05,
    3.27679991e-05, 4.12850932e-05, 5.20159592e-05, 6.55359981e-05,
    8.25701864e-05, 1.04031918e-04, 1.31071996e-04, 1.65140373e-04,
    2.08063837e-04, 2.62143993e-04, 3.30280745e-04, 4.16127674e-04,
    5.24287985e-04, 6.60561491e-04, 8.32255348e-04, 1.04857597e-03,
    1.32112298e-03, 1.66451070e-03, 2.09715194e-03, 2.64224596e-03,
    3.32902139e-03, 4.19430388e-03]
)


#### COSMO-SPECS 4D output variable metadata

metadata_dict = {
    "variables": [
        {
            "name": "t",
            "units": "K",
            "long_name": "Temperature"
        },
        {
            "name": "qv",
            "units": "kg/kg",
            "long_name": "humidity mixing ratio (MR)"
        },
        {
            "name": "ut",
            "units": "unknown",
            "long_name": "wind component"
        },
        {
            "name": "vt",
            "units": "unknown",
            "long_name": "wind component"
        },
        {
            "name": "wt",
            "units": "unknown",
            "long_name": "wind component"
        },
        {
            "name": "qc",
            "units": "kg/kg",
            "long_name": "cloud liquid water MR"
        },
        {
            "name": "qr",
            "units": "kg/kg",
            "long_name": "rain water MR"
        },
        {
            "name": "qi",
            "units": "kg/kg",
            "long_name": "ice water MR"
        },
        {
            "name": "qs",
            "units": "kg/kg",
            "long_name": "snow water MR"
        },
        {
            "name": "nw",
            "units": "#/kg",
            "long_name": "liquid droplet number concentration"
        },
        {
            "name": "qw",
            "units": "kg/kg",
            "long_name": "water mass in liquid droplets MR"
        },
        {
            "name": "qws",
            "units": "kg/kg",
            "long_name": "soluble aerosol mass in liquid droplets MR"
        },
        {
            "name": "qwa",
            "units": "kg/kg",
            "long_name": "total aerosol mass in liquid droplets MR"
        },
        {
            "name": "nf",
            "units": "#/kg",
            "long_name": "mixed-phase droplet number concentration"
        },
        {
            "name": "qf",
            "units": "kg/kg",
            "long_name": "frozen water mass in mixed-phase droplets MR"
        },
        {
            "name": "qfs",
            "units": "kg/kg",
            "long_name": "soluble aerosol mass in mixed-phase droplets MR"
        },
        {
            "name": "qfa",
            "units": "kg/kg",
            "long_name": "total aerosol mass in mixed-phase droplets MR"
        },
        {
            "name": "qfw",
            "units": "kg/kg",
            "long_name": "liquid water mass in mixed-phase droplets MR"
        },
        {
            "name": "ni",
            "units": "#/kg",
            "long_name": "insoluble aerosol particles number concentration"
        },
        {
            "name": "qia",
            "units": "kg/kg",
            "long_name": "insoluble aerosol mass in insoluble aerosol particles MR"
        }
    ]
}



def format_string(s: str) -> str:
    """Formats input string to a certain format."""
    return "INP = {}, FLARE = {}".format(*re.findall(r'[0-9\.e]+', s), 'FEno' if 'FEno' in s else '') 

def bresenham_line(x0, y0, x1, y1):
    """Bresenham's Line Algorithm to compute indices of a line in a 2D grid."""
    points = []
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    err = dx - dy

    while True:
        points.append((x0, y0))
        if x0 == x1 and y0 == y1:
            break
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x0 += sx
        if e2 < dx:
            err += dx
            y0 += sy

    return np.array(points)[:, 0], np.array(points)[:, 1]

def line_indices(domain_lat, domain_lon, start_lat, start_lon, end_lat, end_lon):
    """Compute the indices in a 2D domain that a straight line crosses."""
    lat_res = (domain_lat[-1] - domain_lat[0]) / len(domain_lat)
    lon_res = (domain_lon[-1] - domain_lon[0]) / len(domain_lon)

    # Convert lat-lon to indices
    x0 = int((start_lon - domain_lon[0]) / lon_res)
    y0 = int((start_lat - domain_lat[0]) / lat_res)
    x1 = int((end_lon - domain_lon[0]) / lon_res)
    y1 = int((end_lat - domain_lat[0]) / lat_res)

    return bresenham_line(x0, y0, x1, y1)


# Haversine function to compute the distance between two lat-lon points
def haversine_distance(lat1, lon1, lat2, lon2):
    R = 6371.0e3  # Earth radius in meters
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = (np.sin(dlat / 2.0)**2 +
         np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0)**2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distance = R * c
    return distance


# Define a function to add metadata
def add_metadata(ds):
    # Get all variables that begin with 'd'
    drop_vars = [var for var in ds.variables if var.startswith('d')]

    # Drop these variables from the dataset
    ds = ds.drop_vars(drop_vars)
    
    for var_metadata in metadata_dict['variables']:
        var_name = var_metadata['name']
        if var_name in ds:
            ds[var_name].attrs['units'] = var_metadata['units']
            ds[var_name].attrs['long_name'] = var_metadata['long_name']
    return ds

def read_file_list(filelist: list[str] , varnames: list[str]):
    data = {}

    # Load each file into a separate xarray Dataset and store in the dictionary
    for file in filelist:# + other_files:
        # Use the file name as the key
        key = file.rsplit('/', 1)[-1][:-3]
    
        # Open the dataset and assign to the dictionary
        data[key] = xr.open_mfdataset(file, preprocess=add_metadata)
        if len(varnames) > 0:
            variables_to_drop = [
                var for var in data[key].variables if var not in varnames
                ]
            data[key] = data[key].drop_vars(variables_to_drop)

    return data


def extract_matching_values(data):
    # Define the regex pattern
    pattern = r'^DIACTL\.\d+\.\d+$'
    
    # Extract values for keys that match the pattern
    matching_values = [value for key, value in data.items() if re.match(pattern, key)]
    
    return matching_values


class MultiPanelPlot:
    def __init__(

        self, 
        datasets: Dict[str, xr.Dataset] | list[str], 
        varname: str, 
        mode: str,         
        metadata: Dict[str, Dict] = None, # type: ignore
        timestep0: int | list = None, # type: ignore
        timeframe: str = 'single', 
        latlon: bool = False,
        title: str = '',
        windvectors: bool = False,
        show: bool = False,
        grid: bool = False,
        nrows: int = 2,         ncols: int = 3, 
        vmin: float = 1.0,      vmax: float = 1.0e4, 
        hmin: float = 41,       hmax: float = 44, 
        xmin: float = 1.0e-9,   xmax: float = 1.0e-2, 
        ymin: float = 0.0,      ymax: float = 12.0, 
        idX: int = 12,          idY: int = 12,

    ) -> None:
        """
        Initializer for the MultiPanelPlot class. 

        Args:
            datasets (dict of str: xr.Dataset): The dictionary of data sets for plotting.
            varname (str): The variable name in the data sets.
            mode (str): The plotting mode ('profile' or 'area').
            vmin (float, optional): The minimum value for the color bar. Defaults to 1.0.
            vmax (float, optional): The maximum value for the color bar. Defaults to 1.0e4.
            nrows (int, optional): The number of rows in the subplot. Defaults to 2.
            ncols (int, optional): The number of columns in the subplot. Defaults to 3.
            hmin (float, optional): The minimum model level height. Defaults to 42.
            hmax (float, optional): The maximum model level height. Defaults to 44.
            timeframe (str, optional): The timeframe to plot ('single' or other). Defaults to 'single'.
            xmin (float, optional): The minimum x-value. Defaults to 1.0e-9.
            xmax (float, optional): The maximum x-value. Defaults to 1.0e-2.
            ymin (float, optional): The minimum y-value. Defaults to 0.0.
            ymax (float, optional): The maximum y-value. Defaults to 12.0
            timestep0 (int or list, optional): The time step for the plot. Defaults to 0 or [0, 10].
            metadata (Dict): Metadata from ensemble runs. Defaults None.
            title (str): Suptitle of the figure.
        """
       
        if isinstance(datasets, list):
            self.datasets = read_file_list(datasets, varnames=[varname, 'vt', 'ut', 'wt', 't', 'rho'])

        elif isinstance(datasets, dict):
            self.datasets = copy.deepcopy(datasets)
        else:
            TypeError('dataset can only be list of str or dict of xr.Datasets')
        
        key = next(iter(self.datasets))

        self.varname =  varname
        self.mode = mode
        self.md = metadata
        self.show = show
        self.windvectors = windvectors
        self.grid = grid

        dt0 = datetime.datetime.strptime(self.md[key[3:]]['RUNCTL.ydate_ini'], '%Y%m%d%H')
        self.time = [dt0+datetime.timedelta(seconds=float(15*its)) for its in self.datasets[key].time.values]
        if self.md and latlon:
            latlon00 = extract_matching_values(self.md[key[3:]])
            lat0 = float(latlon00[0]) + float(self.md[key[3:]]['LMGRID.startlat_tot'])
            lon0 = float(latlon00[1]) + float(self.md[key[3:]]['LMGRID.startlon_tot'])
            latN = float(latlon00[0]) - float(self.md[key[3:]]['LMGRID.startlat_tot'])
            lonM = float(latlon00[1]) - float(self.md[key[3:]]['LMGRID.startlon_tot'])
            print(lat0, latN, lon0, lonM)
            self.x = np.linspace(lon0, lonM, self.datasets[key].x.values.size)
            self.y = np.linspace(lat0, latN, self.datasets[key].y.values.size)
        else:
            self.x = self.datasets[key].x.values
            self.y = self.datasets[key].y.values

        self.height = HHLd
        self.xlim = (xmin, xmax)
        self.ylim = (ymin, ymax) 
        self.hlim = (hmin, hmax)
        self.idXY = (idX, idY) 
        self.timeframe = timeframe
        self.bins = np.arange(self.datasets[key][varname].shape[-1])
        self.rgrenz = RGRENZ
        self.unit = self.datasets[key][varname].attrs['units']
        self.cmap = mpl.colors.LinearSegmentedColormap.from_list(
            'mycmap', [
                (0, 'white'), (0.2, 'blue'), (0.4, 'cyan'), (0.5, 'lime'), 
                (0.6, 'yellow'), (0.9, 'red'), (1, 'purple')
                ]
            ) 
        
        self.norm = colors.LogNorm(vmin, vmax)
        self.n_axes = nrows * ncols

        self.pmeshs = []
        self.quivers = []

        my_dpi = 300
        plt.rcParams.update({'font.size': 9})

        self.fig, self.axes = plt.subplots(
            nrows, 
            ncols, 
            figsize=(1920/my_dpi*ncols, 1380/my_dpi*nrows), 
            #figsize=(ncols*4, nrows*3.33),
            dpi=my_dpi, 
            constrained_layout=True)
        if nrows*ncols==1:
            self.axes = np.array(self.axes)
        
        self.timestep = 0 if (timestep0 is None and timeframe == 'single') else (timestep0 if timestep0 is not None else [0, 10])

        self.init_first_plot(self.timestep)
        if title:
            self.fig.suptitle(title, weight='bold')
        


    def set_colorbar(self) -> None:
        """
        Creates and sets a color bar for the figure.
        """
        cbar_ax = self.fig.add_axes([0.91, 0.21, 0.01, 0.6])  # adjust these values as needed
        cbar = self.fig.colorbar(
            cm.ScalarMappable(norm=self.norm, cmap=self.cmap), 
            cax=cbar_ax, 
            orientation='vertical', 
            extend='both', 
            )
        cbar.set_label(self.unit)
        self.colorbar = cbar


    def set_name_tick_params(self, ax: mpl.axes.Axes, name: str) -> None: # type: ignore
        """
        Sets the title and tick parameters for an axes.

        Args:
            ax (mpl.axes.Axes): The axes to configure.
            name (str): The title for the axes.
        """
        #ax.set_title(format_string(name))
        ax.tick_params(which='both', direction='in')
        ax.minorticks_on()
        if self.grid:
            ax.grid(True, which='major', linestyle='--', linewidth='0.11', color='black', alpha=0.5)
        ax.grid(True, which='minor', linestyle=':', linewidth='0.075', color='black', alpha=0.25)


    def init_first_plot(self, timestep: int | list) -> None:
        """
        Initializes the first plot according to the plot mode.

        Args:
            timestep (int): The time step for the plot.
        """
        if self.md is not None:
            for date, ax in zip(self.md.keys(), self.axes.flatten()):
                ax.set_title(f'INP = {self.md[date]["SBM_PAR.dnap_init"]},  FE = {self.md[date]["FLARE_SBM.flare_emission"]}')
                
        if self.mode == 'profile' and isinstance(timestep, int):
            self._init_profile_plot(timestep)
        elif self.mode == 'timeseries' and isinstance(timestep, int):
            self._init_timeseries_plot()         
        elif self.mode == 'area':
            self._init_area_plot(timestep)
        elif self.mode == 'single_spectra':
            self._init_single_spectra_plot(timestep)
        else:
            raise ValueError('Wrong mode given. Available: "profile", "timeseries", "area"')
        
        self.fig.subplots_adjust(hspace=0.3, wspace=0.15, left=0.13, right=0.9, top=0.82, bottom=0.17)  # make space for colorbar
        

    def _init_profile_plot(self, timestep: int) -> None:
        """
        Initializes the profile plot.

        Args:
            timestep (int): The time step for the plot.
        """
        for (name, ds), ax in zip(self.datasets.items(), self.axes.flat):
            data = ds[self.varname].isel(time = timestep, y = 12, x = 12)
            pcm = ax.pcolormesh(self.rgrenz, self.height, data, norm=self.norm, cmap=self.cmap) # type: ignore
            self.pmeshs.append(pcm)
            ax.set(ylim=self.ylim, xlim=self.xlim, xscale='log') # type: ignore
            self.set_name_tick_params(ax, name)
        
        self.fig.text(0.5, 0.06, 'radius [m]', ha='center', va='center')
        self.fig.text(0.1, 0.5, 'height [km]', ha='center', va='center', rotation='vertical')
        self.set_colorbar()      
        self.fig.subplots_adjust(hspace=0.3, wspace=0.15, left=0.13, right=0.9, top=0.82, bottom=0.17)  # make space for colorbar
              


    def _init_timeseries_plot(self) -> None:
        """
        Initializes the timeseries plot.

        Args:
            timestep (int): The time step for the plot.
        """
        for (name, ds), ax in zip(self.datasets.items(), self.axes.flat):
            data = ds[self.varname].isel(y = self.idXY[1], x = self.idXY[0]).sum(dim='bin')
            pcm = ax.pcolormesh(self.time, self.height, data.T, norm=self.norm, cmap=self.cmap) # type: ignore
            self.pmeshs.append(pcm)
            ax.set(ylim=self.ylim) # type: ignore
            self.set_name_tick_params(ax, name)
            ax.xaxis.set_major_locator(md.MinuteLocator(byminute = [0, 30]))
            ax.xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
        

        self.fig.text(0.5, 0.07, 'time [UTC]', ha='center', va='center')
        self.fig.text(0.1, 0.5, 'height [km]', ha='center', va='center', rotation='vertical')
        self.set_colorbar()
        
        self.fig.subplots_adjust(hspace=0.3, wspace=0.15, left=0.13, right=0.9, top=0.82, bottom=0.17)  # make space for colorbar
        


    def _init_area_plot(self, timestep: int | list) -> None:
        """
        Initializes the area plot.

        Args:
            timestep (int): The time step for the plot.
        """
        ds = None
        for (name, ds), ax in zip(self.datasets.items(), self.axes.flat):
            data = self._get_data_area(ds, timestep, self.hlim[0], self.hlim[1])
            pcm = ax.pcolormesh(self.x, self.y, data, norm=self.norm, cmap=self.cmap) # type: ignore
            self.pmeshs.append(pcm)
            self.set_name_tick_params(ax, name)
            if self.windvectors and isinstance(timestep, int):
                self.quivers.append(self.add_windvectors(ds, ax, timestep))
        
        self.fig.text(0.5, 0.09, 'lon', ha='center', va='center')
        self.fig.text(0.1, 0.5, 'lat', ha='center', va='center', rotation='vertical')
        self.set_colorbar()

    
        if self.windvectors and ds is not None:
            U, V = self._get_wind_data(ds, timestep)
            U, V = U.mean(), V.mean()
            mean_wind = np.sqrt(U*U + V*V)
            mean_temp = ds['t'][timestep, self.hlim[0]:self.hlim[1], :, :].mean(('x', 'y', 'z')).values
            self.plot_text = self.fig.text(0.4, 0.04 ,f'MEAN  w = {mean_wind:.2f} m/s   temp = {mean_temp:.2f} °K')

        self.fig.subplots_adjust(hspace=0.3, wspace=0.15, left=0.13, right=0.9, top=0.82, bottom=0.17)  # make space for colorbar
        

    def _init_single_spectra_plot(self, timestep: int | list, it_step: int = 3) -> None:
        """
        Initializes the single_spectra plot.

        Args:
            timestep (int | list): The time step for the plot.
        """

        cmap = plt.cm.Spectral_r
        # Convert to list for uniformity, if only a single timestep is provided
        if not isinstance(timestep, list):
            timestep = [timestep]
        
        z, x, y = self.hlim[0], self.idXY[0], self.idXY[1]  # Hardcoded values; can be changed

        radius = self.rgrenz * 1.0e6  # m to µm
        fac = 1/np.diff(np.log10(radius)) * 1.0e-6  # m3 to cm3

        if isinstance(timestep, list):
            it_start, it_end = timestep
        else:
            it_start, it_end = 0, 1

    
        self.fig.text(0.5, 0.09, 'radius [µm]', ha='center', va='center')
        self.fig.text(0.1, 0.5, 'dN/dlog(D_p) [#/cm]', ha='center', va='center', rotation='vertical')

        ds = None
        for (name, ds), ax in zip(self.datasets.items(), self.axes.flat):
            self.fig.suptitle(f'Temporal Evolution of {ds["nf"].attrs["long_name"]}\n@{self.y[y]:.3f}/{self.x[x]:.3f} and z = {self.height[z]:.3f} [km]',fontsize = 10, fontweight = 'bold' )
            ax.minorticks_on()
            ax.grid(True, which='major', linestyle='-', linewidth='0.5', color='black', alpha=0.5)
            ax.grid(True, which='minor', linestyle=':', linewidth='0.5', color='black', alpha=0.25)


            for it in range(it_start, it_end, it_step):
                irel = (1.*it - it_start) / (it_end - it_start)
                crgb = cmap(irel)  # Assuming there's more to this line...
                chex = colors.to_hex( crgb )
                    
                air_density = ds['rho'][it, z, y, x].values  # kg/m3
                spectra = ds[self.varname].isel(time = it, z = z, x = x, y = y, bin = slice(0, len(fac)))
                spectra = spectra * air_density * fac
                spectra.plot.step(
                    color = chex, 
                    yscale = 'log',
                    xscale = 'log',
                    ax = ax,
                    yincrease = False, 
                    ylim = self.ylim,
                    xlim = self.xlim,
                    alpha = 0.4
                ) # type: ignore
                
                datestr = str( self.time[it] )
                plt.figtext( 0.91, 0.85 - 0.8  * irel, datestr, c = chex, fontsize = 'small')
        
        self.fig.subplots_adjust(hspace=0.3, wspace=0.15, left=0.6, right=0.88, top=0.95, bottom=0.05)  # make space for colorbar


    def _get_data_area(self, ds: xr.Dataset, timestep: int | list, hmin: float, hmax: float) -> xr.DataArray:
        """
        Gets the data for the area plot.

        Args:
            ds (xr.Dataset): The data set.
            timestep (int): The time step.
            hmin (float): The minimum height.
            hmax (float): The maximum height.

        Returns:
            xr.DataArray: The data for the plot.
        """
        if self.timeframe == 'single':
            # Sum over the last two dimensions
            return ds[self.varname].isel(time=timestep, z=slice(hmin, hmax)).sum(dim=('z', 'bin'))
        elif self.timeframe == 'interval' and isinstance(timestep, list):
            # Sum over timeframe and the last two dimensions
            return ds[self.varname].isel(time=slice(timestep[0], timestep[-1]), z=slice(hmin, hmax)).sum(dim=('time', 'z', 'bin'))
        elif self.timeframe == 'all':
            # Sum over time and the last two dimensions
            return ds[self.varname].isel(z=slice(hmin, hmax)).sum(dim=('time', 'z', 'bin'))
        else:
            raise ValueError('timeframe must be in ["single", "interval", "all"]')
        
    def _get_wind_data(self, ds: xr.Dataset, timestep: int | list) -> tuple[np.ndarray, np.ndarray]:
        """
        Gets the data for the area plot.

        Args:
            ds (xr.Dataset): The data set.
            timestep (int): The time step.
            hmin (float): The minimum height.
            hmax (float): The maximum height.

        Returns:
            xr.DataArray: The data for the plot.
        """
        
        z_res = np.abs(np.diff(self.height[self.hlim[0]:self.hlim[1]+1]))
        if self.timeframe == 'single':
            # wind U, V are given as wind * air_density * z_res (kg/s)
            rho = ds['rho'].isel(time=timestep, z=slice(self.hlim[0], self.hlim[1]))
            U = ds['ut'].isel(time=timestep, z=slice(self.hlim[0], self.hlim[1]))
            V = ds['vt'].isel(time=timestep, z=slice(self.hlim[0], self.hlim[1]))
            U = U / z_res[:, np.newaxis, np.newaxis] / rho * 1.0e-6
            V = V / z_res[:, np.newaxis, np.newaxis] / rho * 1.0e-6
            return U.mean('z').values, V.mean('z').values
        else:
            raise ValueError('Only single frames! timeframe="signle"')


    def add_ruler(self, lat_start, lon_start, lat_end, lon_end):

        def _numticks(total_distance, num, alpha=0.8, labels=False):
            kwargsBL = {'alpha': 1, 'linewidth': 0.75, 'color': 'black'}
            kwargsWH = {'alpha': 1, 'linewidth': 1.5, 'color': 'white'}
            
            # Calculate the number of ticks based on 500m intervals
            num_ticks = int(total_distance / num)

            for ax in self.axes.flatten():
                # Draw the ruler to preset figure
                ax.plot([lon_start, lon_end], [lat_start, lat_end], **kwargsWH)
                ax.plot([lon_start, lon_end], [lat_start, lat_end], **kwargsBL)

                # Add ticks at an interval of 500 meters
                for i in range(num_ticks + 1):
                    fraction = i / num_ticks
                    lat_tick = lat_start + fraction * (lat_end - lat_start)
                    lon_tick = lon_start + fraction * (lon_end - lon_start)
                    # Calculate a small perpendicular offset for the tick marks
                    offset_lat = 0.3 * (lon_end - lon_start) / num_ticks
                    offset_lon = -0.3 * (lat_end - lat_start) / num_ticks
                    ax.plot(
                        [lon_tick + offset_lon, lon_tick - offset_lon], 
                        [lat_tick + offset_lat, lat_tick - offset_lat],
                        **kwargsWH
                        )
                    ax.plot(
                        [lon_tick + offset_lon, lon_tick - offset_lon], 
                        [lat_tick + offset_lat, lat_tick - offset_lat],
                        **kwargsBL
                        )
                
                    if labels:
                        lab = ax.text((lon_tick + offset_lon - 0.005), (lat_tick + offset_lat + 0.005),
                                f'{num_ticks-i:.0f}', fontsize=6, fontweight='bold')
                        lab.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])

        # Calculate the total distance between start and end points
        total_distance = haversine_distance(lat_start, lon_start, lat_end, lon_end)
        _numticks(total_distance, 500) # tick every 500 m
        _numticks(total_distance, 1000, labels=True) # tick every 1000 m

    def add_windvectors(self, ds, ax, timestep=0):
        u, v = self._get_wind_data(ds, timestep)
        quiv = ax.quiver(self.x[::2], self.y[::2], u[::2, ::2], v[::2, ::2], 
                         scale=1.0e2, color='black', alpha=0.7)
        return quiv




    def display(self, timestep: int = 0, x: int = 12, y: int = 12, z: int = 41, hmin: float = 42, hmax: float = 44, title: str='') -> None:
        # Displays the plot according to the plot mode.
        """

        Args:
            timestep (int, optional): The time step for the plot. Defaults to 0.
            x (int, optional): The x-value for the plot. Defaults to 12.
            y (int, optional): The y-value for the plot. Defaults to 12.
            hmin (float, optional): The minimum height. Defaults to 42.
            hmax (float, optional): The maximum height. Defaults to 44.
        """
        if title:
            self.fig.suptitle(title, weight='bold')
        if self.mode == 'profile':
            self._display_profile(timestep=timestep, x=x, y=y)
        elif self.mode =='single_spectra':
            self._display_single_spectra(timestep=timestep, x=x, y=y, z=z)
        elif self.mode =='timeseries':
            self._display_timeseries(x=x, y=y)    
        elif self.mode == 'area':
            self._display_area(timestep=timestep, hmin=hmin, hmax=hmax)

        if self.show:
            plt.show()


    def _display_profile(self, timestep: int = 0, x: int = 12, y: int = 12) -> None:
        """
        Displays the profile plot.

        Args:
            timestep (int, optional): The time step for the plot. Defaults to 0.
            x (int, optional): The x-value for the plot. Defaults to 12.
            y (int, optional): The y-value for the plot. Defaults to 12.
        """
        for i, (name, ds) in enumerate(self.datasets.items()):
            if i < self.n_axes and self.varname in ds:
                data = ds[self.varname].isel(time = timestep, y = y, x = x)
                self.pmeshs[i].set_array(data.values.ravel())


    def _display_single_spectra(self, timestep: int = 0, x: int = 12, y: int = 12, z: int = 41) -> None:
        pass


    def _display_area(self, timestep: int | list, hmin: float, hmax: float) -> None:
        """
        Displays the area plot.

        Args:
            timestep (int): The time step for the plot.
            hmin (float): The minimum height.
            hmax (float): The maximum height.
        """

        for i, (name, ds) in enumerate(self.datasets.items()):
            if i < self.n_axes and self.varname in ds:
                data = self._get_data_area(ds, timestep, hmin, hmax)
                self.pmeshs[i].set_array(data.values.ravel())  
                if self.windvectors: 
                    U, V = self._get_wind_data(self.datasets[name], timestep)
                    self.quivers[i].set_UVC(U[::2, ::2], V[::2, ::2])   

            if self.windvectors and i == 0 and self.plot_text:
                mean_wind = np.sqrt(U*U + V*V).mean()
                mean_temp = ds['t'][timestep, self.hlim[0]:self.hlim[1], :, :].mean(('x', 'y', 'z')).values
                self.plot_text.set_text(f'MEAN  w = {mean_wind:.2f} m/s   temp = {mean_temp:.2f} °K')


    def _display_timeseries(self, x: int = 12, y: int = 12) -> None:
        """
        Displays the timeseries plot.

        Args:
            timestep (int): The time step for the plot.
            hmin (float): The minimum height.
            hmax (float): The maximum height.
        """
        for i, (name, ds) in enumerate(self.datasets.items()):
            if i < self.n_axes and self.varname in ds:
                data = ds[self.varname].isel(y = y, x = x).sum(dim='bin').T
                self.pmeshs[i].set_array(data.values.ravel())
                self.axes[i].scatter(x, y, marker='x', s=5, markercolor='black')
 
        

    def interactive(self):
        """
        Creates an interactive plot with sliders according to the plot mode.
        """
        if self.mode in ['profile', 'timeseries']:
            sliders = {
                'timestep': widgets.IntSlider(min=0, max=len(self.time)-1, step=1, value=0) if self.mode == 'profile' else fixed(False),
                'x': widgets.IntSlider(min=0, max=len(self.x)-1, step=1, value=12), 
                'y': widgets.IntSlider(min=0, max=len(self.y)-1, step=1, value=12), 
                'hmin': fixed(False),
                'hmax': fixed(False),
                'savefig': fixed(False)
            }    
        elif self.mode == 'area':
            if self.timeframe == 'single':
                timestep_slider = widgets.IntSlider(min=0, max=len(self.time)-1, step=1, value=self.timestep)
            elif self.timeframe == 'interval':
                timestep_slider = widgets.IntRangeSlider(
                    value=self.timestep,
                    min=0,
                    max=len(self.time),
                    step=1,
                    description='Interval',
                    disabled=False,
                    orientation='horizontal',
                    readout=True,
                    readout_format='d',
                )
            elif self.timeframe == 'all':
                timestep_slider = fixed('all')
            else:
                raise ValueError('Wrong timeframe. Use "single", "interval", or "all"!')

            sliders = {
                'timestep': timestep_slider,
                'x': fixed(False),
                'y': fixed(False),
                'hmin': widgets.IntSlider(min=0, max=len(self.height)-1, step=1, value=42),
                'hmax': widgets.IntSlider(min=0, max=len(self.height)-1, step=1, value=44),
                'savefig': fixed(False),
            }
        else:
            sliders = {}
        
        interact(self.display, **sliders)
        return self


    def save_figure(self, filename: str) -> None:
        """
        Saves the figure to a file.

        Args:
            filename (str): The file name.
        """
        self.fig.savefig(filename, facecolor='white', dpi=300)


