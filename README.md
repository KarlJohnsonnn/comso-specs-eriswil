
# Running the COSMO-SPECS Eriswill Testcase01 on Levante

## Science Background


* Eriswill Testcase01 is a realcase scenario that runs COSMO-SPECS simulations centered around a location in the swiss alps

* IC and BC data have been prepared based on a sequence of nested COSMO runs

* Eriswill is a place where a measurement campaign was conducted by ETH Zurich and TROPOS 
    
    * wintertime fog or low-level stratus which consists mainly of super-cooled liquid droplets was perturbed by ice nucleating particles

    * primary ice formation was investigated with a sophisticated set of measurement devices   

    ![](images/cloudlab.png)



## Prepare and Run Test Case

### Installation
For getting and installing the COSMO-SPECS source code, please follow the instruction here: 
- [COSMO SPECS Installation at DKRZ](../../docs/Installation-at-DKRZ.md)

### Test Case Location

```
cd testcases/RUN_ERISWILL_TESTCASE01/
```

### Get the Additional Data

* *download tar from zenodo* - This is a future option. Upload to zenodo needs to be done.

* shortcut copy it from local dir on levante

```
cp /work/bb1262/data/cosmo-specs/cosmo-specs-eriswill-testcase_ic-bc-cosin-data.tar.gz .
tar xzvf cosmo-specs-eriswill-testcase_ic-bc-cosin-data.tar.gz
```

### Test Run

* get your compiled executable
    ```
    cp ../../build/psbm_fd4_levante .
    ```
* set your account details in the run script, i.e.

    ```
    > grep -in account run_COSMO-SPECS_levante 
    5:#SBATCH --account=bb1262    #### CHANGE THIS TO YOUR ACCOUNT
    ```
* choose a domain size of COSMO:
  * 122x112
  * 42x32
  * 12x12 
  * set it in`run_COSMO-SPECS_levante`: E.g. `CASE='42x32'`
* submit job
    ```
    sbatch run_COSMO-SPECS_levante
    ```

### Output Check List

* Take a look at `*.log` & `*.err` files!

* Is output written into `COS_out/20230125_*x*` and in `Ew1_2023020809_*x*.nc`? 

* Meteogramm data output in `M_Eriswil_SPECS.nc`?


# MultiPanelPlot class
## Introduction
The MultiPanelPlot class is a versatile tool designed to visualize multi-panel plots, especially tailored for visualizing specific data types. The class is part of the vizz module and offers various methods and functionalities to customize and enhance the visual representation of data.

## Dependencies
- numpy
- xarray
- matplotlib
- ... [other libraries from vizz.py]

## Reading in the data
```python
import vizz

root_path = '/work/bb1262/user/schimmel/cosmo-specs-torch/'
test_case = 'cs-eriswil__20230821_172639'
data_path = root_path + f'/cosmo-specs/testcases/RUN_ERISWILL_TESTCASE01/{test_case}/'

with open(meta_pdata_path + f'{test_case}.json'ath) as f:
   metadata = json.load(f)

d3_files = [data_path + '3D_' + file + '.nc' for file in metadata.keys()]

data = vizz.read_file_list(d3_files, varnames=['nf', 'vt', 'ut', 'wt', 't', 'rho', 'qw'])
```

## Initialization
To create a multi-panel plot, initialize an instance of the MultiPanelPlot class as shown below:

```python
import vizz

plot = vizz.MultiPanelPlot(
    data,
    metadata=metadata,
    varname='nf',
    nrows=1,
    ncols=2,
    mode='area',
    vmin=1.0,
    vmax=1.0e4,
    ...
)
```

Parameters:

- `data`: Data to be visualized.
- `metadata`: Metadata associated with the data.
- `varname`: Name of the variable to be plotted.
- `nrows` and ncols: Number of rows and columns for the multi-panel plot.
- `mode`: Visualization mode. E.g., 'area', 'profile'.
- `vmin` and vmax: Minimum and maximum values for the color scale.
- ... [other parameters]

## Features & Methods
### Add Ruler to Plot
You can add a ruler to the plot to measure distances or highlight specific regions:

```python
plot.add_ruler(lat_start, lon_start, lat_end, lon_end)
```

### Display and Save the Plot
To display the plot:

```python
plot.display(timestep=i, title=tit)
```

### To save the plot to a file:

```python
plot.save_figure(f'/path/to/save/{str(i).zfill(3)}_nf.png')
```

## Examples
### Area Mode Visualization
```python
plot1 = vizz.MultiPanelPlot(
    data,
    metadata=metadata,
    varname='nf',
    nrows=1,
    ncols=2,
    mode='area',
    ...
)
plot1.display(timestep=i, title=tit)
```

### Profile Mode Visualization
```python

plot2 = vizz.MultiPanelPlot(
    data,
    metadata=metadata,
    varname='nf',
    nrows=1,
    ncols=2,
    mode='profile',
    ...
)
```