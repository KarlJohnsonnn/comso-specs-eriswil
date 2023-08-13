
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


# Python quicklook and analysis toolbox

# MultiPanelPlot Class

The `MultiPanelPlot` class is designed to facilitate the visualization of multi-panel plots for various types of data, including profiles, time series, and area plots. It provides interactive capabilities to explore the data, and it supports customization of plot parameters.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [API Reference](#api-reference)
- [Contributing](#contributing)
- [License](#license)

## Installation

To use the `MultiPanelPlot` class, you need to import it from the appropriate module in your code. You can also clone this repository to access the class directly.

## Usage

The `MultiPanelPlot` class is initialized with the following parameters:

- `datasets` (dict of str: xr.Dataset): A dictionary of data sets for plotting.
- `varname` (str): The variable name in the data sets.
- `mode` (str): The plotting mode ('profile', 'timeseries', or 'area').
- Other optional parameters for customization, such as `timestep0`, `timeframe`, `vmin`, `vmax`, etc.

Once initialized, you can use the `display()` method to render the plot according to the chosen mode. Additionally, you can use the `interactive()` method to create interactive plots with sliders for exploring the data.

## Examples

Here are a few examples of how to use the `MultiPanelPlot` class:

### Example 1: Profile Plot

```python
# Import necessary modules and classes
import xarray as xr
from your_module import MultiPanelPlot

# Create example datasets
datasets = {
    'Dataset1': xr.Dataset({'variable1': ...}),
    'Dataset2': xr.Dataset({'variable1': ...})
}

# Instantiate the MultiPanelPlot class
mpp = MultiPanelPlot(
    datasets=datasets,
    varname='variable1',
    mode='profile',
    timestep0=0,
    timeframe='single'
)

# Display an interactive profile plot
mpp.interactive()
```

## Example 2: Time Series Plot

```
# Instantiate the MultiPanelPlot class
mpp = MultiPanelPlot(
    datasets=datasets,
    varname='variable1',
    mode='timeseries',
    timestep0=0,
    timeframe='single'
)

# Display an interactive time series plot
mpp.interactive()
```

### Example 3: Area Plot

```
# Instantiate the MultiPanelPlot class
mpp = MultiPanelPlot(
    datasets=datasets,
    varname='variable1',
    mode='area',
    timestep0=0,
    timeframe='single'
)

# Display an interactive area plot
mpp.interactive()
```

---

