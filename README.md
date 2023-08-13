
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

* **TODO**:

    * provide example plots and scripts to generate them

---

