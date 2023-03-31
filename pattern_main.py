# a .py script containing the contents of pattern_main.ipynb sans extensive
# markdown comments.

# ##########################################################
# setup
# ##########################################################
# Load all of the libraries
from matplotlib import pyplot as plt
import pandas as pd


%matplotlib inline
%config InlineBackend.figure_format = 'retina'
plt.rcParams['figure.figsize'] = 12, 6

# easier to read displays in console
pd.set_option('display.max_columns', None)

# source the helper functions
%run helpers.py


# ##########################################################
# Load ESM data
# ##########################################################

fetch_pangeo_table().to_csv('pangeo_table.csv', index=False)
# in same directory
pangeo_data = pd.read_csv('pangeo_table.csv')


esm_list = pangeo_data.model.unique().copy()
print(esm_list)

# Specify the scenarios and variables of interest for scaling
exp_list = ['historical', 'ssp370',
            'ssp585', 'ssp126', 'ssp245',
            'ssp119', 'ssp460', 'ssp434'  ]

var_list = ['tas', 'pr', 'hurs']

# ##########################################################
# Loop to create the annual patterns.
# ##########################################################

OUTPUT_DIR = 'outputs/'

for esm in esm_list:
    for scn in exp_list:
        for variable in var_list:

            # annual patterns
            savename = (esm + '_' + scn + '_' + variable + '_annual_pattern.nc' )
            print(savename)
            try:
                ensemble_ds = do_ps(esm_name = esm,
                                    var_name = variable,
                                    exp_name  = scn,
                                    monthly_or_annual = 'annual',
                                    fit_intercept = True,
                                    save_resids = True,
                                    tgav_DIR = (OUTPUT_DIR + 'tgav'))

                if not ensemble_ds == None:
                    ensemble_ds[0].to_netcdf(OUTPUT_DIR + 'patterns/'+ savename)
                    ensemble_ds[1].to_netcdf(OUTPUT_DIR + 'residuals/' + esm + '_' + scn +
                                             '_' + variable + '_annual_pattern_resids.nc' )

                del(ensemble_ds)




            except:
                print('issue creating ' + savename)
                print('Specific print statement above if data does not exist in pangeo_data.csv.')
                print('And this exception would not have triggered.')
                print('Bigger issue to figure out.')



        # end for over variable
    #end for over experiment
# end for over esm


# ##########################################################
# Loop to create the monthly patterns
# ##########################################################


months = pd.DataFrame(data={
            'month': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            'month_name': ['jan', 'feb', 'mar', 'apr',
                           'may', 'jun', 'jul', 'aug',
                           'sep', 'oct', 'nov', 'dec']
        })

OUTPUT_DIR = 'outputs/'

for esm in esm_list:
    for scn in exp_list:
        for variable in var_list:

            savename = (esm + '_' + scn + '_' + variable + '_monthly_patterns' )
            print(savename)
            try:
                ensemble_ds = do_ps(esm_name = esm,
                                    var_name = variable,
                                    exp_name  = scn,
                                    monthly_or_annual = 'monthly',
                                    fit_intercept = True,
                                    save_resids=True)



                if not ensemble_ds == None:
                    for ind in range(1, 13):
                        month_nm = months.loc[(ind - 1), 'month_name']
                        ensemble_ds[(ind-1)][0].to_netcdf(OUTPUT_DIR + 'patterns/' + savename +
                                                          '_' + month_nm + '.nc')
                        ensemble_ds[(ind-1)][1].to_netcdf(OUTPUT_DIR + 'residuals/' +savename +
                                                          '_' + month_nm + '_resids.nc')

                del(ensemble_ds)




            except:
                print('issue creating ' + savename)
                print('Specific print statement above if data does not exist in pangeo_data.csv.')
                print('And this exception would not have triggered.')
                print('Bigger issue to figure out.')



        # end for over variable
    #end for over experiment
# end for over esm