This code will reproduce the analysis and figures in the Nature Communications, 2018 article, "Generalized Leaky Integrate-And-Fire Models Classify Multiple Neuron Types", by Teeter Et. Al.  This code was written by Corinne Teeter with guest appearances by Vilas Menon and Ramakrishnan Iyer and others on the software team (see the AUTHORS file for more info).  

This repository is provided “as is”.  However, pull requests will be considered on an individual basis with adequate unit tests.  Please inquire whether your desired contribution may be considered.  Of particular interest would be the conversion to python 3 and individual module unit tests.

Testing: there are no unit tests provided within this repository.  However, the scripts create figures which should match the figures in the manuscript.  If the structured data directories are remade, the scripts available in the “sanity_check” folder can be used to insure recreated files and values match the ones provided in this repository and illustrate differences in directory structures.

Although human data is not included in the manuscript, it can be analyzed with the code here by changing the data path from 'mouse_struc_data_dir' to 'human_struc_data_dir' within the files where specified.

Below are the directions for running the code in this repository.

* All units are in International System of Units (SI) unless otherwise noted.

* All options that can be specified will be denoted at the top of the scripts.

------STEP 1:  RUN CODE IN "create_data_dir" (CAN BE SKIPPED)------------------

The final structured data directory which is created by the code within this "create_data_dir" folder is used throughout the rest of the analysis. I include the final structured data directory in the root directory within this repository for convenience ("mouse_struc_data_dir" and "human_struc_data_dir").  For reproducing the exact analysis in the Teeter et al. publication, running the code in this folder to recreate the structured data directories is NOT RECOMMENDED (for the reasons stated below) and can be skipped.  Nonetheless, this code may be useful especially for repeating analysis on additional data that may be released on the Allen Institute Cell Types Database.  

The "mouse_struc_data_dir" or "human_struc_data_dir" structured data directories will contain subdirectories which are named with a specimen id and cre line convention.  Each of the subdirectories will include the files needed for the rest of the analysis code.  These files include: "*_preprocessor_values.json" (1 file), "*_GLIF*_neuron_config.json" (1 file for every available model), "*_GLIF*_exp_var_ratio_10ms.json" (1 file for every available model), "*_GLIF*_subthr_v.json" (1 file for every available model).  For example please see the "mouse_struc_data_dir" or "human_struc_data_dir”. 

Recreating the files in the structured data directory is NOT recommended to simply recreate the manuscript figures for the following reasons: 
* They can take a very long time to run (~one week on a single processor) due to model simulation (we utilized a cluster). 
* The .nwb data files consume substantial space (mouse: 39G, human: 4G). 
* If a cluster is utilized file paths and script may need to be changed.
* If this code is used it may not produce the exact results of the paper because
	**the data in the database may include additional data.
	**you may be using a different version of any of the python packages
* A manual download and structuring of the files for one neuron will be required (as discussed below).  

One neuron assessed in the Teeter et al. publication was not published in the standard Allen Cell Types database (specimen id = 580895033). Although it passed quality control, it was positive for two cre-lines (which are not yet released in the database).  This cell can be accessed at http://download.alleninstitute.org/informatics-archive/september-2017/mouse_cell_types/glif/. Note that the number in this url and the file names in this directory do not correspond to the specimen id. Please follow the directions in step #3 below concerning the placement and renaming the files after download.

Directions:

* An internet connection will be required in order to access the Allen Institute Cell Types Database.

1.  Run "put_neuron_config_in_folder.py" 
	* This code grabs the configuration files and preprocessor files from the Allen Cell Types Database and puts them in a local directory named either "mouse_struc_data_dir" or "human_struc_data_dir".  
	* Note that although human data is is not included in the paper, it can be analyzed here.
	* A file named cell_types_manifest.json will be created in the root directory.

2.  Run "download_nwb_files.py"
	*This grabs the "*.nwb" and "*_ephy_sweeps.json" data files from the Allen Institute Cell Types Database and saves them in a directory named "mouse_nwb" or "human_nwb", which is placed within the root directory. These folders are placed in the root directory because they will be reused by scripts in other folders and therefore reduce the amount of large data files that could be potentially redownloaded
	*Note: in Dec 2017 these downloads consume the following amount of space: 
		mouse: 671 cells for a total of 39G. 
		human: 157 cells for a total of 4G.
	
3. Manually download the neuron utilized in the manuscript but not accessible within the Allen Cell Types Database. 
	* Within the "mouse_nwb" directory (which was created in step 2), create a directory named "specimen_580895033" (the number in this name corresponds to the specimen id as opposed to the other numbers on the download page). Download "580895015_ephys.nwb" from http://download.alleninstitute.org/informatics-archive/september-2017/mouse_cell_types/glif/ and place it within the "specimen_580895033" directory.  RENAME THE FILE TO "ephys.nwb".  Copy the contents of the "ephys_sweeps.json" file from the same website and place it in a file named "ephys_sweeps.json" within the "specimen_580895033" directory.  For examples, look at the other directories within the "mouse_nwb" directory. 
	*In the "mouse_struc_data_dir" directory create a folder named “580895033_Rbp4-Cre_KL100”.  The files can be copied over from the "mouse_struc_data_dir/580895033_Rbp4-Cre_KL100” provided with this repository in the root directory or they can be downloaded from the website as follows:
		* Download one of the "*_preprocessor_values.json" files from the same website (Although there are 4 different "*_preprocessor_values.json" available on the website, they are all the same and only one is necessary) and name it "580895033_Chrna2-Cre_OE25_preprocessor_values.json".  
		* Download "587865629_neuron_config.json" and rename it to "580895033_Chrna2-Cre_OE25_GLIF1_neuron_config.json".
		* Download "587865631_neuron_config.json" and rename it to "580895033_Chrna2-Cre_OE25_GLIF2_neuron_config.json".
		* Download "587865633_neuron_config.json" and rename it to "580895033_Chrna2-Cre_OE25_GLIF3_neuron_config.json".
		* 
		* Download "587865637_neuron_config.json" and rename it to "580895033_Chrna2-Cre_OE25_GLIF5_neuron_config.json".

The resulting directory should look like the “580895033_Chrna2-Cre_OE25” corresponding directory provided in the “mouse_struc_data_dir” directory within the root directory of this repository. 

4. Run "check_sweeps_and_rm_folders.py" to remove all specimen id directories that do not have sufficient or erroneous noise sweeps or model files (note, later code will break if this step is not completed).

5. a) Run "calc_all_exp_var_RUN.py"
	*Runs simulations to calculate the explained variance ratio of noise 1 and noise 2 before and after optimization.
	*The output is stored in the structured data directory in files that end with “_exp_var_ratio_10ms.json”.
	*This file takes many days to run on a single processor.
    b) Check the resulting directory with the "compare_EV_from_calc_and_db.py" script to make sure there are not errors with the explained variance calculations.

6.  Run “exclude_models_from_dir.py
	*Removes neurons from the structured data directory based on exclusion criteria described in the manuscript.

7.  At this point, you should check the resulting directory with "compare_outer_dir_structure.py" to confirm you understand what differences there are in the directory structure (i.e. if more data is now being used).  These scripts are located in the “sanity_check” directory.  If you are satisfied, the directory can now be placed one directory upstream in the root directory for use by the rest of the analysis code.

NOTE: because the subthreshold voltage is not subject to any exclusion criteria, the "*_GLIF*_subthr_v.json" files are not created within this directory.  They are created via the code in "SMfig_15".

------FIGURE 1a: NOT AVAILABLE —————————

Figure 1a made by Stefan Mihalas.

------FIGURE 1b: RUN CODE IN "MTfig_1b"-------------------

1. Run 'make_and_save_model_data.py' this can take ~1 hour or more.  Unless already downloaded, an internet connection will be needed to download two .nwb files which will be placed in a directory named ‘mouse_nwb’.  Note that if the the structured data directories were recreated using the scripts within "create_data_dir" this code will reuse the .nwb files already downloaded. The code will create a folder named “pkl_data” with ~5G of pickled data within it. 


2. Create the figure by running 'mech_zoom_ex_fig_1.py'. This will use the data within the “pkl_data” directory created in step 1.

------FIGURE 2: NOT AVAILABLE----------------------------

Figure 2b made by Jim Berg.

------FIGURE 3a: RUN CODE IN "MTfig_3a_SMfig_7"----------------

Run 'El_vs_thr_and_distributions.py'. 

------FIGURE 3b: RUN CODE IN "MTfig_3b_SMfig_1"-----------------

Run 'spike_cutting.py'.  Close plots to see more generated figures.

------FIGURE 3c and 3d: RUN CODE IN "MTfig_3cd_SMfig_2"-----------------

Run 'linear_regression.py’.

Note that Ramakrishnan Iyer helped write the code involving linear regression and the fitting of the after spike currents.

------FIGURE 3e: RUN CODE IN "MTfig_3e_SMfig_3"------------------------

Run 'multiblip.py'.  Close first set of figures to see figures for second neuron.

------FIGURE 3f: RUN CODE IN "MTfig_3f_SMfig_3"------------------------

Run 'v_comp_of_th.py'

------FIGURE 4 a and b: RUN CODE IN figure_4-------------------------------

Figure 4a may be created by running 'cononical_rasters.py'.

For convenience the data for the explained variance plots are saved in the json_data folder.  The data can be created by running 'create_expVar_data_files.py'. Note that the output file names will be different to reduce confusion.

Figure 4b can be created by running 'plot_expVar_example_curves.py'.

------FIGURE 5: RUN CODE IN "MTfig_5_MTtable_3_SMFigs_8_9_10"-------------


Run "expVar_level_box_plots.py" to make figure from saved data

To remake data in "saved_data" folder, run "expVar_level_calc_stats.py"

Note that this code will print the Friedman values for the screen.  The value for all neurons (first entry in the printed list) is reported in the manuscript.

The distribution of difference values between the excitatory and inhibitory neurons as found in the main text of the article can be recalculated by running expVar_stats_of_differences.py.  In the output that is printed to the screen, the "level" refers to the differences between the levels specified. Statistical tests are described in the publication.  Note that the data file being read in can be bipassed by 

------FIGURE 6:-------------

VILAS

------FIGURE 7:-------------

VILAS

------SUPPLEMENTARY FIGURE 1: RUN CODE IN "MTfig_3e_SMfig_3"-------------------

As mentioned above

------SUPPLEMENTARY FIGURE 2: RUN CODE IN "MTfig_3cd_SMfig_2"-------------------

As mentioned above

------SUPPLEMENTARY FIGURE 3: RUN CODE IN "MTfig_3e_SMfig_3"-------------------

As mentioned above

------SUPPLEMENTARY FIGURE 4: RUN CODE IN "SMfig_4"-------------------

Run "aspects_of_spiking_vs_exp_var.py".  This file uses the "spikecut_standard_err.csv" file where the standard error generated via "spikecut_calc_err.py" is saved.

To recreate the "spikecut_standard_err.csv" run "spikecut_calc_err.py" and make sure to follow directions in the file.  This code takes ~ 2 to 3 hours to run.

------SUPPLEMENTARY FIGURE 5: NOT AVAILABLE-------------------

------SUPPLEMENTARY FIGURE 6: RUN CODE IN "SMfig_6"---------------------------

Note that in the manuscript I forgot to convert the explained variance ratio on the y-axis to percentage by multiplying by 100.

Run "expVar_level_box_plots.py" to make figure from saved data

To remake data in "saved_data" folder, run "expVar_level_calc_stats.py".  Note that this resulting output file will create the explained variance ratio to be in percentages.

------SUPPLEMENTARY FIGURE 7: RUN CODE IN "MTfig_3a_SMfig_7"-------------------

As mentioned above

------SUPPLEMENTARY FIGURE 8: RUN CODE IN "MTfig_5_MTtable_3_SMFigs_8_9_10"-------------------

As mentioned above

------SUPPLEMENTARY FIGURE 9: RUN CODE IN "MTfig_5_MTtable_3_SMFigs_8_9_10"-------------------

As mentioned above

------SUPPLEMENTARY FIGURE 10: RUN CODE IN "MTfig_5_MTtable_3_SMFigs_8_9_10"-------------------

As mentioned above

------SUPPLEMENTARY FIGURE 11: RUN CODE IN "SMfig_11"----------------------------

For convenience, the output created by the files in this folder are saved in 'aic_spike_times_noise1.csv' and 'aic_subthresh_v_noise1.csv'.

To make the figure run 'plot_aic_diff.py' with the corresponding data file desired specified in the code.

To make the 'aic_spike_times_noise1.csv' run 'aic_spike_times_noise1.py'.  Note that this will take several hours and will download all of the .nwb files (if not previously downloaded).

To remake the 'aic_subthresh_v_noise1.csv' run 'aic_subthresh_v_noise1.py'.  Note that the output file will be saved with a different name to avoid confusion.

------SUPPLEMENTARY FIGURE 12:---------------------

Vilas

------SUPPLEMENTARY FIGURE 13:---------------------

Vilas

------SUPPLEMENTARY FIGURE 14:---------------------

Vilas

------SUPPLEMENTARY FIGURE 15: RUN CODE IN "SMfig_15"--------------------

Run "sub_v_plots.py" to remake plots using data saved in "subthreshold_data.pkl".  

Run "calc_all_subthresh_v_diff.py" or "calc_all_subthresh_v_diff_RUN.py" (example of how to run on cluster) to remake data saved in "subthreshold_data.pkl".

------SUPPLEMENTARY FIGURE 16: RUN CODE IN "SMfig_16"--------------------

Run "contributions_to_clustering_performance.py".  Note that this code used the data saved in 'MTfig_5_MTtable_3_SMFigs_8_9_10/saved_data/stats_out.pkl')

------SUPPLEMENTARY FIGURE 17---------------------

Vilas

--------TABLE 1: NOT AVALABLE, HAND MADE -----------------------------------

--------TABLE 3: NOT AVALABLE, HAND MADE -----------------------------------

--------TABLE 3: RUN CODE IN "MTfig_5_MTtable_3_SMFigs_8_9_10"--------------

Run code as described under Figure 5.

The output of "expVar_level_calc_stats.py" printed to the screen is table 5 in LaTeX format.

--------SUPPLEMENTARY TABLE 1:  RUN CODE IN "SMtables_1_2_3"--------------------

Run "param_dists_table_and_cononical.py".  Table in LaTeX format will be printed to screen. Be Descriptions of best model for various percentile requirements will also be printed to the screen.  Note that in the Supplimentary Material the explained variance was not converted to percentages; however, they are in this code.

--------SUPPLEMENTARY TABLE 2 (2 AND 3 ARE TWO PARTS TO SAME TABLE:  RUN CODE IN "SMtables_1_2_3"--------------------

Run "param_dists_table_and_cononical.py".  Table in LaTeX format will be printed to screen.

--------SUPPLEMENTARY TABLE 4 (4 AND 5 ARE TWO PARTS TO SAME TABLE:  RUN CODE IN "SMtables_4_5_6_7"--------------------

Run "cluster_GLIFparam_table.py". Table in LaTeX format will be printed to screen.

--------SUPPLEMENTARY TABLE 6 (6 AND 7 ARE TWO PARTS TO SAME TABLE:  RUN CODE IN "SMtables_4_5_6_7"--------------------

Run "cluster_feature_table.py". Table in LaTeX format will be printed to screen.

--------Queries of biophysical explained variances reported on page 9 within the Disscussion within the article--------------------

Run "query_biophys_expVar.py" and outputs will be printed to the screen.

