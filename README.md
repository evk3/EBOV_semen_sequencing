# EBOV_semen_sequencing
Collection of scripts used for "Active Ebola virus Replication and Heterogeneous Evolutionary Rates in EVD Survivors." 

The samogitia.py script was originally created by Gytis Dudas.  Please visit his repository for the original version (https://github.com/evogytis).

Two modified samogitia.py scripts were used to calculate the EBOV rate estimates for EBOV semen and acute sequences provided by EVD survivors.  For the Western African cohort, we used samogitia_subrate_blood_lineage_V1.py and for the US EVD cohort we used samogitia_subrate_intermingle_lineage_V2_no_comments.py.

Usage and overview (Western African Cohort):

	python ./samogitia_subrate_blood_lineage_V1.py -t <logcombined_with_burn-in_removed_tree_file> -a subrate --substr p9_semen p19_semen p41_semen p74_semen p4_semen p90_semen -o <output_text_file>  

This script will find the shared most recent common ancestor for sequences containing the provided substrings ("leaves").  It will trace the lineage path from the provided leaf to tMRCA for all leaves with the provided substring.  It will remove duplicate paths to tMRCA and calculate the evolutionary subrate on unique branches of a clade.  Example: We generated two viral sequences from survivor 74.  This script will find their shared MRCA and ignore any other branches that might be on that subclade.  It will only calculate the evolutionary rate on unique branches shared between tMRCA and the two survivor 74 sequences.  Confused?  Shoot me an email: evk3@cdc.gov.    

Arguments: 

	-t <logcombined_with_burn-in_removed_tree_file> 
This script accepts a tree file from beast as input.  We used branch lengths in units of time and not subs/site.  We also logcombined and removed the burn-in from several beast runs to make the tree file (using Andrew Rambaut's logcombiner). 

	-a subrate --substr p9_semen p19_semen p41_semen p74_semen p4_semen p90_semen
This argument indicates that we want to calculate the evolutionary subrate for sequence names containing the listed substrings (--substr).  This script was modified from the original samogitia script to accept regular expressions as input.  Thus, all sequence names from p9 contain the substring "p9_semen"...etc.     

	-o <output_text_file>  This argument specifies the name of the output file.  
	
Output will be a tab-delimited text file containing columns for each of the provided substrings.  The rows are the states from the tree file.  

Using the substring above, you will have the following columns: 

	state p19_semen p41_semen p4_semen  p74_semen p90_semen p9_semen  Combo Other Full  
	10000 0.001205 0.000766 0.000603 0.000810 0.001122 0.00792 0.000877 0.000973 0.00096.  
	20000 ...     
	
"Combo" column is the combined rate estimate from all of the semen substrings.
"Other" column is the overall rate estimate minus the rate from the semen samples.
"Full" column is the overall evolutionary rate estimate for all branches in the tree.

Note: We used a modified version of Gytis' ipython notebook to read the output text file and generate the rate distributions and 95% HPD tails (*.ipnyb script included as supplmenentary files with Park et al. (2015) Cell 161(7):1516-26).



Usage and overview (US EVD Cohort):

	python ./samogitia_subrate_intermingle_lineage_V2_no_comments.py -t <logcombined_with_burn-in_removed_tree_file> -a subrate --substr LBR_pA_plasma LBR_pA_semen SLE_pC_blood "SLE_pC_semen|SLE_pC_urine" "GIN_pE_blood|GIN_pE_plasma" GIN_pE_semen -o <output_text_file>
  
This script will find the shared most recent common ancestor for blood and semen sequences.  It will trace the lineage path from the provided leaf to tMRCA for all leaves with the provided substring (ex - "LBR_pA_plasma" and "LBR_pA_semen").  It will remove shared pathes between the acute and semen sequences and assign any shared paths to the acute sequences (ie - the assumption is that unique mutations arose first in the blood and were transmitted onto the semen sequences during viral persistence).  

A note of caution before using this script: The substring order matters!  Substrings need to be entered as <blood_p1> <semen_p1> <blood_p2> <semen_p2>.  Regular expressions for sequence names is accepted, just put the regular expressions in double quotes.  Since we couldn't untagle the mutations shared between the semen and urine sequences, we calculated the evolutionary rate on unique branches shared between the semen and urine as a single rate.  Also, the calculations for the evolutionary rate for the eye sequence is hard-coded and not entered as an argument.  Confused?  Shoot me an email: evk3@cdc.gov.  
  
Arguments:

	-t <logcombined_with_burn-in_removed_tree_file> 
This script accepts a tree file from beast as input.  We used branch lengths in units of time and not subs/site.  We also logcombined and removed the burn-in from several beast runs to make the tree file (using Andrew Rambaut's logcombiner).

	-a subrate --substr LBR_pA_plasma LBR_pA_semen SLE_pC_blood "SLE_pC_semen|SLE_pC_urine" "GIN_pE_blood|GIN_pE_plasma" GIN_pE_semen
This argument indicates that we want to calculate the evolutionary subrate for sequence names containing the listed substrings (--substr).  This script was modified from the original samogitia script to accept regular expressions as input.  Thus, all acute sequence names from pA contain the substring "LBR_pA_plasma"...etc.
    
	-o <output_text_file>  
This argument specifies the name of the output file.  Output will be a tab-delimited text file containing columns for each of the provided substrings.  The rows are the states from the tree file.  Using the substring above, you will have the following columns:

	state GIN_pE_blood|GIN_pE_plasma GIN_pE_semen LBR_pA_plasma  LBR_pA_semen SLE_pC_blood SLE_pC_semen|SLE_pC_urine	Eye_branch_rate  Blood_Combo	Semen_Combo Acute_rate Full	
	10000	0.001205	0.000766	0.000603	0.000810	0.001122	0.00792	0.000877	0.000973	0.00096.	
	20000	...			
		
"Eye_branch_rate" column is the evolutionary rate for the single eye viral sequence.
"Blood_combo" column is the overall rate estimate for all of the acute viral sequences.
"Semen_ombo" column is the overall evolutionary rate estimate for all of the persistent viral sequences.
"Full" column is the overall evolutionary rate estimate for all sequences in the tree.

Note: We used a modified version of Gytis' ipython notebook to read the output text file and generate the rate distributions and 95% HPD tails (*.ipnyb script included as supplmenentary files with Park et al. (2015) Cell 161(7):1516-26).
