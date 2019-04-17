# first step is to make the blueprint file. Did that by manualling editing the example. Used the aversaged folded sfs from 20 bootstraps.

# Then generate the batch file
java -cp stairway_plot_v2/stairway_plot_es Stairbuilder raso_two-epoch_fold.blueprint

# then run the batch file
bash raso_two-epoch_fold.blueprint.sh
 
