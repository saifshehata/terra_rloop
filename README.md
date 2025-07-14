# terra_rloop

Run scripts in the bash directory one by one numerically in the order they are named in.
Then run the scripts in the R directory also one by one to generate the figures.

The R scripts mainly use the outputs of the bash scripts numbered 09-12.

You will need to update my paths to the ones for your own directories in the bash and R scripts (e.g. update the 'project_dir' variable in the bash scripts to the path to your own project directory, and update the 'setwd()' path in the R scripts to your own working directory).

Hope you enjoy reproducing this project!
