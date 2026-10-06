# Running this script runs the whole simulation study presented in --------
# Please note that the simulation requires a significant amount of runtime

source("1_data_simulation.R")
source("2_design_creation.R")
source("3_evaluation.R")
source("4_plots.R")

# Running the following lines create --------------------------------------
# 1.) The results in Section 5 (Study Example)

# The suppressWarnings() is used to prevent the warning stating that
# the folders containing the results have already been created beforehand
suppressWarnings(
  source("example_S.5.R")
)

# 2.) Figure 4 in Appendix A.1

source("appendix_A.1.R")

# 3.) Figure 10 in Appendix A.5

source("appendix_A.5.R")

# 4.) Table 6 in Appendix A.7

source("appendix_A.7.R")
