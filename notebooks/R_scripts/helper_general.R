
`%nin%` <- Negate('%in%')

# Running
# srun --job-name "InteractiveGPUJob" --gres=gpu:1 --cpus-per-task 1 --mem 5G --time 1:00:00 --pty bash

# library(devtools)
# library(glue)
# # Sys.setenv(DOWNLOAD_STATIC_LIBV8=1)
# library(ROpenCVLite)

# ROpenCVLite::installOpenCV()

# mamba env config vars set LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lib64



# mamba activate spatialomics_wb
# mamba install conda-forge::ocl-icd-system
#  install conda-forge::libv8

# devtools::install_github("swarm-lab/ROpenCVLite")


# Sys.getenv('LD_LIBRARY_PATH')


# home_dir <- "/central/groups/MazmanianLab/joeB"
# Sys.setenv(PATH = paste(
#     glue("{home_dir}/software/mambaforge/envs/spatialomics/bin"),
#     Sys.getenv("PATH"),
#     sep = ":"
# ))



# export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH


# devtools::install_github("tractatus/wholebrain")
# devtools::install_github("mjin1812/SMART")

# WARNING: No ICDs were found. Either,
# - Install a conda package providing a OpenCL implementation (pocl, oclgrind, intel-compute-runtime, beignet) or 
# - Make your system-wide implementation visible by installing ocl-icd-system conda package. 



# ROpenCVLite::isOpenCVInstalled()



# library(rstan)
# model<-'data{
# int N;
# real y[N];
# }
# parameters{
# real mu;
# real sigma;
# }
# model{
# y ~ normal(mu, sigma);
# }'


# model_data <- list( y = rnorm(10), N = 10  )
# fit <- stan(model_code = model, data = model_data, iter = 4000, chains =4)
# la <- extract(fit)
# hist(la$mu)


# library(rgl)
# options(rgl.printRglwidget = TRUE)
# rgl::plot3d(x=rnorm(1000), y=rnorm(1000), z=rnorm(1000))
