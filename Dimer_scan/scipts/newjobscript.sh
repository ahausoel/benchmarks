#!/bin/bash
#@ job_type = MPICH
#@ class = micro
#@ node = 1
#@ island_count=1
#@ total_tasks= 56
#@ wall_clock_limit = 02:40:00
##@ wall_clock_limit = 01:20:00
##@ wall_clock_limit = 00:40:00
##@ wall_clock_limit = 00:20:00
##@ wall_clock_limit = 00:10:00
#@ job_name = mytest
#@ network.MPI = sn_all,not_shared,us
#####@ initialdir = $(home)/mydir
#@ output = job.$(schedd_host).$(jobid).err
#@ error =  job.$(schedd_host).$(jobid).out
#@ notification=never
####@ notify_user=youremail_at_yoursite.xx
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup of environment
#export PYTHONPATH=/home/hpc/pr94vu/di73miv/opt_newpython/lib/python2.7/site-packages
#module load python/2.7_intel
#module unload mpi.ibm
#module load mpi.intel
#mpiexec -n 56 /gpfs/work/pr94vu/di73miv/w2dynamics___patrik_alexander_merge/DMFT.py
source ~/.bash_profile
#mpirun -np 56 /gpfs/work/pr94vu/di73miv/w2dynamics___patrik_alexander_merge/cthyb

#export LD_LIBRARY_PATH=/home/hpc/pr94vu/di73miv/opt/lib:$LD_LIBRARY_PATH

#export PYTHONPATH=/home/hpc/pr94vu/di73miv/work/w2dynamics_github_cmake___original:$PYTHONPATH
export PYTHONPATH=/home/hpc/pr94vu/di73miv/work/w2dynamics_github___neu___cmake___REPRODAGAIN___KLON:$PYTHONPATH
export LD_PRELOAD="/home/hpc/pr94vu/di73miv/CODES/w2dynamics_superstatesampling___likegit/build/nfft_local/lib/libnfft3.so"

### for 1000 i need 20 minutes

time=8000
#time=4000
#time=2000
#time=1000
#time=500
#time=250
#time = 125

#mpiexec -n 56 ./cthyb -t $time
mpiexec -n 56 ./w2dyn_cthyb -t $time
