#! /bin/bash
echo Job $SLURM_JOB_ID
echo Array index $SLURM_ARRAY_TASK_ID

mapfile -t files <  $1

#mapfile -t trees < $2

index=$(( SLURM_ARRAY_TASK_ID - 1 ))
k=${files[$index]}
#t=${trees[$index]}
if [ ! -f $k.BUSTED.json ]; then
    (echo 1; echo $k; echo y;echo 1;) | ~/bin/hyphy/HYPHYMP  ~/BUSTED-SRV/HBL/BUSTED-SRV.bf 
fi
