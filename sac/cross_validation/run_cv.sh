reps=$1
param=$2
# julia cross_val.jl $param $rep

for r in $(seq 1 $reps);
do
    julia cross_val.jl $param $r &

done