reps=$1

cd process_G 
for r in $(seq 1 $reps);
do
    julia make_tin_bosonic.jl ../in_files/hchain_beta2048 $r

    julia cross_val.jl free 1
done

