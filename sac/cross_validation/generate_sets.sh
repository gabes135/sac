folder=$1
reps=$2

cd process_G 
for r in $(seq 1 $reps);
do
    julia make_tin_cv.jl "../in_files/$folder" $r
done
