#!/bin/bash

basename=benzene
mols=('c1ccccc1.smi')
nmols=(4)
spg=14
diag=false

symprec=0.1

nstruc=50
nstrucbatch=50
factor=1.1
tfactor=1.0

qc_method=DFTB
opt_method=LBFGS
opt_maxcycle=30
opt_maxstepsize=0.01
opt_fmax=0.1

system=ito
qsub_cmd=pjsub
time="1:00:00"

#jobtype=initonly
#jobtype=optonly
jobtype=all

. ${HOME}/miniconda3/etc/profile.d/conda.sh

conda activate py38-cryspy

nbatch=`expr ${nstruc} / ${nstrucbatch}`
fac=`printf "%.2f" ${factor}`
tfac=`printf "%.2f" ${tfactor}`
fmx=`printf "%.4f" ${opt_fmax}`
mstep=`printf "%.5f" ${opt_maxstepsize}`

fmols=""
fnmols=""
for ((i=0; i < ${#mols[@]}; i++)); do
    fmols="${fmols},'${mols[$i]}'"
    fnmols="${fnmols},${nmols[$i]}"
done
fmols="[${fmols:1}]"
fnmols="[${fnmols:1}]"

diag=${diag,,}
diagsub=""
if [ ${diag} = 'true' ]; then
    diagsub="-diag"
fi

in_file=${basename}_spg${spg}${diagsub}_factor${fac}_tfactor${tfac}_${qc_method}_${opt_method}_opt_fmax${fmx}_maxstep${mstep}.yaml
gendir=factor${fac}_tfactor${tfac}
spgdir=${basename}_spg${spg}${diagsub}
wrkdir=./${spgdir}/${gendir}/${qc_method}/${opt_method}/opt_fmax${fmx}_maxsteps${mstep}
if [ -d ${wrkdir} ]; then
    echo rm -rf ${wrkdir}
    rm -rf ${wrkdir}
fi

cat <<EOF > ${in_file}
basename: ${basename}
mols: ${fmols}
nmols: ${fnmols}
spg: ${spg}
diag: ${diag}

symprec: ${symprec}

nstruc: ${nstruc}
factor: ${factor}
tfactor: ${tfactor}

qc_method: ${qc_method}
opt_method: ${opt_method}
opt_fmax: ${opt_fmax}
opt_maxcycle: ${opt_maxcycle}
opt_maxstepsize: ${opt_maxstepsize}
EOF

script_file_template=${system}_${qsub_cmd}_template.sh

echo jobtype: ${jobtype}
if [ ${jobtype} = 'initonly' ]; then
    runopt="-noopt"
    script_file=${jobtype}_${basename}_spg${spg}${diagsub}_factor${fac}_tfactor${tfac}.sh
    out_file=${jobtype}_${basename}_spg${spg}${diagsub}_factor${fac}_tfactor${tfac}.log
    cat ${script_file_template} \
	| sed  -e "s/%TIME%/${time}/g" -e "s/%IN_FILE%/${in_file}/g" -e "s/%OUT_FILE%/${out_file}/g" \
	-e "s/%RUNOPT%/${runopt}/g" \
	> ${script_file}
    echo ${qsub_cmd} ${script_file}
    ${qsub_cmd} ${script_file}
elif [ ${jobtype} = 'optonly' ]; then
    for ((i=0; i < ${nbatch}; i++)); do
	istbgn=$((i*${nstrucbatch}+1))
	istend=$((istbgn+${nstrucbatch}-1))
	istbgnf=`printf "%.5d" ${istbgn}`
	istendf=`printf "%.5d" ${istend}`
	runopt="-nogen -istbgn ${istbgn} -istend ${istend}"
	script_file=${jobtype}_${basename}_spg${spg}${diagsub}_factor${fac}_tfactor${tfac}_${qc_method}_${opt_method}_opt_fmax${fmx}_maxstep${mstep}_nstruc${istbgnf}-${istendf}.sh
	out_file=${jobtype}_${basename}_spg${spg}${diagsub}_factor${fac}_tfactor${tfac}_${qc_method}_${opt_method}_opt_fmax${fmx}_maxstep${mstep}_nstruc${istbgnf}-${istendf}.log
	cat ${script_file_template} \
	    | sed  -e "s/%TIME%/${time}/g" -e "s/%IN_FILE%/${in_file}/g" -e "s/%OUT_FILE%/${out_file}/g" \
	    -e "s/%RUNOPT%/${runopt}/g" \
	    > ${script_file}
	echo ${qsub_cmd} ${script_file}
	${qsub_cmd} ${script_file}
    done
else
    runopt=""
    script_file=${jobtype}_${basename}_spg${spg}${diagsub}_factor${fac}_tfactor${tfac}_${qc_method}_${opt_method}_opt_fmax${fmx}_maxstep${mstep}.sh
    out_file=${jobtype}_${basename}_spg${spg}${diagsub}_factor${fac}_tfactor${tfac}_${qc_method}_${opt_method}_opt_fmax${fmx}_maxstep${mstep}.log
    cat ${script_file_template} \
	| sed  -e "s/%TIME%/${time}/g" -e "s/%IN_FILE%/${in_file}/g" -e "s/%OUT_FILE%/${out_file}/g" \
	-e "s/%RUNOPT%/${runopt}/g" \
	> ${script_file}
    echo ${qsub_cmd} ${script_file}
    ${qsub_cmd} ${script_file}
fi

conda deactivate

exit 0
