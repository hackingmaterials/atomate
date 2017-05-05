for i in `ls -d interpolation*`
do
    echo $i
    cd $i
    rm *orig
    mkdir inputs
    mkdir outputs
    for j in INCAR KPOINTS POSCAR POTCAR
    do
        mv $j inputs/.
    done
    for j in vasprun.xml CONTCAR OUTCAR WAVECAR CHGCAR
    do
        mv $j outputs/.
    done
    cd ../
done
