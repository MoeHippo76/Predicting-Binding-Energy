#!/bin/bash
rm *.agr

cd src
make clean
make
cp Hfbtho_r ../
cd ..

cp ./src/Hfbtho_r .

# --- Default parameters --- #
Nosc=4
Z=12
N=12
test= true


# --- Read parameters ---#
# run providing number of N Z Nosc as input and eventually -test for local
# ./run.sh -z Z -n N -e Nosc
while getopts z:n:e:t:m: option
do
    case "${option}"
    in
        z) Z=${OPTARG};;
        n) N=${OPTARG};;
        e) Nosc=${OPTARG};;
        t) test=true;;
        m) module=${OPTARG};; # not implemented for HFBTHO
    esac
done


#module load intelcuda
module load intel
export OMP_STACKSIZE=1g


start=0
stop=38

# need A-dep limits for q-moment
#0.7*A^(2.67)/3300
A=$Z+$N


#SC="$(echo "scale=10; 0.7*e(l($A)*2.6667)/3300*0.1*0.5"| bc -l)"

beta="$(echo "scale=10; 0.3"| bc)"

#SC="$(echo "scale=10; $beta*e(l($A)*1.3333)/110"| bc -l)"

SC="$(echo "scale=10; $beta*e(l($A)*1.3333)/(5.5*$stop)"| bc -l)"

#echo $SC
#echo "$(echo "scale=10; 0.7*e(l(48)*2.67)/3300"| bc -l)"

echo "HFBTHO run Z=$Z, N=$N, Nosc=$Nosc test=$test"

for (( II=start; II<=stop; II++ ))
do
  
    {
    echo "&HFBTHO_FILE"
    echo ' OUTFILE = "out'$II'.out" /'
    echo "&HFBTHO_GENERAL"
    echo " number_of_shells = "$Nosc", "
    echo " oscillator_length = -1.0, "
    echo " basis_deformation = 0.01,"
    echo " proton_number = $Z, neutron_number = $N, type_of_calculation = 1 /"
    echo "&HFBTHO_ITERATIONS"
    echo " number_iterations = 400, accuracy = 1.E-9, restart_file = 1 /"
    echo "&HFBTHO_FUNCTIONAL"
    echo " functional = 'UNE1', add_initial_pairing = T,"
    echo " type_of_coulomb = 2/"
    echo "&HFBTHO_PAIRING"
    echo " user_pairing = T, vpair_n = 0.1, vpair_p = 0.1,"
    echo " pairing_cutoff = 1460.0, pairing_feature = 1.0 /"
    echo "&HFBTHO_CONSTRAINTS"
    echo " lambda_values = 1, 2, 3, 4, 5, 6, 7, 8,"
    echo " lambda_active = 0, 1, 0, 0, 0, 0, 0, 0,"
    #echo " expectation_values = 0.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /"
    echo "expectation_values = 0.0, $(echo "scale=4; $SC*($II-13)"| bc ), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /  "
    echo "&HFBTHO_BLOCKING"
    echo " proton_blocking = 0, 0, 0, 0, 0, neutron_blocking = 0, 0, 0, 0, 0 /"
    echo "&HFBTHO_PROJECTION"
    echo " switch_to_THO = 0, projection_is_on = 0,"
    echo " gauge_points = 1, delta_Z = 0, delta_N = 0 /"
    echo "&HFBTHO_TEMPERATURE"
    echo " set_temperature = T, temperature = 0.05 /"
    echo "&HFBTHO_DEBUG"
    echo " number_Gauss = 40, number_Laguerre = 40, number_Legendre = 80,"
    echo " compatibility_HFODD = F, number_states = 1500, force_parity = T,"
    echo " print_time = 0 /"
    } > input$II.inp

    if $test; then
        # echo "./Hfbtho_r < input$II.inp > worker_${II}.out &:"
        # echo "Z=$Z, N=$N, Nosc=$Nosc test=$test"
        time ./Hfbtho_r < input$II.inp > worker_${II}.out &
    else
        time srun -p lu -A ${ACCOUNT} --reservation=${ACCOUNT} -t 00:59:00 ./Hfbtho_r < input$II.inp > worker_${II}.out &
    fi

done



echo ""
wait

echo ''
echo 'Collecting output in one file'
#rm evsb.agr
#rm qvsb.agr
#rm dft_mesh.dat
#rm ene_beta2.agr
#rm q4_beta2.agr
#rm q6_beta2.agr
echo ''

echo "r2      ","q2       ", "beta    ", "q4     ", "q6     ", "q8     ", "ene     " >> dft_mesh.dat


for (( II=start; II<=stop; II++ ))
do
    
    x="$(grep "beta2" out$II.out | awk -F"beta2......................................" '{print $NF}')"
    y1="$(grep "tEnergy" out$II.out | awk -F"(qp)...................................." '{print $NF}')"
    y2="$(grep "quadrupole moment" out$II.out | awk -F"[b]...................................." '{print $NF}')"
    y3="$(grep "hexadecapole moment" out$II.out | awk -F"........................................................." '{print $NF}')"
    y4="$(grep "q6" out$II.out | awk -F"..........................................................." '{print $NF}')"
    y5="$(grep "q8" out$II.out | awk -F"..........................................................." '{print $NF}')"
    y6="$(grep "rms-radius .........." out$II.out | awk -F"..........................................................." '{print $NF}')"



    #echo  $x  $y1 >> evsb.agr
    #echo  $x  $y2 >> qvsb.agr
    #echo  $x  $y3 >> qvse4.agr
    #echo  $x  $y4 >> qvse6.agr
    #echo  $x  $y5 >> qvse8.agr
    #echo  $x  $y6 >> qvsr2.agr

    echo $y6 $y2 $x $y3 $y4 $y5 $y1 "0 1 0 0 0 0" >> dft_mesh.dat

    echo $x $y1 >> ene_beta2.agr

    echo $x $y3 >> q4_beta2.agr

    echo $x $y4 >> q6_beta2.agr

done

mkdir ./data/

mv dft_mesh.dat ./data/.
mv ene_beta2.agr ./data/.


#echo ""
#for (( II=start; II<=stop; II++ ))
#do
#grep "beta2" out$II.out | awk -F"beta2......................................" '{print $NF}'
#done

#gfortran post.f90
#./a.out


#less eq2.out
