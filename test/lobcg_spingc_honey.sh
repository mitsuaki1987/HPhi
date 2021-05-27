#!/bin/sh -e

mkdir -p lobcg_spingc_honey/
cd lobcg_spingc_honey

cat > stan.in <<EOF
W = 2
L = 3
model = "SpinGC"
method = "CG"
lattice = "Honeycomb"
J0x = -1.0
J0y =  2.0
J0z =  3.3
J1x =  1.0
J1y = 1.0
J1z =  0.0
J2x =  3.0
J2y =  2.0
J2z = 1.0
2S=1
h=-0.1
exct = 3
EOF

${MPIRUN} ../../src/HPhi++ -s stan.in

# Check value

cat > reference.dat <<EOF
 0
  -10.6884732749673663
  0.0000000000000000
  0.0141890312794323

 1
  -9.9987732950471511
  0.0000000000000000
  0.2032909253553725

 2
  -9.4973428762800314
  0.0000000000000000
  0.0040230516306462
EOF
paste output/zvo_energy.dat reference.dat > paste1.dat
diff=`awk '
BEGIN{diff=0.0}
{diff+=sqrt(($2-$3)*($2-$3))}
END{printf "%8.6f", diff}
' paste1.dat`
echo "Diff output/zvo_energy.dat : " ${diff}
test "${diff}" = "0.000000"

# Check one-body G

cat > reference.dat <<EOF
   0.4988175807 0.0000000000
   0.5011824193 0.0000000000
   0.4988175807 0.0000000000
   0.5011824193 0.0000000000
EOF
paste output/zvo_cisajs_eigen0.dat reference.dat > paste2.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($5-$7)*($5-$7)+($6-$8)*($6-$8))}
END{printf "%8.6f", diff/NR}
' paste2.dat`
echo "Diff output/zvo_cisajs_eigen0.dat : " ${diff}
test "${diff}" = "0.000000"

cat > reference.dat <<EOF
   0.4830590896 0.0000000000
   0.5169409104 0.0000000000
   0.4830590896 0.0000000000
   0.5169409104 0.0000000000
EOF
paste output/zvo_cisajs_eigen1.dat reference.dat > paste3.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($5-$7)*($5-$7)+($6-$8)*($6-$8))}
END{printf "%8.6f", diff/NR}
' paste3.dat`
echo "Diff output/zvo_cisajs_eigen1.dat : " ${diff}
test "${diff}" = "0.000000"

cat > reference.dat <<EOF
   0.4996647729 0.0000000000
   0.5003352271 0.0000000000
   0.4996647198 0.0000000000
   0.5003352802 0.0000000000
EOF
paste output/zvo_cisajs_eigen2.dat reference.dat > paste4.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($5-$7)*($5-$7)+($6-$8)*($6-$8))}
END{printf "%8.6f", diff/NR}
' paste4.dat`
echo "Diff output/zvo_cisajs_eigen2.dat : " ${diff}
test "${diff}" = "0.000000"

# Check two-body G

cat > reference.dat <<EOF
   0.4988175807 0.0000000000
   0.0000000000 0.0000000000
   0.1521301487 0.0000000000
   0.3466874321 0.0000000000
   0.2626521347 0.0000000000
   0.2361654460 0.0000000000
   0.2271325282 0.0000000000
   0.2716850525 0.0000000000
   0.3357833439 0.0000000000
   0.1630342368 0.0000000000
   0.1655249521 0.0000000000
   0.3332926286 0.0000000000
   0.2589835612 0.0000000000
   0.2398340195 0.0000000000
   0.2365579851 0.0000000000
   0.2622595956 0.0000000000
   0.3357833439 0.0000000000
   0.1630342368 0.0000000000
   0.0379153164 0.0000000000
   0.4609022643 0.0000000000
   0.2589835612 0.0000000000
   0.2398340195 0.0000000000
   0.2368368785 0.0000000000
   0.2619807022 0.0000000000
   0.4988175807 0.0000000000
   -0.1228020806 -0.0000000000
   0.0582427024 0.0000000000
   -0.1071075403 -0.0000000000
   0.0791503097 0.0000000000
   -0.0926361158 0.0000000000
   0.0926788613 0.0000000000
   -0.0860742771 -0.0000000000
   0.0791503097 0.0000000000
   -0.4006078544 0.0000000000
   0.0926788613 0.0000000000
   -0.0582576025 -0.0000000000
   0.5011824193 0.0000000000
   -0.1228020806 0.0000000000
   0.0582427024 -0.0000000000
   -0.1071075403 0.0000000000
   0.0791503097 -0.0000000000
   -0.0926361158 -0.0000000000
   0.0926788613 -0.0000000000
   -0.0860742771 0.0000000000
   0.0791503097 -0.0000000000
   -0.4006078544 -0.0000000000
   0.0926788613 -0.0000000000
   -0.0582576025 0.0000000000
   0.0000000000 0.0000000000
   0.5011824193 0.0000000000
   0.3466874321 0.0000000000
   0.1544949872 0.0000000000
   0.2361654460 0.0000000000
   0.2650169733 0.0000000000
   0.2716850525 0.0000000000
   0.2294973668 0.0000000000
   0.1630342368 0.0000000000
   0.3381481825 0.0000000000
   0.3332926286 0.0000000000
   0.1678897907 0.0000000000
   0.2398340195 0.0000000000
   0.2613483997 0.0000000000
   0.2622595956 0.0000000000
   0.2389228237 0.0000000000
   0.1630342368 0.0000000000
   0.3381481825 0.0000000000
   0.4609022643 0.0000000000
   0.0402801549 0.0000000000
   0.2398340195 0.0000000000
   0.2613483997 0.0000000000
   0.2619807022 0.0000000000
   0.2392017171 0.0000000000
   0.1521301487 0.0000000000
   0.3466874321 0.0000000000
   0.4988175807 0.0000000000
   0.0000000000 0.0000000000
   0.2271325282 0.0000000000
   0.2716850525 0.0000000000
   0.2626521347 0.0000000000
   0.2361654460 0.0000000000
   0.0379153164 0.0000000000
   0.4609022643 0.0000000000
   0.3357833439 0.0000000000
   0.1630342368 0.0000000000
   0.2368368785 0.0000000000
   0.2619807022 0.0000000000
   0.2589835612 0.0000000000
   0.2398340195 0.0000000000
   0.1655249521 0.0000000000
   0.3332926286 0.0000000000
   0.3357833439 0.0000000000
   0.1630342368 0.0000000000
   0.2365579851 0.0000000000
   0.2622595956 0.0000000000
   0.2589835612 0.0000000000
   0.2398340195 0.0000000000
   -0.1228020806 0.0000000000
   0.4988175807 0.0000000000
   -0.1071075403 -0.0000000000
   0.0582427024 -0.0000000000
   -0.4006078544 -0.0000000000
   0.0791503097 -0.0000000000
   -0.0582576025 0.0000000000
   0.0926788613 -0.0000000000
   -0.0926361158 -0.0000000000
   0.0791503097 -0.0000000000
   -0.0860742771 0.0000000000
   0.0926788613 0.0000000000
   -0.1228020806 -0.0000000000
   0.5011824193 0.0000000000
   -0.1071075403 0.0000000000
   0.0582427024 0.0000000000
   -0.4006078544 0.0000000000
   0.0791503097 0.0000000000
   -0.0582576025 -0.0000000000
   0.0926788613 0.0000000000
   -0.0926361158 0.0000000000
   0.0791503097 0.0000000000
   -0.0860742771 -0.0000000000
   0.0926788613 -0.0000000000
   0.3466874321 0.0000000000
   0.1544949872 0.0000000000
   0.0000000000 0.0000000000
   0.5011824193 0.0000000000
   0.2716850525 0.0000000000
   0.2294973668 0.0000000000
   0.2361654460 0.0000000000
   0.2650169733 0.0000000000
   0.4609022643 0.0000000000
   0.0402801549 0.0000000000
   0.1630342368 0.0000000000
   0.3381481825 0.0000000000
   0.2619807022 0.0000000000
   0.2392017171 0.0000000000
   0.2398340195 0.0000000000
   0.2613483997 0.0000000000
   0.3332926286 0.0000000000
   0.1678897907 0.0000000000
   0.1630342368 0.0000000000
   0.3381481825 0.0000000000
   0.2622595956 0.0000000000
   0.2389228237 0.0000000000
   0.2398340195 0.0000000000
   0.2613483997 0.0000000000
EOF
paste output/zvo_cisajscktalt_eigen0.dat reference.dat > paste5.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($9-$11)*($9-$11)+($10-$12)*($10-$12))}
END{printf "%8.6f", diff/NR}
' paste5.dat`
echo "Diff output/zvo_cisajscktalt_eigen0.dat : " ${diff}
test "${diff}" = "0.000000"

cat > reference.dat <<EOF
   0.4830590896 0.0000000000
   0.0000000000 0.0000000000
   0.1619768145 0.0000000000
   0.3210822751 0.0000000000
   0.2396884134 0.0000000000
   0.2433706762 0.0000000000
   0.2154021372 0.0000000000
   0.2676569524 0.0000000000
   0.2760350460 0.0000000000
   0.2070240435 0.0000000000
   0.2071300657 0.0000000000
   0.2759290238 0.0000000000
   0.2324423534 0.0000000000
   0.2506167361 0.0000000000
   0.2359534778 0.0000000000
   0.2471056118 0.0000000000
   0.2760350460 0.0000000000
   0.2070240435 0.0000000000
   0.0987613820 0.0000000000
   0.3842977076 0.0000000000
   0.2324423534 0.0000000000
   0.2506167361 0.0000000000
   0.2315833498 0.0000000000
   0.2514757398 0.0000000000
   0.4830590896 0.0000000000
   -0.2106169908 -0.0000000000
   0.1499570583 0.0000000000
   -0.1806799489 -0.0000000000
   0.1530331032 0.0000000000
   -0.1481810562 -0.0000000000
   0.1657767292 0.0000000000
   -0.1565344772 -0.0000000000
   0.1530331032 0.0000000000
   -0.3497358810 -0.0000000000
   0.1657767292 0.0000000000
   -0.1475302823 -0.0000000000
   0.5169409104 0.0000000000
   -0.2106169908 0.0000000000
   0.1499570583 -0.0000000000
   -0.1806799489 0.0000000000
   0.1530331032 -0.0000000000
   -0.1481810562 0.0000000000
   0.1657767292 -0.0000000000
   -0.1565344772 0.0000000000
   0.1530331032 -0.0000000000
   -0.3497358810 0.0000000000
   0.1657767292 -0.0000000000
   -0.1475302823 0.0000000000
   0.0000000000 0.0000000000
   0.5169409104 0.0000000000
   0.3210822751 0.0000000000
   0.1958586354 0.0000000000
   0.2433706762 0.0000000000
   0.2735702343 0.0000000000
   0.2676569524 0.0000000000
   0.2492839581 0.0000000000
   0.2070240435 0.0000000000
   0.3099168669 0.0000000000
   0.2759290238 0.0000000000
   0.2410118866 0.0000000000
   0.2506167361 0.0000000000
   0.2663241743 0.0000000000
   0.2471056118 0.0000000000
   0.2698352987 0.0000000000
   0.2070240435 0.0000000000
   0.3099168669 0.0000000000
   0.3842977076 0.0000000000
   0.1326432029 0.0000000000
   0.2506167361 0.0000000000
   0.2663241743 0.0000000000
   0.2514757398 0.0000000000
   0.2654651707 0.0000000000
   0.1619768145 0.0000000000
   0.3210822751 0.0000000000
   0.4830590896 0.0000000000
   0.0000000000 0.0000000000
   0.2154021372 0.0000000000
   0.2676569524 0.0000000000
   0.2396884134 0.0000000000
   0.2433706762 0.0000000000
   0.0987613820 0.0000000000
   0.3842977076 0.0000000000
   0.2760350460 0.0000000000
   0.2070240435 0.0000000000
   0.2315833498 0.0000000000
   0.2514757398 0.0000000000
   0.2324423534 0.0000000000
   0.2506167361 0.0000000000
   0.2071300657 0.0000000000
   0.2759290238 0.0000000000
   0.2760350460 0.0000000000
   0.2070240435 0.0000000000
   0.2359534778 0.0000000000
   0.2471056118 0.0000000000
   0.2324423534 0.0000000000
   0.2506167361 0.0000000000
   -0.2106169908 0.0000000000
   0.4830590896 0.0000000000
   -0.1806799489 -0.0000000000
   0.1499570583 -0.0000000000
   -0.3497358810 0.0000000000
   0.1530331032 -0.0000000000
   -0.1475302823 0.0000000000
   0.1657767292 -0.0000000000
   -0.1481810562 0.0000000000
   0.1530331032 -0.0000000000
   -0.1565344772 0.0000000000
   0.1657767292 -0.0000000000
   -0.2106169908 -0.0000000000
   0.5169409104 0.0000000000
   -0.1806799489 0.0000000000
   0.1499570583 0.0000000000
   -0.3497358810 -0.0000000000
   0.1530331032 0.0000000000
   -0.1475302823 -0.0000000000
   0.1657767292 0.0000000000
   -0.1481810562 -0.0000000000
   0.1530331032 0.0000000000
   -0.1565344772 -0.0000000000
   0.1657767292 0.0000000000
   0.3210822751 0.0000000000
   0.1958586354 0.0000000000
   0.0000000000 0.0000000000
   0.5169409104 0.0000000000
   0.2676569524 0.0000000000
   0.2492839581 0.0000000000
   0.2433706762 0.0000000000
   0.2735702343 0.0000000000
   0.3842977076 0.0000000000
   0.1326432029 0.0000000000
   0.2070240435 0.0000000000
   0.3099168669 0.0000000000
   0.2514757398 0.0000000000
   0.2654651707 0.0000000000
   0.2506167361 0.0000000000
   0.2663241743 0.0000000000
   0.2759290238 0.0000000000
   0.2410118866 0.0000000000
   0.2070240435 0.0000000000
   0.3099168669 0.0000000000
   0.2471056118 0.0000000000
   0.2698352987 0.0000000000
   0.2506167361 0.0000000000
   0.2663241743 0.0000000000
EOF
paste output/zvo_cisajscktalt_eigen1.dat reference.dat > paste6.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($9-$11)*($9-$11)+($10-$12)*($10-$12))}
END{printf "%8.6f", diff/NR}
' paste6.dat`
echo "Diff output/zvo_cisajscktalt_eigen1.dat : " ${diff}
test "${diff}" = "0.000000"

cat > reference.dat <<EOF
   0.4996647729 0.0000000000
   0.0000000000 0.0000000000
   0.0935013949 0.0000000000
   0.4061633780 0.0000000000
   0.3933203688 0.0000000000
   0.1063444041 0.0000000000
   0.0992809214 0.0000000000
   0.4003838515 0.0000000000
   0.3989496968 0.0000000000
   0.1007150760 0.0000000000
   0.1018391298 0.0000000000
   0.3978256431 0.0000000000
   0.3902232024 0.0000000000
   0.1094415705 0.0000000000
   0.1078640760 0.0000000000
   0.3918006969 0.0000000000
   0.3989496967 0.0000000000
   0.1007150762 0.0000000000
   0.0230618638 0.0000000000
   0.4766029091 0.0000000000
   0.3902232060 0.0000000000
   0.1094415669 0.0000000000
   0.1072224778 0.0000000000
   0.3924422951 0.0000000000
   0.4996647729 0.0000000000
   -0.0783664026 -0.0000000148
   0.0202135403 -0.0000000009
   -0.0902616774 -0.0000000162
   0.0308400971 0.0000000033
   -0.0153038564 -0.0000000147
   0.0344133921 0.0000000050
   -0.0061669390 -0.0000000150
   0.0308401049 -0.0000000042
   -0.2687060154 0.0000000078
   0.0344133657 -0.0000000006
   -0.0141481335 -0.0000000011
   0.5003352271 0.0000000000
   -0.0783664026 0.0000000148
   0.0202135403 0.0000000009
   -0.0902616774 0.0000000162
   0.0308400971 -0.0000000033
   -0.0153038564 0.0000000147
   0.0344133921 -0.0000000050
   -0.0061669390 0.0000000150
   0.0308401049 0.0000000042
   -0.2687060154 -0.0000000078
   0.0344133657 0.0000000006
   -0.0141481335 0.0000000011
   0.0000000000 0.0000000000
   0.5003352271 0.0000000000
   0.4061633250 0.0000000000
   0.0941719022 0.0000000000
   0.1063444070 0.0000000000
   0.3939908201 0.0000000000
   0.4003838007 0.0000000000
   0.0999514264 0.0000000000
   0.1007150762 0.0000000000
   0.3996201509 0.0000000000
   0.3978255907 0.0000000000
   0.1025096364 0.0000000000
   0.1094415708 0.0000000000
   0.3908936564 0.0000000000
   0.3918006371 0.0000000000
   0.1085345901 0.0000000000
   0.1007150749 0.0000000000
   0.3996201522 0.0000000000
   0.4766028555 0.0000000000
   0.0237323717 0.0000000000
   0.1094415645 0.0000000000
   0.3908936626 0.0000000000
   0.3924422401 0.0000000000
   0.1078929871 0.0000000000
   0.0935013949 0.0000000000
   0.4061633250 0.0000000000
   0.4996647198 0.0000000000
   0.0000000000 0.0000000000
   0.0992809201 0.0000000000
   0.4003837997 0.0000000000
   0.3933203127 0.0000000000
   0.1063444071 0.0000000000
   0.0230618664 0.0000000000
   0.4766028534 0.0000000000
   0.3989496423 0.0000000000
   0.1007150776 0.0000000000
   0.1072224839 0.0000000000
   0.3924422359 0.0000000000
   0.3902231471 0.0000000000
   0.1094415727 0.0000000000
   0.1018391271 0.0000000000
   0.3978255927 0.0000000000
   0.3989496464 0.0000000000
   0.1007150734 0.0000000000
   0.1078640744 0.0000000000
   0.3918006455 0.0000000000
   0.3902231541 0.0000000000
   0.1094415657 0.0000000000
   -0.0783664026 0.0000000148
   0.4996647198 0.0000000000
   -0.0902616673 0.0000000046
   0.0202135434 -0.0000000003
   -0.2687060152 0.0000000039
   0.0308401053 -0.0000000016
   -0.0141481365 0.0000000104
   0.0344133670 0.0000000036
   -0.0153038380 0.0000000017
   0.0308400945 -0.0000000014
   -0.0061669180 -0.0000000128
   0.0344133868 0.0000000005
   -0.0783664026 -0.0000000148
   0.5003352802 0.0000000000
   -0.0902616673 -0.0000000046
   0.0202135434 0.0000000003
   -0.2687060152 -0.0000000039
   0.0308401053 0.0000000016
   -0.0141481365 -0.0000000104
   0.0344133670 -0.0000000036
   -0.0153038380 -0.0000000017
   0.0308400945 0.0000000014
   -0.0061669180 0.0000000128
   0.0344133868 -0.0000000005
   0.4061633780 0.0000000000
   0.0941719022 0.0000000000
   0.0000000000 0.0000000000
   0.5003352802 0.0000000000
   0.4003838557 0.0000000000
   0.0999514245 0.0000000000
   0.1063444094 0.0000000000
   0.3939908708 0.0000000000
   0.4766029066 0.0000000000
   0.0237323736 0.0000000000
   0.1007150782 0.0000000000
   0.3996202019 0.0000000000
   0.3924422892 0.0000000000
   0.1078929910 0.0000000000
   0.1094415660 0.0000000000
   0.3908937142 0.0000000000
   0.3978256445 0.0000000000
   0.1025096357 0.0000000000
   0.1007150728 0.0000000000
   0.3996202074 0.0000000000
   0.3918006961 0.0000000000
   0.1085345841 0.0000000000
   0.1094415637 0.0000000000
   0.3908937164 0.0000000000
EOF
paste output/zvo_cisajscktalt_eigen2.dat reference.dat > paste7.dat
diff=`awk '
BEGIN{diff=0.0} 
{diff+=sqrt(($9-$11)*($9-$11)+($10-$12)*($10-$12))}
END{printf "%8.6f", diff/NR}
' paste7.dat`
echo "Diff output/zvo_cisajscktalt_eigen2.dat : " ${diff}
test "${diff}" = "0.000000"

exit $?
