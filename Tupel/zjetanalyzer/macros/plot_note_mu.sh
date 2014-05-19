ismu=1
declare -a arr=(
mjjjj
n_jet_30_unfold
MZ
Zpt
Zrapidity
Zphi
ljet_rapidity
ljet_eta
sljet_rapidity
sljet_eta
nvtx
jetpt_1
jetpt_2
jetY_1
jetY_2
lpt
leta
lphi
)
declare -a title=(
"M(jjjj)"
"N_{jet}"
"M_{#mu#mu}[GeV]"
"P_{T}^{#mu#mu}[GeV]"
"Y_{#mu#mu}"
"#phi_{#mu#mu}"
"ljet_rapidity"
"ljet_eta"
"sljet_rapidity"
"sljet_eta"
"nvtx"
"P_{T}^{jet}[GeV]"
"P_{T}^{jet}[GeV]"
"#eta_{jet}"
"#eta_{jet}"
"P_{T}^{#mu}[GeV]"
"#eta_{#mu}"
"#phi_{#mu}"
)


for ((n =0 ; n<${#arr[@]} ;n++ ))
do

                echo "${arr[$n]}"
                echo "${title[$n]}"
                cat plot_note.C | sed "s/__HISTNAME__/${arr[$n]}/g
					s/__ISMU__/$ismu/g"> macro_tmp2.C 
if [ $n == 3 ] || [ $n == 0 ] || [ $n == 15 ] ; then
		cat macro_tmp2.C | sed "s/__NBIN__/4/g"> macro_tmp3.C
		cat macro_tmp3.C | sed "s/__TITLE__/${title[$n]}/g"> macro_tmp.C
else 

	                     cat macro_tmp2.C | sed "s/__NBIN__/1/g"> macro_tmp3.C
                             cat macro_tmp3.C | sed "s/__TITLE__/${title[$n]}/g"> macro_tmp.C		
       
fi
		root -l -q macro_tmp.C 

rm macro_tmp2.C
rm macro_tmp3.C
rm macro_tmp.C
done
