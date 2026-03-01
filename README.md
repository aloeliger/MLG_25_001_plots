# Environment

```
cmsrel CMSSW_15_1_0_pre1
cd CMSSW_15_1_0_pre1/src/
git clone https://github.com/aloeliger/MLG_25_001_plots.git
cd MLG_15_001_plots/
python3 -m pip venv plot_env
source plot_env/bin/activate
python3 -m pip install -r requirements.txt
```

# Running all plots

## Panel/Figure 3 Plot A
`python makeObjMultPlots.py --object L1Jet --input inputs/hist_result_plotA_plotB_plotC.pkl --output outputs/mult_L1Jet_mult`

## Panel/Figure 3 Plot B
`python makeObjMultPlots.py --object L1EG  --input inputs/hist_result_plotA_plotB_plotC.pkl --output outputs/mult_L1EG_mult`

## Panel/Figure 3 Plot C
`python makeObjMultPlots.py --object L1Mu  --input inputs/hist_result_plotA_plotB_plotC.pkl --output outputs/mult_L1Mu_mult`

## Panel/Figure 3 Plot D
`python makeL1DistPlot.py --input inputs/hist_result_plotD_plotE.pkl --output outputs/l1_ht_efficiency --observable ht`

## Panel/Figure 3 Plot E
`python makeL1DistPlot.py --input inputs/hist_result_plotD_plotE.pkl --output outputs/l1_met_efficiency --observable met`

## Panel/Figure 3 Plot F
`python makeHTPurityPlot.py --input inputs/hist_result_plotF.pkl --output outputs/l1_ht_purity`

## Panel/Figure 3 Plot G
`python makeDimuonPlot.py --input inputs/hist_result_plotG.pkl --output outputs/dimuon_mass`