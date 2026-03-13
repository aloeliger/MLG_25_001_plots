rule all:
   input:
      "outputs/AXOL1TL_v4_axo_style_score_plot.pdf",
      "outputs/AXOL1TL_v4_axo_style_score_plot.png",
      "outputs/CICADA_2024_axo_style_score_plot.pdf",
      "outputs/CICADA_2024_axo_style_score_plot.png",
      "outputs/1D_correlation_plot.pdf",
      "outputs/1D_correlation_plot.png",
      "outputs/L1Jet_mult.pdf",
      "outputs/L1Jet_mult.png",
      "outputs/L1EG_mult.pdf",
      "outputs/L1EG_mult.png",
      "outputs/L1Mu_mult.pdf",
      "outputs/L1Mu_mult.png",
      "outputs/l1_ht_dist.pdf",
      "outputs/l1_ht_dist.png",
      "outputs/l1_met_dist.pdf",
      "outputs/l1_met_dist.png",
      "outputs/l1_ht_purity.pdf",
      "outputs/l1_ht_purity.png",
      "outputs/L1Jet_mult_nPV10.pdf",
      "outputs/L1Jet_mult_nPV10.png",
      "outputs/L1EG_mult_nPV10.pdf",
      "outputs/L1EG_mult_nPV10.png",
      "outputs/L1Mu_mult_nPV10.pdf",
      "outputs/L1Mu_mult_nPV10.png",
      "outputs/l1_ht_dist_nPV10.pdf",
      "outputs/l1_ht_dist_nPV10.png",
      "outputs/l1_met_dist_nPV10.pdf",
      "outputs/l1_met_dist_nPV10.png",
      "outputs/l1_ht_purity_nPV10.pdf",
      "outputs/l1_ht_purity_nPV10.png",
      "outputs/dimuon_mass.pdf",
      "outputs/dimuon_mass.png",
      "outputs/dimuon_mass_nPV10.pdf",
      "outputs/dimuon_mass_nPV10.png",

rule axo_style_score_plots:
   input:
      "inputs/CICADA2024_CICADAScore_plot_info.pkl",
      "inputs/axol1tl_v4_AXOScore_plot_info.pkl",
      "make_axo_style_score_plots.py",
   output:
      "outputs/AXOL1TL_v4_axo_style_score_plot.pdf",
      "outputs/AXOL1TL_v4_axo_style_score_plot.png",
      "outputs/CICADA_2024_axo_style_score_plot.pdf",
      "outputs/CICADA_2024_axo_style_score_plot.png",
   shell:
      "python3 make_axo_style_score_plots.py --output outputs/"

rule correlation_plots:
   input:
      "inputs/correlation_dict.pkl",
      "make_correlation_plots.py",
   output:
      "outputs/1D_correlation_plot.pdf",
      "outputs/1D_correlation_plot.png",
   shell:
      "python3 make_correlation_plots.py --input inputs/correlation_dict.pkl --output outputs/"

rule obj_mult_plots:
   input:
      "inputs/hists_plotA_plotB_plotC.root",
      "makeObjMultPlots.py",
   output:
      "outputs/L1Jet_mult.pdf",
      "outputs/L1Jet_mult.png",
      "outputs/L1EG_mult.pdf",
      "outputs/L1EG_mult.png",
      "outputs/L1Mu_mult.pdf",
      "outputs/L1Mu_mult.png",
   shell:
      "python3 makeObjMultPlots.py --object L1Jet --input inputs/hists_plotA_plotB_plotC.root --output outputs/L1Jet_mult && "
      "python3 makeObjMultPlots.py --object L1EG  --input inputs/hists_plotA_plotB_plotC.root --output outputs/L1EG_mult && "
      "python3 makeObjMultPlots.py --object L1Mu  --input inputs/hists_plotA_plotB_plotC.root --output outputs/L1Mu_mult"

rule obj_mult_nPV10:
   input:
      "inputs/hists_plotA_plotB_plotC_nPV10.root",
      "makeObjMultPlots.py",
   output:
      "outputs/L1Jet_mult_nPV10.pdf",
      "outputs/L1Jet_mult_nPV10.png",
      "outputs/L1EG_mult_nPV10.pdf",
      "outputs/L1EG_mult_nPV10.png",
      "outputs/L1Mu_mult_nPV10.pdf",
      "outputs/L1Mu_mult_nPV10.png",
   shell:
      "python3 makeObjMultPlots.py --object L1Jet --input inputs/hists_plotA_plotB_plotC_nPV10.root --output outputs/L1Jet_mult_nPV10 && "
      "python3 makeObjMultPlots.py --object L1EG  --input inputs/hists_plotA_plotB_plotC_nPV10.root --output outputs/L1EG_mult_nPV10 && "
      "python3 makeObjMultPlots.py --object L1Mu  --input inputs/hists_plotA_plotB_plotC_nPV10.root --output outputs/L1Mu_mult_nPV10"

rule l1_dist_plots:
   input:
      "inputs/hists_plotD_plotE.root",
      "makeL1DistPlot.py",
   output:
      "outputs/l1_ht_dist.pdf",
      "outputs/l1_ht_dist.png",
      "outputs/l1_met_dist.pdf",
      "outputs/l1_met_dist.png",
   shell:
      "python3 makeL1DistPlot.py --input inputs/hists_plotD_plotE.root --output outputs/l1_ht_dist --observable ht && "
      "python3 makeL1DistPlot.py --input inputs/hists_plotD_plotE.root --output outputs/l1_met_dist --observable met"

rule l1_dist_plots_nPV10:
   input:
      "inputs/hists_plotD_plotE_nPV10.root",
      "makeL1DistPlot.py",
   output:
      "outputs/l1_ht_dist_nPV10.pdf",
      "outputs/l1_ht_dist_nPV10.png",
      "outputs/l1_met_dist_nPV10.pdf",
      "outputs/l1_met_dist_nPV10.png",
   shell:
      "python3 makeL1DistPlot.py --input inputs/hists_plotD_plotE_nPV10.root --output outputs/l1_ht_dist_nPV10 --observable ht && "
      "python3 makeL1DistPlot.py --input inputs/hists_plotD_plotE_nPV10.root --output outputs/l1_met_dist_nPV10 --observable met"


rule ht_purity_plot:
   input:
      "inputs/hists_plotF.root",
      "makeHTPurityPlot.py",
   output:
      "outputs/l1_ht_purity.pdf",
      "outputs/l1_ht_purity.png",
   shell:
      "python3 makeHTPurityPlot.py --input inputs/hists_plotF.root --output outputs/l1_ht_purity"

rule ht_purity_plot_nPV10:
   input:
      "inputs/hists_plotF_nPV10.root",
      "makeHTPurityPlot.py",
   output:
      "outputs/l1_ht_purity_nPV10.pdf",
      "outputs/l1_ht_purity_nPV10.png",
   shell:
      "python3 makeHTPurityPlot.py --input inputs/hists_plotF_nPV10.root --output outputs/l1_ht_purity_nPV10"

rule dimuon_mass_plot:
   input:
      "inputs/hists_plotG.root",
      "makeDimuonPlot.py",
   output:
      "outputs/dimuon_mass.pdf",
      "outputs/dimuon_mass.png",
   shell:
      "python3 makeDimuonPlot.py --input inputs/hists_plotG.root --output outputs/dimuon_mass"

rule dimuon_mass_plot_nPV10:
   input:
      "inputs/hists_plotG_nPV10.root",
      "makeDimuonPlot.py",
   output:
      "outputs/dimuon_mass_nPV10.pdf",
      "outputs/dimuon_mass_nPV10.png",
   shell:
      "python3 makeDimuonPlot.py --input inputs/hists_plotG_nPV10.root --output outputs/dimuon_mass_nPV10"
