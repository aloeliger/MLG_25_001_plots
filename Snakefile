rule all:
   input:
      "outputs/AXOL1TL_v3_axo_style_score_plot.pdf",
      "outputs/AXOL1TL_v3_axo_style_score_plot.png",
      "outputs/CICADA_2024_axo_style_score_plot.pdf",
      "outputs/CICADA_2024_axo_style_score_plot.png",
      "outputs/1D_correlation_plot.pdf",
      "outputs/1D_correlation_plot.png",

rule axo_style_score_plots:
   input:
      "inputs/CICADA2024_CICADAScore_plot_info.pkl",
      "inputs/axol1tl_v3_AXOScore_plot_info.pkl",
      "make_axo_style_score_plots.py",
   output:
      "outputs/AXOL1TL_v3_axo_style_score_plot.pdf",
      "outputs/AXOL1TL_v3_axo_style_score_plot.png",
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
      "inputs/hist_result_plotA_plotB_plotC.pkl",
      "makeObjMultPlots.py",
   output:
      "outputs/L1Jet_mult.pdf",
      "outputs/L1Jet_mult.png",
      "outputs/L1EG_mult.pdf",
      "outputs/L1EG_mult.png",
      "outputs/L1Mu_mult.pdf",
      "outputs/L1Mu_mult.png",
   shell:
      "python3 makeObjMultPlots.py --object L1Jet --input inputs/hist_result_plotA_plotB_plotC.pkl --output outputs/L1Jet_mult",
      "python3 makeObjMultPlots.py --object L1EG  --input inputs/hist_result_plotA_plotB_plotC.pkl --output outputs/L1EG_mult",
      "python3 makeObjMultPlots.py --object L1Mu  --input inputs/hist_result_plotA_plotB_plotC.pkl --output outputs/L1Mu_mult",

rule l1_dist_plots:
   input:
      "inputs/hist_result_plotD_plotE.pkl",
      "makeL1DistPlot.py",
   output:
      "outputs/l1_ht_efficiency.pdf",
      "outputs/l1_ht_efficiency.png",
      "outputs/l1_met_efficiency.pdf",
      "outputs/l1_met_efficiency.png",
   shell:
      "python3 makeL1DistPlot.py --input inputs/hist_result_plotD_plotE.pkl --output outputs/l1_ht_efficiency --observable ht",
      "python3 makeL1DistPlot.py --input inputs/hist_result_plotD_plotE.pkl --output outputs/l1_met_efficiency --observable met",

rule ht_purity_plot:
   input:
      "inputs/hist_result_plotF.pkl",
      "makeHTPurityPlot.py",
   output:
      "outputs/l1_ht_purity.pdf",
      "outputs/l1_ht_purity.png",
   shell:
      "python3 makeHTPurityPlot.py --input inputs/hist_result_plotF.pkl --output outputs/l1_ht_purity"

rule dimuon_mass_plot:
   input:
      "inputs/hist_result_plotG.pkl",
      "makeDimuonPlot.py",
   output:
      "outputs/dimuon_mass.pdf",
      "outputs/dimuon_mass.png",
   shell:
      "python3 makeDimuonPlot.py --input inputs/hist_result_plotG.pkl --output outputs/dimuon_mass"