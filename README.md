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