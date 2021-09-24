# banded_repeater_analysis
Estimating incompleteness due to banded repeaters 

This repository contains the data and analysis scripts used in the paper ([here](https://arxiv.org/abs/2108.04474)). 

Here is the layout for this repository:

* [`notebooks/`](https://github.com/KshitijAggarwal/banded_repeater_analysis/tree/main/notebooks): consists of the following analysis notebooks and scripts:

| Notebook | Description |
| ----------- | ----------- |
| [121102_energies.ipynb](https://github.com/KshitijAggarwal/banded_repeater_analysis/blob/main/notebooks/121102_energies.ipynb) | Analysis comparing FRB121102 bursts detected with FAST and Arecibo to simulated bursts.  |
| [fast_121102.ipynb](https://github.com/KshitijAggarwal/banded_repeater_analysis/blob/main/notebooks/bandlimited_spectra_vis.ipynb) | Comparing FAST FRB121102 burst energies calculated using two different methods. |        
| [spectra_param_pdf.ipynb](https://github.com/KshitijAggarwal/banded_repeater_analysis/blob/main/notebooks/spectra_param_pdf.ipynb) | Visualizing intrinsic and recovered distribution of spectra parameters. |
| [variable_power_law_slope.ipynb](https://github.com/KshitijAggarwal/banded_repeater_analysis/blob/main/notebooks/variable_power_law_slope.ipynb) | Effect of different intrinsic power law slopes on incompleteness. |
| [bandlimited_spectra_vis.ipynb](https://github.com/KshitijAggarwal/banded_repeater_analysis/blob/main/notebooks/bandlimited_spectra_vis.ipynb) | Visualizing burst spectra and its effect on detectability.|
| [power_laws_and_pdf.ipynb](https://github.com/KshitijAggarwal/banded_repeater_analysis/blob/main/notebooks/power_laws_and_pdf.ipynb) | Energy distributions of banded bursts detected at varying fluence thresholds. | 
| [subband_search.ipynb](https://github.com/KshitijAggarwal/banded_repeater_analysis/blob/main/notebooks/subband_search.ipynb) | Analysis of using subband searches on bandlimited repeating FRBs. |

* [`data/`](https://github.com/KshitijAggarwal/banded_repeater_analysis/tree/main/data): Consists of burst properties reported by [Aggarwal et al 2021](https://ui.adsabs.harvard.edu/abs/2021arXiv210705658A/abstract) and [Li et al 2021](https://ui.adsabs.harvard.edu/abs/2021arXiv210708205L/abstract). Please cite these two papers if you make use of data in this directory.  

* [`plotting.py`](https://github.com/KshitijAggarwal/banded_repeater_analysis/blob/main/plotting.py): Function to set figure dimensions to avoid scaling in LaTeX
* [`utils.py`](https://github.com/KshitijAggarwal/banded_repeater_analysis/blob/main/utils.py): Helper functions (fitting power laws, converting between fluence and energy, etc)
* [`simulate.py`](https://github.com/KshitijAggarwal/banded_repeater_analysis/blob/main/simulate.py): Functions to simulate properties of bandlimited repeater bursts and then run a mock search on them.

# Citation

Please cite the following paper if you use the data and analysis scripts present in this repository: 

```bash
@ARTICLE{2021arXiv210804474A,
       author = {{Aggarwal}, Kshitij},
        title = "{Observational effects of banded repeating FRBs}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2021,
        month = aug,
          eid = {arXiv:2108.04474},
        pages = {arXiv:2108.04474},
archivePrefix = {arXiv},
       eprint = {2108.04474},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021arXiv210804474A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
