LVB
---
***Phylogenetic Software Package***

LVB is a maximum-parsimony (MP) phylogenetic inference tool currently under development by the [*Barker Lab Group*](https://www.ed.ac.uk/profile/daniel-barker) at the University of Edinburgh. LVB utilises a simulated annealing algorithm to aid searches within tree space.

Current Developers:\
[*Joseph Guscott*](https://github.com/josephguscott)\
[*Daniel Barker*](https://www.ed.ac.uk/profile/daniel-barker)

[![](https://img.shields.io/badge/Build-Passing-brightgreen)](https://github.com/phylolvb/lvb/releases/tag/3.5)
[![](https://img.shields.io/badge/Core%20Tests-Passing-brightgreen)]()
[![](https://img.shields.io/badge/Current%20Release-3.5-blue)](https://github.com/phylolvb/lvb/releases/tag/3.5)
[![](https://img.shields.io/badge/Release%20Date-02%2F2019-blue)](https://github.com/phylolvb/lvb/releases/tag/3.5)\
[![](https://img.shields.io/badge/DOI%3A-https%3A%2F%2Fdoi.org%2F10.1093%2Fbioinformatics%2Fbtg402-blue)](https://doi.org/10.1093/bioinformatics/btg402)

---

Installation
---

***Downloading source code:***

~~~~
git clone https://github.com/phylolvb/lvb 
~~~~

***Building LVB:***

Serial array version:
~~~~
make
~~~~

Hash version:
~~~~
make -f Makefile.HASH
~~~~

MapReduce version:
~~~~
make -f Makefile.MAPREDUCE
~~~~

---

Version History
---
***4.0***\
Currently under development

***3.5***\
*February 2019*\
This update enabled LVB to carry out TBR branch-swapping, in addition to the already previously implemented NNI and SPR. It also saw some development to algorithms that tried to identify the most advantageous combinations in which to run the different branch-swapping procedures. Other changes included refinement of the simulated annealing starting temperature, more detailed output statistics, and removal of the MSF input format, bootstrapping and weighting.   

---

Publications
---

- Strobl, M.A. and Barker, D., 2016. On simulated annealing phase transitions in phylogeny reconstruction. Molecular Phylogenetics and Evolution, 101, pp.46-55.\
DOI: [*10.1016/j.ympev.2016.05.001*](https://www.sciencedirect.com/science/article/pii/S1055790316300823?via%3Dihub)


- Barker, D., 2004. LVB: parsimony and simulated annealing in the search for phylogenetic trees. Bioinformatics, 20, pp.274-275.\
DOI: [*10.1093/bioinformatics/btg402*](https://academic.oup.com/bioinformatics/article/20/2/274/204936)

---
