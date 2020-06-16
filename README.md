# mocktecarlo
Academic example written to validate part of the code used in [this publication](https://www.sciencedirect.com/science/article/pii/S037026931830337X?via%3Dihub)

Result were obtained using:
   - Cuba Library (4.2)
   - LHAPDF (6.1.6)
   - Fastjet (3.2.1)

The code can be readily compiled with
```
scons .
```

Although the code served its purposed it as left in github as it serves as an example on how to write a MC for HEP as it has several examples of all necessary pieces:

- specific phase space generator
- matrix elements for LO and NLO
- subtraction terms for NLO
- convolution with the PDFs
- cuts on the final state particles

Any questions, feel free to contact me!

If it served you in some way (let me know!) and please cite:
```
@article{Cruz-Martinez:2018rod,
    author = "Cruz-Martinez, J. and Gehrmann, T. and Glover, E.W.N. and Huss, A.",
    title = "{Second-order QCD effects in Higgs boson production through vector boson fusion}",
    eprint = "1802.02445",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "CERN-TH-2018-020, IPPP-18-8, ZU-TH-05-18",
    doi = "10.1016/j.physletb.2018.04.046",
    journal = "Phys. Lett. B",
    volume = "781",
    pages = "672--677",
    year = "2018"
}
```
