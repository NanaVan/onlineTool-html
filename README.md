# onlineTool-html
A html for [online tool](https://online-tool-html.vercel.app/) of beam experiment
Constant values: CODATA 2018
Nuclear data: NUBASE2020
Database: SQLite([SQL.js](https://github.com/sql-js/sql.js))

## Ion Number Estimation
A small tool for calculating transfer efficiency based on information of the primary beam and the secondary beam.

## Harmonic Calculation
A tool for calculatiing the location and harmonic of the ion peaks in the span of the observed frequency domain.

## Storage Lifetime Estimation
A tool for estimating storage 1/e lifetime of the ion in storage rings.

Beam Loss Mechanism:
* Radioactive Electron Capture (REC) with E-COOLing: 
	* REC@ H.Poth. Electron cooling: Theory, experiment, application. Physics Reports. 196 (1990): 135-297
* Elastic Scattering (ES) with residual gas: 
	* single-scattering@ H.Poth. Electron cooling: Theory, experiment, application. Physics Reports. 196 (1990): 135-297
	* multiple-scattering@ B.Franzke. Interaction of stored ion beams with the residual gas. CAS-CERN Accelerator School: 4th advanced accelerator physics course. (1992): 100-119
* Electron Capture (EC) with residual gas:
	* single-EC@ I.S. Dmitriev, et al. On the target thickness to attain equilibrium charge distribution in a beam of fast ions. Nucl. Instr. Meth. Phys. Res. B. 14. (1986): 515-526
	* EC@ A.S. Schlacher, et al. Electron capture for fast highly charged ions in gas targets: an empirical scaling rule. Phys. Rev. A. 27. (1983) 11: 3372
	* single-EC@ B. Franzke. Vacuum requirements for heavy ion synchrotrons. IEEE Trans. Nucl. Sci. NS-28 (1981): 2116-2118
* Electron Stripping (EL) with residual gas:
	* EL@ D. Habs, et al. First experiments with the Heidelberg test storage ring TSR. Nucl. Instr. Meth. Phys. Res. B. 43. (1989): 390-410
	* single-EL@ B. Franzke. Vacuum requirements for heavy ion synchrotrons. IEEE Trans. Nucl. Sci. NS-28 (1981): 2116-2118
* Nuclear Reaction (NR) with residual gas:
	* NR @ X.Y. Zhang. Study of Lifetime of Light Ion Beam stored in CSRm for Internal Target Experiment. (2005) MSc thesis
