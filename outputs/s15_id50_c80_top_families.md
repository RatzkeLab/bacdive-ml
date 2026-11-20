# Top-weighted protein families per label (logreg on all QC strains)

C=0.003; Swiss-Prot best hit (e<=1e-5); frac = share of +/- strains carrying it.


## l-arabinose|growth/assimilation

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.073 | PpSQ1_00400 | 6-deoxy-6-sulfo-D-gluconate dehydratase | 56 | 0.55 | 0.21 | 1836 |
| +0.069 | ytfR | Galactofuranose transporter ATP-binding protein YtfR | 52 | 0.38 | 0.06 | 1559 |
| +0.065 | yjfF | Inner membrane ABC transporter permease protein YjfF | 48 | 0.37 | 0.06 | 1511 |
| +0.065 | kdaD | 2-dehydro-3-deoxy-L-arabinonate dehydratase | 62 | 0.36 | 0.10 | 810 |
| +0.062 | ytfT | Galactofuranose transporter permease protein YtfT | 60 | 0.36 | 0.05 | 1385 |
| +0.056 | xylH | Xylose transport system permease protein XylH | 42 | 0.32 | 0.07 | 1391 |
| +0.056 | araA | L-arabinose isomerase | 63 | 0.26 | 0.08 | 2120 |
| +0.052 | xylD | D-xylonate dehydratase | 67 | 0.40 | 0.16 | 1043 |
| +0.050 | sbpA | Multiple sugar-binding periplasmic protein SbpA | 66 | 0.32 | 0.07 | 1401 |
| +0.045 | araA | L-arabinose 1-dehydrogenase (NAD(P)(+)) | 61 | 0.36 | 0.08 | 743 |
| +0.039 | ytfQ | Galactofuranose-binding protein YtfQ | 52 | 0.21 | 0.03 | 757 |
| +0.039 | araB | Ribulokinase | 51 | 0.14 | 0.04 | 1024 |
| -0.025 | msrB | Peptide methionine sulfoxide reductase MsrB | 74 | 0.51 | 0.55 | 4294 |
| -0.022 | styD | Phenylacetaldehyde dehydrogenase | 59 | 0.23 | 0.23 | 1546 |
| -0.022 | prpC | 2-methylcitrate synthase | 76 | 0.18 | 0.29 | 1109 |
| -0.022 | kdsD | Arabinose 5-phosphate isomerase KdsD | 54 | 0.43 | 0.58 | 2984 |

## l-arabinose|acid/gas (fermentation)

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.137 | araD | L-ribulose-5-phosphate 4-epimerase | 78 | 0.69 | 0.24 | 2074 |
| +0.132 | araA | L-arabinose isomerase | 63 | 0.67 | 0.09 | 2120 |
| +0.105 | araA | L-arabinose isomerase | 68 | 0.16 | 0.02 | 293 |
| +0.087 | Abf43 | Alpha-L-arabinofuranosidase 43 | 61 | 0.42 | 0.03 | 779 |
| +0.086 |  | (no hit) |  | 0.16 | 0.05 | 482 |
| +0.076 | xylB | Xylulose kinase | 22 | 0.17 | 0.02 | 423 |
| +0.074 | xylA2 | Xylose isomerase 2 | 66 | 0.30 | 0.13 | 1287 |
| +0.072 | xylT | D-xylose transporter | 43 | 0.10 | 0.01 | 157 |
| +0.061 | araB | Ribulokinase | 51 | 0.44 | 0.03 | 1024 |
| +0.054 | xylH | Xylose transport system permease protein XylH | 42 | 0.17 | 0.05 | 1391 |
| +0.053 | sbpA | Multiple sugar-binding periplasmic protein SbpA | 66 | 0.17 | 0.05 | 1401 |
| +0.051 | gpr | L-glyceraldehyde 3-phosphate reductase | 61 | 0.79 | 0.32 | 4083 |
| -0.034 | rplN | Large ribosomal subunit protein uL14 | 78 | 0.72 | 0.65 | 4987 |
| -0.031 | rplT | Large ribosomal subunit protein bL20 | 78 | 0.71 | 0.61 | 5741 |
| -0.026 | clpX | ATP-dependent Clp protease ATP-binding subunit ClpX | 55 | 0.96 | 0.97 | 8823 |
| -0.025 | rpsO | Small ribosomal subunit protein uS15 | 72 | 0.79 | 0.87 | 7062 |

## nitrate|reduction/e- acceptor

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.401 | narH | Nitrate reductase beta chain | 64 | 0.57 | 0.08 | 2682 |
| +0.382 | narG | Nitrate reductase alpha subunit | 67 | 0.55 | 0.08 | 2583 |
| +0.198 | narV | Respiratory nitrate reductase 2 gamma chain | 53 | 0.31 | 0.02 | 1344 |
| +0.140 | narX | Nitrate reductase-like protein NarX | 58 | 0.17 | 0.05 | 939 |
| +0.119 | narK | Nitrate/nitrite antiporter NarK | 54 | 0.22 | 0.01 | 1068 |
| +0.086 | narK1 | Probable nitrate/nitrite antiporter NarK1 | 36 | 0.13 | 0.01 | 441 |
| +0.078 | narP | Nitrate/nitrite response regulator protein NarP | 71 | 0.25 | 0.02 | 1235 |
| +0.078 | napA | Periplasmic nitrate reductase | 84 | 0.17 | 0.02 | 857 |
| +0.077 | ubiU | Ubiquinone biosynthesis hydroxylase UbiU | 57 | 0.33 | 0.05 | 1626 |
| +0.057 | mog | Molybdopterin adenylyltransferase | 59 | 0.31 | 0.08 | 1763 |
| +0.057 | dmsB | Anaerobic dimethyl sulfoxide reductase chain B | 34 | 0.03 | 0.01 | 171 |
| +0.056 | napC | Cytochrome c-type protein NapC | 54 | 0.10 | 0.01 | 435 |
| -0.044 | kdpA | Potassium-transporting ATPase potassium-binding subunit | 72 | 0.33 | 0.28 | 3100 |
| -0.038 | phrB | (6-4) photolyase | 54 | 0.09 | 0.12 | 837 |
| -0.037 | kdpB | Potassium-transporting ATPase ATP-binding subunit | 70 | 0.39 | 0.38 | 3894 |
| -0.036 | TT_C0602 | Phenylacetate-coenzyme A ligase | 51 | 0.27 | 0.27 | 2675 |

## nitrite|reduction/e- acceptor

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.059 | nosZ | Nitrous-oxide reductase | 73 | 0.34 | 0.03 | 589 |
| +0.053 | ubiU | Ubiquinone biosynthesis hydroxylase UbiU | 57 | 0.62 | 0.11 | 1626 |
| +0.045 | norB | Nitric oxide reductase subunit B | 29 | 0.16 | 0.03 | 392 |
| +0.040 | norB | Nitric oxide reductase subunit B | 69 | 0.27 | 0.03 | 570 |
| +0.038 | nirF | Protein NirF | 60 | 0.24 | 0.01 | 316 |
| +0.037 | nirJ | Pre-heme d1 synthase | 56 | 0.24 | 0.02 | 318 |
| +0.036 | nirN | Dihydro-heme d1 dehydrogenase | 56 | 0.22 | 0.01 | 303 |
| +0.035 | nirG | Siroheme decarboxylase NirG subunit | 50 | 0.23 | 0.01 | 299 |
| +0.033 | nirH | Siroheme decarboxylase NirH subunit | 54 | 0.22 | 0.01 | 293 |
| +0.032 | aniA | Copper-containing nitrite reductase | 66 | 0.11 | 0.02 | 268 |
| +0.032 | nirS | Nitrite reductase | 75 | 0.21 | 0.01 | 281 |
| +0.032 | narH | Nitrate reductase beta chain | 64 | 0.53 | 0.28 | 2682 |
| -0.036 | dnaK | Chaperone protein DnaK | 82 | 0.84 | 0.70 | 6506 |
| -0.030 | ychF | Ribosome-binding ATPase YchF | 73 | 0.89 | 0.89 | 7982 |
| -0.028 | trpB | Tryptophan synthase beta chain | 78 | 0.82 | 0.77 | 7086 |
| -0.028 | pstB | Phosphate import ATP-binding protein PstB | 80 | 0.71 | 0.85 | 7620 |

## sulfate|reduction/e- acceptor

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.008 | aprA | Adenylylsulfate reductase subunit alpha | 58 | 1.00 | 0.09 | 57 |
| +0.008 | hdrA2 | CoB--CoM heterodisulfide reductase iron-sulfur subunit A 2 | 42 | 1.00 | 0.09 | 56 |
| +0.007 | tusE | Sulfurtransferase TusE | 58 | 1.00 | 0.12 | 171 |
| +0.007 | dsrA | Sulfite reductase, dissimilatory-type subunit alpha | 59 | 1.00 | 0.12 | 58 |
| +0.007 | dsrB | Sulfite reductase, dissimilatory-type subunit beta | 59 | 1.00 | 0.12 | 58 |
| +0.007 |  | (no hit) |  | 0.61 | 0.03 | 33 |
| +0.007 | Y2-aiiA | N-acyl homoserine lactonase | 30 | 0.83 | 0.00 | 25 |
| +0.006 | cydB | Cytochrome bd ubiquinol oxidase subunit 2 | 39 | 0.89 | 0.03 | 41 |
| +0.006 |  | (no hit) |  | 0.78 | 0.00 | 29 |
| +0.006 | atpC | ATP synthase epsilon chain | 71 | 0.83 | 0.00 | 35 |
| +0.006 | rlmE | Ribosomal RNA large subunit methyltransferase E | 71 | 0.83 | 0.00 | 33 |
| +0.006 |  | (no hit) |  | 0.83 | 0.00 | 32 |
| -0.011 | purE | N5-carboxyaminoimidazole ribonucleotide mutase | 83 | 0.11 | 0.78 | 7197 |
| -0.010 | eno | Enolase | 93 | 0.44 | 0.94 | 5141 |
| -0.008 | hisF | Imidazole glycerol phosphate synthase subunit HisF | 69 | 0.22 | 0.78 | 6453 |
| -0.008 | aq_141 | Uncharacterized HIT-like protein aq_141 | 56 | 0.00 | 0.53 | 2150 |

## decanoate|growth/assimilation

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.027 | dkgB | Methylglyoxal reductase DkgB | 61 | 0.58 | 0.09 | 1216 |
| +0.023 | At3g03980 | NADPH-dependent aldehyde reductase-like protein, chloroplastic | 38 | 0.34 | 0.09 | 873 |
| +0.020 | crgA | HTH-type transcriptional regulator CrgA | 50 | 0.64 | 0.09 | 831 |
| +0.020 |  | (no hit) |  | 0.42 | 0.20 | 1279 |
| +0.019 | yhjC | Uncharacterized HTH-type transcriptional regulator YhjC | 44 | 0.23 | 0.04 | 335 |
| +0.019 | cueR | HTH-type transcriptional regulator CueR | 32 | 0.21 | 0.08 | 482 |
| +0.019 | cdhR | HTH-type transcriptional regulator CdhR | 30 | 0.42 | 0.13 | 1090 |
| +0.018 | cydA | Cytochrome bd-I ubiquinol oxidase subunit 1 | 56 | 0.55 | 0.30 | 2105 |
| +0.018 | baeR | Transcriptional regulatory protein BaeR | 47 | 0.22 | 0.07 | 522 |
| +0.018 | astB | N-succinylarginine dihydrolase | 84 | 0.66 | 0.13 | 1612 |
| +0.017 | TTHA0550 | Putative enoyl-CoA hydratase | 34 | 0.84 | 0.23 | 1523 |
| +0.017 | mmsA | Methylmalonate-semialdehyde dehydrogenase [acylating] | 71 | 0.80 | 0.15 | 1222 |
| -0.023 | ycf16 | Probable ATP-dependent transporter ycf16 | 61 | 0.17 | 0.74 | 6790 |
| -0.023 | pyrG | CTP synthase | 84 | 0.51 | 0.88 | 7437 |
| -0.023 | rpsO | Small ribosomal subunit protein uS15 | 72 | 0.81 | 0.81 | 7062 |
| -0.022 | ung | Uracil-DNA glycosylase | 65 | 0.12 | 0.30 | 3090 |

## d-xylose|acid/gas (fermentation)

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.155 | xylA2 | Xylose isomerase 2 | 66 | 0.49 | 0.07 | 1287 |
| +0.068 | xylB | Xylulose kinase | 70 | 0.23 | 0.04 | 416 |
| +0.065 | araA | L-arabinose isomerase | 63 | 0.50 | 0.17 | 2120 |
| +0.053 | araD | L-ribulose-5-phosphate 4-epimerase | 78 | 0.49 | 0.18 | 2074 |
| +0.051 | Abf43 | Alpha-L-arabinofuranosidase 43 | 61 | 0.29 | 0.04 | 779 |
| +0.050 | xylA | Xylose isomerase | 65 | 0.24 | 0.03 | 1505 |
| +0.048 | xynB | Beta-xylosidase | 48 | 0.15 | 0.02 | 285 |
| +0.044 | xylH | Xylose transport system permease protein XylH | 53 | 0.19 | 0.01 | 889 |
| +0.041 | xylR | Xylose operon regulatory protein | 99 | 0.16 | 0.01 | 634 |
| +0.041 |  | (no hit) |  | 0.17 | 0.06 | 482 |
| +0.040 | xylF | D-xylose-binding periplasmic protein | 59 | 0.20 | 0.02 | 953 |
| +0.040 | xylE | D-xylose-proton symporter | 54 | 0.13 | 0.01 | 385 |
| -0.030 | tcyB | L-cystine transport system permease protein TcyB | 57 | 0.28 | 0.22 | 2111 |
| -0.028 | rplW | Large ribosomal subunit protein uL23 | 55 | 0.21 | 0.44 | 2143 |
| -0.027 | rpsQ | Small ribosomal subunit protein uS17 | 85 | 0.63 | 0.67 | 5299 |
| -0.026 | argG | Argininosuccinate synthase | 79 | 0.27 | 0.46 | 2141 |

## l-rhamnose|acid/gas (fermentation)

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.089 | rhaM | L-rhamnose mutarotase | 68 | 0.67 | 0.04 | 769 |
| +0.081 | rhaA | L-rhamnose isomerase | 68 | 0.29 | 0.05 | 590 |
| +0.054 | rhaB | Rhamnulokinase | 94 | 0.56 | 0.01 | 511 |
| +0.050 | rhaT | L-rhamnose-proton symporter | 81 | 0.53 | 0.01 | 481 |
| +0.042 | yesR | Unsaturated rhamnogalacturonyl hydrolase YesR | 34 | 0.42 | 0.02 | 519 |
| +0.041 | rhaD | Rhamnulose-1-phosphate aldolase | 61 | 0.54 | 0.01 | 494 |
| +0.040 | rhaD | Rhamnulose-1-phosphate aldolase | 69 | 0.11 | 0.04 | 239 |
| +0.039 | iolF | Minor myo-inositol transporter IolF | 55 | 0.09 | 0.02 | 229 |
| +0.036 | rhaB | Rhamnulokinase | 63 | 0.08 | 0.03 | 158 |
| +0.034 | kduD | 2-dehydro-3-deoxy-D-gluconate 5-dehydrogenase | 59 | 0.69 | 0.22 | 2332 |
| +0.032 | araA | L-arabinose isomerase | 63 | 0.72 | 0.21 | 2120 |
| +0.030 | fucO | Lactaldehyde reductase | 62 | 0.49 | 0.07 | 757 |
| -0.032 | tuf1 | Elongation factor Tu | 95 | 0.87 | 0.94 | 8277 |
| -0.031 | fhs | Formate--tetrahydrofolate ligase | 89 | 0.07 | 0.30 | 1430 |
| -0.028 | infA | Translation initiation factor IF-1 | 94 | 0.97 | 0.97 | 8612 |
| -0.028 | rpsL | Small ribosomal subunit protein uS12 | 98 | 0.91 | 0.90 | 7916 |

## d-galacturonic acid|acid/gas (fermentation)

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.010 | dgcM | Diguanylate cyclase DgcM | 62 | 0.75 | 0.00 | 298 |
| +0.010 | silS | Probable sensor kinase SilS | 61 | 0.91 | 0.03 | 328 |
| +0.010 | ytfQ | Galactofuranose-binding protein YtfQ | 90 | 0.93 | 0.10 | 792 |
| +0.008 | exuT | Hexuronate transporter | 66 | 0.95 | 0.12 | 679 |
| +0.008 | uxaB | Altronate oxidoreductase | 90 | 0.98 | 0.14 | 512 |
| +0.008 | pptA | Tautomerase PptA | 73 | 0.77 | 0.03 | 323 |
| +0.007 | yecJ | Uncharacterized protein YecJ | 80 | 0.95 | 0.12 | 428 |
| +0.007 | ytfR | Galactofuranose transporter ATP-binding protein YtfR | 52 | 0.98 | 0.17 | 1559 |
| +0.007 | ytfT | Galactofuranose transporter permease protein YtfT | 60 | 0.98 | 0.17 | 1385 |
| +0.007 | yjfF | Inner membrane ABC transporter permease protein YjfF | 48 | 0.98 | 0.17 | 1511 |
| +0.007 | yesO | Putative ABC transporter substrate-binding protein YesO | 27 | 0.45 | 0.05 | 236 |
| +0.007 | yjcB | Uncharacterized protein YjcB | 64 | 0.91 | 0.12 | 402 |
| -0.009 | rplM | Large ribosomal subunit protein uL13 | 73 | 0.34 | 0.80 | 7508 |
| -0.009 | pstB | Phosphate import ATP-binding protein PstB | 80 | 0.14 | 0.83 | 7620 |
| -0.009 | hemB | Delta-aminolevulinic acid dehydratase | 87 | 0.07 | 0.63 | 1849 |
| -0.007 | rplE | Large ribosomal subunit protein uL5 | 86 | 0.98 | 1.00 | 8106 |

## tryptophan|growth/assimilation

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.110 | tnaA | Tryptophanase | 64 | 0.72 | 0.03 | 636 |
| +0.047 | TVAG_054490 | Probable beta-eliminating lyase | 62 | 0.15 | 0.01 | 206 |
| +0.035 | atoD | Acetate CoA-transferase subunit alpha | 54 | 0.34 | 0.06 | 744 |
| +0.028 | murP | PTS system N-acetylmuramic acid-specific EIIBC component | 91 | 0.37 | 0.04 | 398 |
| +0.025 | uroc1 | Uurocanate hydratase | 48 | 0.18 | 0.05 | 841 |
| +0.024 | rbsA | Ribose import ATP-binding protein RbsA | 76 | 0.65 | 0.11 | 908 |
| +0.024 | pcaJ | 3-oxoadipate CoA-transferase subunit B | 52 | 0.47 | 0.58 | 5254 |
| +0.021 | tnaA1 | Tryptophanase 1 | 35 | 0.06 | 0.00 | 29 |
| +0.020 | mdeA | L-methionine gamma-lyase | 72 | 0.11 | 0.05 | 405 |
| +0.020 | mtr | Tryptophan-specific transport protein | 44 | 0.22 | 0.02 | 183 |
| +0.019 | Swol_1936 | Crotonyl-CoA hydratase | 56 | 0.11 | 0.04 | 339 |
| +0.018 | rbsB | Ribose import binding protein RbsB | 59 | 0.66 | 0.14 | 1273 |
| -0.035 | pyrG | CTP synthase | 84 | 0.30 | 0.75 | 7437 |
| -0.034 | atpD | ATP synthase subunit beta | 77 | 0.62 | 0.77 | 7158 |
| -0.034 | rmlC | dTDP-4-dehydrorhamnose 3,5-epimerase | 73 | 0.38 | 0.39 | 3006 |
| -0.032 | rpsM | Small ribosomal subunit protein uS13 | 94 | 0.96 | 0.87 | 8146 |

## urea|hydrolysis/degradation

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.177 | ureC | Urease subunit alpha | 95 | 0.82 | 0.30 | 3656 |
| +0.175 | ureA | Urease subunit gamma | 89 | 0.83 | 0.33 | 3712 |
| +0.139 | ureG | Urease accessory protein UreG | 84 | 0.45 | 0.17 | 2056 |
| +0.103 | ureB | Urease subunit beta | 75 | 0.60 | 0.28 | 2913 |
| +0.032 |  | (no hit) |  | 0.17 | 0.04 | 535 |
| +0.031 | ureB | Urease subunit beta | 62 | 0.05 | 0.02 | 321 |
| +0.027 | kdpA | Potassium-transporting ATPase potassium-binding subunit | 72 | 0.61 | 0.32 | 3100 |
| +0.026 | MJ0412 | Uncharacterized ABC transporter ATP-binding protein MJ0412 | 60 | 0.46 | 0.23 | 2589 |
| +0.024 | ramA | HTH-type transcriptional activator RamA | 44 | 0.20 | 0.04 | 424 |
| +0.024 | paoB | Aldehyde oxidoreductase FAD-binding subunit PaoB | 59 | 0.41 | 0.19 | 1995 |
| +0.023 | yurT | Uncharacterized protein YurT | 28 | 0.32 | 0.10 | 1155 |
| +0.023 | ureE | Urease accessory protein UreE | 68 | 0.03 | 0.00 | 49 |
| -0.037 | rplO | Large ribosomal subunit protein uL15 | 56 | 0.60 | 0.53 | 5059 |
| -0.035 | clpP | ATP-dependent Clp protease proteolytic subunit | 96 | 0.55 | 0.82 | 7239 |
| -0.033 | eno | Enolase | 93 | 0.38 | 0.53 | 5141 |
| -0.033 | clpB | Chaperone protein ClpB | 84 | 0.56 | 0.78 | 6768 |

## citrate|growth/assimilation

| coef | gene | protein | %id | frac in + | frac in - | n genomes |
|---|---|---|---|---|---|---|
| +0.046 |  | Uncharacterized 52.8 kDa protein in TAR-I ttuC' 3'region | 39 | 0.31 | 0.06 | 1041 |
| +0.043 | tctD | Transcriptional regulatory protein TctD | 48 | 0.33 | 0.06 | 1149 |
| +0.038 |  | (no hit) |  | 0.27 | 0.05 | 880 |
| +0.030 | citN | Citrate transporter | 85 | 0.17 | 0.06 | 812 |
| +0.029 | iaaA | Isoaspartyl peptidase | 50 | 0.20 | 0.08 | 897 |
| +0.026 | pmrB | Sensor protein kinase PmrB | 32 | 0.12 | 0.01 | 249 |
| +0.026 | leuD | 3-isopropylmalate dehydratase small subunit | 70 | 0.13 | 0.09 | 777 |
| +0.025 |  | (no hit) |  | 0.17 | 0.15 | 1219 |
| +0.025 | epsL | Uncharacterized sugar transferase EpsL | 58 | 0.15 | 0.18 | 1499 |
| +0.024 |  | (no hit) |  | 0.15 | 0.12 | 954 |
| +0.024 | MIMI_R617 | N-glycosidase R617 | 42 | 0.20 | 0.13 | 1055 |
| +0.024 |  | (no hit) |  | 0.12 | 0.01 | 257 |
| -0.031 | clpP | ATP-dependent Clp protease proteolytic subunit | 96 | 0.72 | 0.77 | 7239 |
| -0.028 | pcl1 | Fe(2+)/Mn(2+) transporter pcl1 | 54 | 0.23 | 0.23 | 1880 |
| -0.027 | rpmH | Large ribosomal subunit protein bL34 | 89 | 0.23 | 0.27 | 2203 |
| -0.027 | vapI | Virulence-associated protein I | 57 | 0.04 | 0.07 | 512 |
