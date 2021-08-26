1. 2018_08_21_Bayley_allPW
big PW-community experiment. What was grown: 7 species alone, all pairwise combinations of 7 species, and the complete 7-member community: all of these were grown both on CCA medium at pH 5 and at pH 7. 

Columns:
-spec: growing species
-cond: name of pairwise coculture partner species, or "alone" if monoculture or "community" if grown in complete community
-pH: medium starting pH (either pH 5 or pH 7)
-rep: replicate number (all time points/coculture conditions within a replicate were inoculated from the same cell stock dilution)
-day: number of days growing (day of inoculation is day 0).
-count: colonies counted on a plate
-dilution: log10 serial dilution plate (except for day zero where dilution plates were not 10-fold)
-cfus: back-calculated colony forming units (CFUs) in entire 96-well plug biofilm
-measured.pH: pH of sample homogenized into 1 mL 1X PBS (note these values are skewed by a buffer and do not accurately represent the pH of the biofilm/medium).
-class: kingdom of growing species (either bacteria or fungi).

--------------------
2. 2020_11_09_LeaveOneOut
community leave-one-out experiment done on CCA pH 5. What was grown: complete 7-member community, and few 6-member communities excluding one member (either Penicillium sp. st. JBC, Diutina catenulata 135E, or Staph xylosus BC10), and one 5-member community excluding two early deacidifiers ("eD"), so excluding both Diutina catenulata 135E & Staph xylosus BC10.

Columns:
-species: growing species
-condition: name of community condition: "comm" for complete 7-member community, "comm-JBC" for community excluding Penicillium sp. st. JBC, "comm-135E" when excluding Diutina, "comm-BC10" when excluding S. xylosus, "comm-eD" when excluding both Diutina and S. xylosus.
-media: (same info as "pH" column in dataset #1) starting condition of CCA (pH 5 for all of these)
-replicate, day, count, dilution: same as for dataset 1.
-pH: pH of sample, with biofilm measured directly (so these measurements DO represent actual pH of medium/biofilm).
-total.CFUs: (same info as "cfus" column in dataset #1) back-calculated colony forming units (CFUs) in entire 96-well plug biofilm

--------------------
3. 2020_11_17_pH7LOO
community leave-one-out experiment done on CCA pH 7. What was grown: 7 species alone, the complete 7-member community, and a 6-member community excluding Penicillium sp. st. JBC.

Columns:
-species: growing species
-condition: name of community condition ("comm" for complete 7-member community, "comm-JBC" for community excluding Penicillium JBC), or "alone" if monoculture.
-media: (same info as "pH" column in dataset #1) starting condition of CCA (pH 7 for all of these)
-replicate, day, count, dilution: same as for dataset 1.
-pH: (same info as "pH" column of dataset #2)
-total.CFUs: (same info as "cfus" column in dataset #1, and "total.CFUs" column in dataset #2)
-order of pH measure: safely ignore this, it was used only for me to assess any drift in pH probe measurements over time in use.