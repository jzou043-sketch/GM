GM extra traits ready for server rerun
======================================

Mouse full GM input files are provided as both CSV and XLSX.
Pig full GM input files are provided as CSV because the pig genotype matrix has 52,843 SNP columns, exceeding Excel's 16,384-column worksheet limit.

Use CSV files on the server.

Mouse traits:
  Obesity.BodyLength
  Biochem.Glucose

Pig traits:
  t1
  t2

Main server inputs:
  mouse/mouse_obesity_bodylength_gm_dat.csv
  mouse/mouse_biochem_glucose_gm_dat.csv
  pig/PIC_t1_gm_dat.csv
  pig/PIC_t2_gm_dat.csv
