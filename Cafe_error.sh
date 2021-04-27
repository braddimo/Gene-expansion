#!/home/bad7666/miniconda2/bin/cafe
tree ((past:205.27,ssid:205.27):57.86,((ofav:106.13,mcav:106.13):99.14,cnat:205.27):57.86);
load -i Orthogroups.GeneCount_no_zeros.txt -t 10 -l reports/run6_caferror_files/cafe_final_log.txt
errormodel -model reports/run6_caferror_files/cafe_errormodel_1.220703125e-05.txt -sp past
errormodel -model reports/run6_caferror_files/cafe_errormodel_1.220703125e-05.txt -sp ofav
errormodel -model reports/run6_caferror_files/cafe_errormodel_1.220703125e-05.txt -sp cnat
errormodel -model reports/run6_caferror_files/cafe_errormodel_1.220703125e-05.txt -sp ssid
errormodel -model reports/run6_caferror_files/cafe_errormodel_1.220703125e-05.txt -sp mcav
lambda -s -t (1,((1,1)1,(1,1)1)1);
report reports/run6_caferror_files/cafe_final_report
