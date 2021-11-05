#/home/bad7666/miniconda3/bin/cafe
#version
#date
load -i Orthogroup_gene_counts_no_Zeros_revision.txt -t 5 -l reports/report1.txt
tree ((past:205.27,ssid:205.27):57.86,((ofav:106.13,mcav:106.13):99.14,cnat:205.27):57.86);
errormodel -model cafe_errormodel_0.1.txt -all
lambda -s -t (1,((1,1)1,(1,1)1)1);
report reports/cafe_report_1