##my_dict= {"my_key": ('1','2',3),"second_key": '4', "third_key": "5"}
##print[number for number, name in my_dict.items() if name != 0

#my_dict = {"one": 1,"two":2,"three":3,"four":4}

#for key in my_dict.values():
#    print("Key : {} , Value : {}".format(key,my_dict[key]))


ortho_file = open("Ssid_Orthogroups.txt", "r").read().split("\r\n")

first_dict = {}
for line in ortho_file[1:]:
	orthogroup = line.split("\t")[0]
	contig_list = line.split("\t")[1:]
	first_dict[orthogroup] = []
	if len(contig_list[0]) < 1:
		first_dict[orthogroup].append("0")
	else:
		for contig in contig_list:
			if len(contig) > 1:
				first_dict[orthogroup].append(contig)

outfile = open("Ssid_gene_to_orthogroup.csv", "w")
outfile.write("Orthogroup" + "," + "Contig" + "\n")
for key, value in first_dict.items():
	for contig in value:
		outfile.write(key + "," + contig + "\n")

