# in the downloaded repertoire files, separate the HIP files to one folder
# and Keck to another

import os


directory = '../../data/emerson_download'
file_list = []
for entry in os.scandir(directory):
    file_list += [entry.path]

file_list.sort()

for file_path in file_list:
    file_name = file_path.split("/")[-1]
    if file_name.split(".")[-1] == "tsv":
        if file_name[:1] == "P":
            cmd = 'mv ../../data/emerson_download/'+file_name+\
                  ' ../../data/emerson-2017-natgen/HIP_folder/'+file_name
            os.system(cmd)
        else:
            cmd = 'mv ../../data/emerson_download/'+file_name+\
                  ' ../../data/emerson-2017-natgen/Keck_folder/'+file_name
            os.system(cmd)
