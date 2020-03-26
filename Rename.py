# coding=UTF-8
import os
import json

#打开json文件夹；
## 根据自己情况修改*号为自己的路经
file = open('*/TCGA/metadata.cart.2020-03-20.json', encoding='utf-8') 
#读取json文件夹；
json_precess = json.loads(str(file.read()))
#创建空字典；
dict ={}
for i in json_precess:
    print(i['file_name'])
    print(i['associated_entities'][0]['entity_submitter_id'])
    dict[str(i['file_name']).strip('.gz')] =i['associated_entities'][0]['entity_submitter_id']
print(dict)


#mainfest文件路径；
path ='*/TCGA/mainfest'
filelist = os.listdir(path)  
#mainfest子路径下的所有文件列表；
for file_one in filelist:
    file = path + '/' +file_one
    file_tar=file_one + '.FPKM.txt.gz'
    key=file_one +'.FPKM.txt'
    print(file_tar)
    print(file)
    list = os.listdir(file)[0]
    print(list)
    if '.gz' in list:

        olddir = os.path.join(path, file_tar)
        #print(olddir)
        filename =os.path.splitext(file_tar)[0]
        newdir = os.path.join(path,dict[key])
        os.rename(olddir,newdir)
