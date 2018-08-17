
import os,sys,re


  
fout = open('manifest.txt', 'w') # 合并内容到该文件 

nfile=0
fout.write('sample-id,absolute-filepath,direction\n')
for root, dirs, files in os.walk(sys.argv[1]): 
	if len(dirs)==0: 
		for fl in files: 
			info = "$PWD/%s/%s" % (root,fl) 
			
			if info[-3:] =='.gz': # 只将后缀名为py的文件内容合并 
				nfile+=1
				sn = re.search('F.P.\.[\d]+',fl).group()
				if re.search("2.fq.gz$",fl):
					fout.write("%s,%s,reverse\n"%(sn,info))
				else:
					fout.write("%s,%s,forward\n"%(sn,info))
					



fout.close() 
emm = "    %s fastq file were writed" % (nfile) 
print(emm)


