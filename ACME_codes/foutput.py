def foutput(content, output_path):
    content = str(content)
    fout =open(output_path,"a+")
    fout.write(content+"\n")
    fout.close()   