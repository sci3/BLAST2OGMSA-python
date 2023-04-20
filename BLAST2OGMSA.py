import sys

USAGE = "\nusage: python BLAST2OGMSA.py -method=[Gblocks|trimAl|BMGE|noisy] <file.aln> <seqdump.txt> <output.fasta>\n"

method = ''
aln = sys.argv[1]
seqdump = sys.argv[2]
out = sys.argv[3]

for paras in sys.argv:
    if '-help' in paras or '-h' in paras:
        print(USAGE)
        sys.exit()
    elif 'method' in paras:
        method = paras.split('=')[1]

if not aln:
    print(USAGE)
    print("Please provide the raw blast alignment file download from MSA viewer!\n")
    sys.exit()
if not seqdump:
    print(USAGE)
    print("Please provide the raw sequence file download from blastn results.\n")
    sys.exit()
if not out:
    print(USAGE)
    print("Please provide the name of output file.\n")
    sys.exit()
##############################################################################
list = []
with open(seqdump, "r") as file:
    for line in file:
        if line.startswith(">"):
            y = line[1:].split("|")
            if len(y) > 2:
                w = y[1].split()
            else:
                w = y[0].split()
            id = f"{w[0]}-{w[1]}-{w[2]}-{w[3]}"
            list.append(id)
############################################################################
seq = {}
seq_number = 0

with open(aln, "r") as infile:
    for line in infile:
        line = line.rstrip()
        if line.startswith(">"):
            seq_number = 0
            sid = line[1:]
            w = sid.split("|")
            sid = w[1]
        else:
            seq_number += 1
            compare_number2 = len(re.findall("\w", line))
            if sid not in seq:
                seq[sid] = [""] * seq_number + [line]
            else:
                compare_number1 = len(re.findall("\w", seq[sid][seq_number-1]))
                if compare_number2 >= compare_number1:
                    seq[sid][seq_number-1] = line

with open(aln + ".temp", "w") as outfile:
    for sid, sequences in seq.items():
        outfile.write(">" + sid + "\n")
        for outseq in sequences:
            outfile.write(outseq + "\n")
##############################################################################
import os
import glob

trimed = glob.glob("*.temp")
for trimed_file in trimed:
    if method == "Gblocks":
        os.system(f"./bin/Gblocks {trimed_file} out")
    if method == "trimAl":
        os.system(f"./bin/trimal -in {trimed_file} -out {trimed_file}-gb -fasta -htmlout {trimed_file}.html -automated1")
    if method == "BMGE":
        os.system(f"java -jar ./bin/BMGE.jar -i {trimed_file} -t DNA -s YES -of {trimed_file}-gb -oh {trimed_file}.html")
    if method == "noisy":
        os.system(f"./bin/noisy {trimed_file}")
    os.unlink(trimed_file)
########################################################################
import glob
import os

trimed = glob.glob("*.temp")

for trimed_file in trimed:
    if method == "Gblocks":
        gb_files = glob.glob("*.temp-gb")
        for gb_file in gb_files:
            delete = 0
            with open(gb_file, 'r') as gb, open(gb_file+'.out', 'w') as gbout:
                for line in gb:
                    if line.startswith('>'):
                        gbout.write(line)
                    elif line[0:10].isalnum():
                        delete += 1
                        gbout.write(line.replace(" ", ""))
            os.remove(gb_file)
            if delete == 0:
                os.remove(gb_file+'.out')
            temp_name2 = gb_file+'.out'
            temp_name2 = temp_name2.replace('temp-gb.out', 'fasta')
            os.rename(gb_file+'.out', temp_name2)
            
    elif method == "trimAl":
        gb_files = glob.glob("*.temp-gb")
        for gb_file in gb_files:
            delete = 0
            with open(gb_file, 'r') as gb, open(gb_file+'.out', 'w') as gbout:
                for line in gb:
                    if line.startswith('>'):
                        gbout.write(line)
                    elif line.replace('-', '').strip().isalpha():
                        delete += 1
                        gbout.write(line)
            os.remove(gb_file)
            if delete == 0:
                os.remove(gb_file+'.out')
            temp_name2 = gb_file.replace('temp-gb', 'fasta')
            os.rename(gb_file+'.out', temp_name2)
    
    elif method == "BMGE":
        gb_files = glob.glob("*.temp-gb")
        for gb_file in gb_files:
            delete = 0
            with open(gb_file, 'r') as gb, open(gb_file+'.out', 'w') as gbout:
                for line in gb:
                    if line.startswith('>'):
                        gbout.write(line)
                    elif line.replace('-', '').strip().isalpha():
                        delete += 1
                        gbout.write(line)
            os.remove(gb_file)
            if delete == 0:
                os.remove(gb_file+'.out')
            temp_name2 = gb_file.replace('temp-gb', 'fasta')
            os.rename(gb_file+'.out', temp_name2)
    
    elif method == "noisy":
        gb_files = glob.glob("*.fas")
        for gb_file in gb_files:
            delete = 0
            with open(gb_file, 'r') as gb, open(gb_file+'.out', 'w') as gbout:
                for line in gb:
                    if line.startswith('>'):
                        gbout.write(line.replace(" ", ""))
                    elif line.replace('-', '').strip().isalpha():
                        delete += 1
                        gbout.write(line)
            os.remove(gb_file)
            if delete == 0:
                os.remove(gb_file+'.out')
            temp_name2 = gb_file.replace('_out.fas', '.fasta')
            os.rename(gb_file+'.out', temp_name2)
####################################################
seq2 = {}
sid2 = ""
with open(aln.fasta) as INI:
    for line in INI:
        if line.startswith(">"):
            sid2 = line.strip()[1:]
        else:
            seq2[sid2] = seq2.get(sid2, "") + line.strip()

query = list(seq2.keys())
end_name = ""
with open(out, "w") as OUT:
    for id in list:
        new = id.split("-")
        if "Query" in new[0]:
            OUT.write(f">{new[1]}_{new[2]}_{new[3]}\n")
        else:
            OUT.write(f">{new[0]}_{new[1]}_{new[2]}\n")
        OUT.write(seq2.get(new[0], "") + "\n")
