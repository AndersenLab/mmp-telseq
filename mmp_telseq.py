#!/usr/bin/python

#SBATCH --job-name=mmp_tel
#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out


#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem=16384

#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/mmp_telseq/sra

import os, sys
import glob
import re
import subprocess
from subprocess import PIPE, Popen
from datetime import datetime

def file_exists(filename):
    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        return True
    else:
        return False

class EAV:
    """
    Very simple Entity-Attribute-Value Object
    """

    def __init__(self):
        self.entity = ""
        self.sub_entity = ""
        self.attribute = ""
        self.sub_attribute = ""
        self.value = ""
        self.timestamp = datetime.now()
        self.comment = ""
        self.file = None

    def __repr__(self):
        return "\nEntity:{self.entity}\n\
                Entity:{self.sub_entity}\n\
                Attribute:{self.attribute}\n\
                Sub-Attribute:{self.sub_attribute}\n\
                Value:{self.value}\n\
                timestamp:{self.timestamp}\n".format(**locals())

    def save(self):
        if self.file is None:
            raise Exception("No Log File Set")
        if not file_exists(self.file):
            write_header = True
        else:
            write_header = False
        with(open(self.file, "a")) as f:
            if write_header is True:
                f.write("entity\tsub_entity\tattribute\tsub_attribute\tvalue\tcomment\ttimestamp\n")
            line = '\t'.join(map(str,[self.entity,
                              self.sub_entity,
                              self.attribute,
                              self.sub_attribute,
                              self.value,
                              self.comment,
                              self.timestamp]))
            f.write(line + "\n")

def get_contigs(bam):
    header, err = subprocess.Popen(["samtools","view","-H",bam], stdout=PIPE, stderr=PIPE).communicate()
    if err != "":
        raise Exception(err)
    # Extract contigs from header and convert contigs to integers
    contigs = {}
    for x in re.findall("@SQ\WSN:(?P<chrom>[A-Za-z0-9_]*)\WLN:(?P<length>[0-9]+)", header):
        contigs[x[0]] = int(x[1])
    return contigs

def coverage(bam, mtchr = None):
    # Check to see if file exists
    if os.path.isfile(bam) == False:
        raise Exception("Bam file does not exist")
    contigs = get_contigs(bam)

    # Guess mitochondrial chromosome
    mtchr = [x for x in contigs if x.lower().find("m") == 0]
    if len(mtchr) != 1:
        mtchr = None
    else:
        mtchr = mtchr[0]

    coverage_dict = {}
    for c in contigs.keys():
        command = "samtools depth -r %s %s | awk '{sum+=$3;cnt++}END{print cnt \"\t\" sum}'" % (c, bam)
        coverage_dict[c] = {}
        coverage_dict[c]["Bases Mapped"], coverage_dict[c]["Sum of Depths"] = map(int,subprocess.Popen(command, stdout=PIPE, shell = True).communicate()[0].strip().split("\t"))
        coverage_dict[c]["Breadth of Coverage"] = coverage_dict[c]["Bases Mapped"] / float(contigs[c])
        coverage_dict[c]["Depth of Coverage"] = coverage_dict[c]["Sum of Depths"] / float(contigs[c])
        coverage_dict[c]["Length"] = int(contigs[c])

    # Calculate Genome Wide Breadth of Coverage and Depth of Coverage
    genome_length = float(sum(contigs.values()))
    coverage_dict["genome"] = {}
    coverage_dict["genome"]["Length"] = int(genome_length)
    coverage_dict["genome"]["Bases Mapped"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k != "genome"])
    coverage_dict["genome"]["Sum of Depths"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k != "genome"])
    coverage_dict["genome"]["Breadth of Coverage"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)
    coverage_dict["genome"]["Depth of Coverage"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k != "genome"]) / float(genome_length)

    if mtchr != None:
        # Calculate nuclear breadth of coverage and depth of coverage
        ignore_contigs = [mtchr, "genome", "nuclear"]
        coverage_dict["nuclear"] = {}
        coverage_dict["nuclear"]["Length"] = sum([x["Length"] for k,x in coverage_dict.iteritems() if k not in ignore_contigs ])
        coverage_dict["nuclear"]["Bases Mapped"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs])
        coverage_dict["nuclear"]["Sum of Depths"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs])
        coverage_dict["nuclear"]["Breadth of Coverage"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs]) / float(coverage_dict["nuclear"]["Length"])
        coverage_dict["nuclear"]["Depth of Coverage"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs]) / float(coverage_dict["nuclear"]["Length"])

        # Calculate the ratio of mtDNA depth to nuclear depth
        coverage_dict["genome"]["mt_ratio"] = coverage_dict[mtchr]["Depth of Coverage"] / float(coverage_dict["nuclear"]["Depth of Coverage"])

    # Flatten Dictionary
    coverage = []
    for k,v in coverage_dict.items():
        for x in v.items():
            coverage += [(k,x[0], x[1])]
    return coverage

line_num = int(sys.argv[1]) - 1


f=open('../strain_info.txt')
lines = f.readlines()
lines = [x.strip().split("\t") for x in lines]

line = lines[line_num]
strain_name = line[0].split(" ")[2]
strain_bp = line[0].split(" ")[4]
reference = "/lscr2/andersenlab/dec211/pyPipeline/genomes/WS245/c_elegans.PRJNA13758.WS245.genomic.fa.gz"

# Download sra files
"""for line in lines:
    for i in line[1:]:
        strain = line[0].split(" ")[2]
        length = line[0].split(" ")[4]
        if len(glob.glob("../telseq/{strain}.{length}*".format(**locals()))) == 0:
            print strain, length
            i06 = i[0:6]
            i09 = i[0:9]
            loc_string = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{i06}/{i09}/{i}.sra"
            loc_string = loc_string.format(**locals())
            print "downloading " +  loc_string.format(**locals())
            print "curl {loc_string}".format(**locals())
            os.system("curl {loc_string} > {i}.sra".format(**locals()))
"""

# Process SRA Files
for i in line[1:]:
    os.system("fastq-dump --split-files --gzip {i} ".format(i=i))
    os.system("rm /exports/people/andersenlab/dec211/ncbi/public/{i}.sra".format(i=i))
    #os.system("rm {i}".format(**locals()))
    # Generate read group
    RG = r'@RG\tID:{i}\tSM:{strain_name}'.format(**locals())
    # Align
    os.system(r"bwa mem -R '{RG}' -t 4 {reference} {i}_1.fastq.gz {i}_2.fastq.gz > ../bam/{i}.tmp.bam".format(i=i.replace(".sra",""), RG=RG, reference=reference))
    # Sort
    os.system("samtools sort -O bam -T ../bam/{i}.TEMP.bam -@ 4 ../bam/{i}.tmp.bam > ../bam/{i}.sorted.bam && samtools index ../bam/{i}.sorted.bam".format(**locals()))
    # Remove temporary BAM and fastq
    os.system("rm {i}_1.fastq.gz && rm {i}_2.fastq.gz".format(i=i))
    os.system("rm ../bam/{i}.tmp.bam".format(i=i))


# Combine processed BAM Files.
SRA_files = ' '.join(["../bam/" + x.replace(".sra","") + ".sorted.bam" for x in line[1:]])
if len(["../bam/" + x.replace(".sra","") + ".sorted.bam" for x in line[1:]]) > 1:
    merge_command = "samtools merge -f -@ 4 ../bam/{strain_name}.{strain_bp}.bam {SRA_files} && samtools index ../bam/{strain_name}.{strain_bp}.bam".format(**locals())
    os.system(merge_command)
else:
    os.system("mv {SRA_files} ../bam/{strain_name}.{strain_bp}.bam".format(**locals()))
    os.system("samtools index ../bam/{strain_name}.{strain_bp}.bam".format(**locals()))

for i in line[1:]:
    os.system("rm ../bam/{i}.sorted.bam && rm ../bam/{i}.sorted.bam.bai".format(i=i))

# Produce Coverage Statistics Here
bam = "../bam/{strain_name}.{strain_bp}.bam".format(**locals())

eav = EAV()
eav.file = "../eav.txt"
eav.entity = strain_name
eav.sub_entity = strain_bp


for contig, k,v in coverage(bam, "MtDNA"):
    eav.sub_attribute = contig + " (" + k + ")"
    eav.value = v
    eav.save()

# Run Telseq Here
#telseq -z 'AATCCG' -u $file.bam -o $file.telseq_elegans.AATCCG.noreadgroup.txt
os.system("telseq -m -z 'TTAGGC' -u ../bam/{strain_name}.{strain_bp}.bam -o ../telseq/{strain_name}.{strain_bp}.telseq_elegans.TTAGGC.noreadgroup.txt".format(**locals()))
os.system("telseq -m -z 'GTATGC' -u ../bam/{strain_name}.{strain_bp}.bam -o ../telseq/{strain_name}.{strain_bp}.telseq_elegans.GTATGC.noreadgroup.txt".format(**locals()))
#telseq -z 'GTCTAG' -u $file.bam -o $file.telseq_elegans.GTCTAG.noreadgroup.txt

# Delete sra file 
for i in line[1:]:
    #os.system("rm {i}.sra ".format(i=i))
    pass

# Delete bam file
os.system("rm ../bam/{strain_name}.{strain_bp}.bam".format(**locals()))