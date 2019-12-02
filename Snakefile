import os, argparse, subprocess, re, sys, time

#from python.scan import scan_me

#subprocess.call("module load aocc", shell = True)
shell.prefix("module load aocc; ")

#set directory path for files to be analyzed
direc = "/home/swisotsky/BUSTED-SRV/32_subset/"

allFiles = list()

numbers = re.compile ("^\.[0-9]+$")

#generate list of files that don't already have BUSTED-SRV.json results file

for root, dirs, files in os.walk(direc):
  for each_file in files:
    name, ext = os.path.splitext(each_file)
    if len(ext) > 0 and ext in ['.mt', '.nex'] or numbers.match (ext)  :
      existing = os.path.join (direc, name + ext + ".BUSTED-SRV.json")
      if not os.path.isfile (existing):
        file = os.path.join(root, name + ext)
        allFiles.append(file)


#print(allFiles)
#files, = glob_wildcards("test/{file}")

#need this so it runs for each file in the list
rule all:
        input:
                expand("{file}.BUSTED-SRV.json", file=allFiles)



### this rule runs the HYPHY analysis with the specified parameters and input and output files ###
rule BS32: 
    input: 
        "{file}"
    output: 
        "{file}.BUSTED-SRV.json"
    shell: '/home/swisotsky/bin/hyphy-dev/hyphy/HYPHYMP LIBPATH=/home/swisotsky/bin/hyphy-dev/hyphy/res/  busted --alignment {input} --srv Yes --grid-size 2000 --starting-points 10 --output {output}
