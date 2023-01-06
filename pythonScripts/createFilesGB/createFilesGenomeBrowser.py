## Run this file in a folder inside the 'allBigWig' folder that has a 'hg38' folder inside

import os

print('**********STARTED**********\n')
# ************************************************************************************
# ****************************** Paths and sample names ******************************
# ************************************************************************************

folderName = os.getcwd().split('/')[-1]

print('folderName: ', folderName)

myPath = os.getcwd() + '/hg38'

print('myPath: ', myPath, '\n')

sampleNames = [f.replace('_F.bw', '') for f in os.listdir(myPath) if '_F.bw' in f]
sampleNames.sort()
print('sampleNames', sampleNames, '\n')

email = 'analibanomonteiro@medicina.ulisboa.pt'

# ************************************************************************************
# **************************** Folder/files organization ****************************
# ************************************************************************************

"""
 myFolder\
   hg38\
      *.bw
      trackDb.txt
   genomes.txt
   hub.txt
"""

# ************************************************************************************
# ***************************** Files contents/templates *****************************
# ************************************************************************************

genomes_content =\
'genome hg38\n\
trackDb hg38/trackDb.txt'

hub_template =\
'hub {%folderName%}\n\
shortLabel {%folderName%}\n\
longLabel {%folderName%}\n\
genomesFile genomes.txt\n\
email {%email%}'

trackDb_template =\
'track {%sampleName%}\n\
container multiWig\n\
shortLabel {%sampleName%}\n\
longLabel {%sampleName%}\n\
type bigWig 0 60000\n\
viewLimits -10000:10000\n\
visibility full\n\
maxHeightPixels 160:120:11\n\
aggregate solidOverlay\n\
showSubtrackColorOnUi on\n\
windowingFunction maximum\n\
alwaysZero on\n\
priority 1.4\n\
configurable on\n\
autoScale on\n\
\n\
track {%sampleName%}_F\n\
bigDataUrl {%sampleName%}_F.bw\n\
shortLabel {%sampleName%}_F\n\
longLabel {%sampleName%}_F\n\
parent {%sampleName%}\n\
type bigWig\n\
color 0,0,255\n\
\n\
track {%sampleName%}_R\n\
bigDataUrl {%sampleName%}_R.bw\n\
shortLabel {%sampleName%}_R\n\
longLabel {%sampleName%}_R\n\
parent {%sampleName%}\n\
type bigWig\n\
color 255,0,0\n\
\n\n'

# ************************************************************************************
# *********************************** Write Files ************************************
# ************************************************************************************

with open('genomes.txt', 'w') as myFile:
    myFile.write(genomes_content)

print('File genomes.txt written in {}'.format(folderName), '\n')

with open('hub.txt', 'w') as myFile:
    myFile.write(hub_template.replace('{%folderName%}', folderName).replace('{%email%}', email))

print('File hub.txt written in {}'.format(folderName), '\n')

trackDb_content = ''
for sample in sampleNames:
    trackDb_content += trackDb_template.replace('{%sampleName%}', sample)

with open('hg38/trackDb.txt', 'w') as myFile:
    myFile.write(trackDb_content)

print('File trackDb.txt written in {}'.format(folderName + '/hg38'), '\n')


print('**********FINISHED!!**********')




