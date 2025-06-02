import sys,os,glob,random,argparse
from collections import Counter

def detachExtensions(fileName, depth=0, baseName=True):
    if baseName: fileName = os.path.basename(fileName)
    if depth == 0:
        return fileName.split('.')[0]
    elif depth < 0:
        return '.'.join(fileName.split('.')[:depth]) ## remain 2 components
    else:
        return '.'.join(fileName.split('.')[:-depth]) ## remove 2 components

def change2num(string):
    if type(string) != str: return string
    if string.strip('0123456789-.') != '' or string.count('.') > 1 \
       or string.count('-') > 1 or string == '':
        return string
    elif "." == string or "-" == string:
        return string
    elif '.' in string:
        return float(string)
    else:
        return int(string)

def parseGFF(gffFN):
    gffFP = open(gffFN,'r')
    geneList = []
    for line in gffFP:
#        print(line)
        if line.startswith(">"): break
        if line.startswith("#"): continue
        contigID, annot, category, gStart, gEnd, etc1, gDir, etc2, geneDec  = list(map(change2num,line.strip().split('\t')))
        if annot.startswith("prokka"): continue
        locusTag = geneDec.split(';')[0].split('=')[1]
        if "gene=" in geneDec:
            geneID = geneDec.split("gene=")[1]
            geneID = geneID[:geneID.find(';')]
        else: geneID = ''
        if 'product=' in geneDec: geneDec = geneDec.split("product=")[1]
        else: geneDec = ""
        geneList.append((contigID, gStart, gEnd, gDir, category, locusTag, geneID, geneDec))
    return geneList   

def main():
    parser = argparse.ArgumentParser(description="Join Prokka GFF and KO annotation to generate a merged TSV file.")
    parser.add_argument('--gff', type=str, required=True, help="Input .gff file from Prokka annotation")
    parser.add_argument('--ko', type=str, required=True, help="Input .ko file with locus_tag and KO ID per line")
    parser.add_argument('--output_tsv', type=str, required=True, help="Path to output merged .tsv file")

    args = parser.parse_args()
    gffFN, koFN, outFN = args.gff, args.ko, args.output_tsv  
    wFP = open(outFN,'w')
    wFP.write("Contig\tgeneStart\tgeneEnd\tgeneDir\tcategory\tLocusTag\tgeneID\tgeneDec\tkoID\n")
    geneList = parseGFF(gffFN)
    koD = {line.strip().split('\t')[0]:line.strip().split('\t')[1].split(',') for line in open(koFN) if len(line.strip().split('\t')) > 1}
    print(os.path.basename(gffFN),"#Gene",len(geneList),"#KOassingedGene",len(koD))
    for contigID, gStart, gEnd, gDir, category, locusTag,geneID, geneDec in geneList:
        koIDstr = ','.join(koD.get(locusTag,[]))
        wFP.write(f"{contigID}\t{gStart}\t{gEnd}\t{gDir}\t{category}\t{locusTag}\t{geneID}\t{geneDec}\t{koIDstr}\n")
    wFP.close()


if __name__ == "__main__": main()
