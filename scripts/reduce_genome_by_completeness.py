import os,sys,glob,random,argparse
from collections import Counter
from statistics import mean 

def breakGenome(gffkoTsvFN, targetCompleteness, originalCompleteness, wFN):
    if targetCompleteness > originalCompleteness: raise ValueError("Target completeness is higher than original completeness, but this script can't extrapolate gene contents")
    gffkoTsvFP = open(gffkoTsvFN,'r')
    columns = gffkoTsvFP.readline().split('\t')
    contigs, locusTag2ko = {}, {}
    for line in gffkoTsvFP:
        #print(line)
        contigID, geneStart, geneEnd, geneDir, category, locusTag, geneID, geneDec, koID = line.rstrip('\n').split('\t')
        #print([contigID, geneStart, geneEnd, geneDir, category, locusTag, geneID, geneDec, koID])
        contigs.setdefault(contigID,[])
        contigs[contigID].append(locusTag)
        locusTag2ko[locusTag] = koID 
    nGene = sum([len(contig) for contig in contigs.values()])
    targetNGene = int(nGene * targetCompleteness/originalCompleteness)
    print(f"{targetNGene} genes sampled from {nGene} genes")
    selected = []
    contigIDs = list(contigs.keys())
    random.shuffle(contigIDs)
    for contigID in contigIDs:
        contig = contigs[contigID]
        if len(selected) + len(contig) == targetNGene: break
        elif len(selected) + len(contig) > targetNGene:
            nMoreGene = targetNGene - len(selected)
            startIDX = random.randint(0,len(contig) - nMoreGene)
            selected.extend(contig[startIDX:startIDX+nMoreGene])
            break
        else: selected.extend(contig)

    koFP = open(wFN,'w')
    for locusTag in sorted(selected):
        if locusTag2ko[locusTag] != "":
            koFP.write(f"{locusTag}\t{locusTag2ko[locusTag]}\n")
    gffkoTsvFP.close(); koFP.close()

def main():
    parser = argparse.ArgumentParser(description="Reduce genome gene list based on target completeness using merged KO + annotation file.")
    parser.add_argument('--input_tsv', type=str, required=True, help="Input merged annotation file (.tsv) with locus_tag and KO columns.")
    parser.add_argument('--output', type=str, required=True, help="Path of output ko file")
    parser.add_argument('--target_completeness', type=float, required=True, help="Target genome completeness percentage (e.g., 70 for 70%)")
    parser.add_argument('--original_completeness', type=float, required=True,help="Genome completeness of query genome (e.g., 70 for 70%)")
    parser.add_argument('--seed', type=int, default=12345, help="Optional random seed for reproducibility")
    args = parser.parse_args()

    random.seed(args.seed)
    breakGenome(args.input_tsv, args.target_completeness, args.original_completeness, args.output)

if __name__ == "__main__": main()
