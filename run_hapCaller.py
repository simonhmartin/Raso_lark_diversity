#!/usr/bin/python


import sys, time, os, argparse, subprocess


parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bam", help="bam Name", action = "store", required = True)
parser.add_argument("-r", "--reference", help="reference fasta file", action='store', required = True)
parser.add_argument("--runName", help="Run name", action="store", default = "HC")

parser.add_argument("--threads", help="Threads for -nct option", type=int, action='store', default = 1)

parser.add_argument("--outDir", help="Output directory", action='store', default = "./")
parser.add_argument("--gatkDir", help="GATK .jar directory", action='store')
parser.add_argument("--tmpDir", help="Temporary directory for Java memory overflow", action='store', default = "./")
parser.add_argument("--javaRam", help="Amount of ram for java to use", action='store', default="2g")


###HaplotypeCaller specific arguments
parser.add_argument("-ploidy", help="Sample ploidy for HaplotypeCaller", type=int, action='store', default = 2)
parser.add_argument("--minReadsPerAlignmentStart", help="min Reads Per Alignment Start", type=int, action='store', default = 10)
parser.add_argument("--maxReadsInRegionPerSample", help="max Reads In Region Per Sample", type=int, action='store', default = 10000)

parser.add_argument("--dontUseSoftClippedBases", help="Hapcaller wil ignore soft clipped bases", action='store_true')
parser.add_argument("--heterozygosity", help="Heterozygosity prior for Hapcaller", type = float, action='store', default = 0.02)
parser.add_argument("--outputMode", help="Output mode for Hapcaller", action='store', choices = ["EMIT_ALL_SITES","EMIT_ALL_CONFIDENT_SITES","EMIT_VARIANTS_ONLY"], default = "EMIT_ALL_SITES")
parser.add_argument("-stand_call_conf", help="Confidence threshold for calling in HaplotypeCaller", type=int, action='store', default = 30)
parser.add_argument("--regions", help="Intervals to analyse", action='store')

###flow control arguments
parser.add_argument("--test", help="Just print commands", action="store_true", required = False)



args = parser.parse_args()

bam = args.bam
ref = args.reference


threads = args.threads


outDir = args.outDir + "/"
gatkDir = args.gatkDir + "/"
tmpDir = args.tmpDir + "/"


javaRam = args.javaRam


runName = args.runName


###HaplotypeCaller specific arguments
ploidy = args.ploidy
het = args.heterozygosity
outputMode = args.outputMode
callConf = args.stand_call_conf
if args.dontUseSoftClippedBases:
    softClipped = "--dontUseSoftClippedBases"
else:
    softClipped = ""

if args.regions:
    regions = "-L " + args.regions
else:
    regions = ""


###flow control
test = args.test

fileName = bam.rstrip(".bam")

#############################################################################################################################


hapCallerCommand = " ".join(["java", "-Xmx" + javaRam, "-jar", gatkDir + "GenomeAnalysisTK.jar",
                            "-T", "HaplotypeCaller", "-nct", str(threads),
                            "-R", ref, "-I", bam, "-o", outDir + ".".join([fileName,runName,"g","vcf.gz"]),
                            "-ploidy", str(ploidy),  softClipped, "--heterozygosity", str(het),
                            "--emitRefConfidence GVCF", "--max_alternate_alleles 6", "-variant_index_type LINEAR",
                            "-variant_index_parameter 128000", "--output_mode", outputMode,
                            "-stand_call_conf", str(callConf), 
                            "--minReadsPerAlignmentStart", str(args.minReadsPerAlignmentStart),
                            "--maxReadsInRegionPerSample", str(args.maxReadsInRegionPerSample),
                            regions, "2>", outDir + ".".join([fileName,runName,"log"]) ])

print >> sys.stderr, "\n","HaplotypeCaller command:","\n",hapCallerCommand,"\n"

if not test:
    
    print >> sys.stderr, "Starting HaplotypeCaller at", time.strftime("%d/%m/%Y"), time.strftime("%H:%M:%S")
    
    os.system(hapCallerCommand)
    
    #check for success
    if subprocess.check_output(["tail", "-1", outDir + ".".join([fileName,runName,"log"])]).split(".")[0] != "Runtime":
        print >> sys.stderr, "\nHaplotypeCaller failed - check log file.\n"
        sys.exit()
    
    print >> sys.stderr, "HaplotypeCaller complete at", time.strftime("%d/%m/%Y"), time.strftime("%H:%M:%S")


