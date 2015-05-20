"""merge counts and compare individual and the personalized genome results"""

import scipy as SP
import io
import sys
import os
import pdb
import re
import time
import glob
import string
import time
import h5py
import re
import cPickle
import pysam
import pandas

if __name__ == '__main__':

    reference = './../reference'


    gff_file = './../reference/b37_g1k_phase2/personal_transcriptome/gencode.v19.annotation_trunc.gtf'
    gff_cache = gff_file+'.pickle'

    exons = 'exons' in sys.argv

    if exons:
        out_GRCH37 = './../reference/b37_g1k_phase2/counts/exons'
        out_SNP_maternal = './../reference/1kgp3-svs-pass_NA12878_hg19_150109_snpsIndels/counts/exons_maternal'
        out_SNP_paternal = './../reference/1kgp3-svs-pass_NA12878_hg19_150109_snpsIndels/counts/exons_paternal'
        out_SV_maternal = './../reference/1kgp3-svs-pass_NA12878_hg19_150109_snpsIndelsSVs/counts/exons_maternal'
        out_SV_paternal = './../reference/1kgp3-svs-pass_NA12878_hg19_150109_snpsIndelsSVs/counts/exons_paternal'
        out_dir = 'out/exons'
    else:
        out_GRCH37 = './../reference/b37_g1k_phase2/counts/genes'
        out_SNP_maternal = './../reference/1kgp3-svs-pass_NA12878_hg19_150109_snpsIndels/counts/genes_maternal'
        out_SNP_paternal = './../reference/1kgp3-svs-pass_NA12878_hg19_150109_snpsIndels/counts/genes_paternal'
        out_SV_maternal = './../reference/1kgp3-svs-pass_NA12878_hg19_150109_snpsIndelsSVs/counts/genes_maternal'
        out_SV_paternal = './../reference/1kgp3-svs-pass_NA12878_hg19_150109_snpsIndelsSVs/counts/genes_paternal'
        out_dir = 'out/genes'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    [Rgene,Rexon] = cPickle.load(open(gff_cache,'rb'))
    if exons:
        gff = Rexon
    else:
        gff = Rgene
       
    element_ids = gff.keys()
    RV = []
    RV_file_exist = []
    RV_file = []
    #element_ids = element_ids[0:10]
    for element_id in element_ids: 
        #check that result files exist for all alignments
        chrom = gff[element_id][1]
        result_file = element_id+'.pickle'
        count_file_GRCH37   = os.path.join(out_GRCH37,'chr%s' % chrom,result_file)
        count_file_SNP_maternal = os.path.join(out_SNP_maternal,'chr%s_maternal' % chrom,result_file)
        count_file_SNP_paternal = os.path.join(out_SNP_paternal,'chr%s_paternal' % chrom,result_file)
        count_file_SV_maternal = os.path.join(out_SV_maternal,'chr%s_maternal' % chrom,result_file)
        count_file_SV_paternal = os.path.join(out_SV_paternal,'chr%s_paternal' % chrom,result_file)

        if ((not os.path.exists(count_file_GRCH37)) or (not os.path.exists(count_file_SNP_maternal)) or (not os.path.exists(count_file_SNP_paternal)) or (not os.path.exists(count_file_SV_maternal)) or (not os.path.exists(count_file_SV_paternal))):
            print "skip: %s" % element_id
            RV_file_exist.append([element_id,os.path.exists(count_file_GRCH37),os.path.exists(count_file_SNP_maternal),os.path.exists(count_file_SNP_paternal),os.path.exists(count_file_SV_maternal),os.path.exists(count_file_SV_paternal)])
            RV_file.append([element_id,count_file_GRCH37,count_file_SNP_maternal,count_file_SNP_paternal,count_file_SV_maternal,count_file_SV_paternal])
            continue
        #1. load lists
        count_GRCH37 = cPickle.load(open(count_file_GRCH37,'rb'))
        count_SNP_maternal = cPickle.load(open(count_file_SNP_maternal,'rb'))
        count_SNP_paternal = cPickle.load(open(count_file_SNP_paternal,'rb'))
        count_SV_maternal = cPickle.load(open(count_file_SV_maternal,'rb'))
        count_SV_paternal = cPickle.load(open(count_file_SV_paternal,'rb'))
        
        count_SNP = SP.union1d(count_SNP_maternal,count_SNP_paternal)
        count_SV = SP.union1d(count_SV_maternal,count_SV_paternal)
        count_intersect_GRCH37_SNP  = SP.intersect1d(count_SNP,count_GRCH37)
        count_intersect_GRCH37_SV  = SP.intersect1d(count_SV,count_GRCH37)
        count_intersect_SNP_SV  = SP.intersect1d(count_SNP,count_SV)

        count_ex_GRCH37_SNP = SP.setdiff1d(count_GRCH37,count_SNP)
        count_ex_GRCH37_SV = SP.setdiff1d(count_GRCH37,count_SV)
        count_ex_SNP_GRCH37 = SP.setdiff1d(count_SNP,count_GRCH37)
        count_ex_SV_GRCH37 = SP.setdiff1d(count_SV,count_GRCH37)
        count_ex_SNP_SV = SP.setdiff1d(count_SNP,count_SV)
        count_ex_SV_SNP = SP.setdiff1d(count_SV,count_SNP)
    
        #store a couple of things
        rv = []
        rv = {'element_id': element_id,'count_ref': len(count_GRCH37),'count_SNP_maternal':len(count_SNP_maternal),'count_SNP_paternal':len(count_SNP_paternal),'count_SV_maternal':len(count_SV_maternal),'count_SV_paternal':len(count_SV_paternal),'count_SNP':len(count_SNP),'count_SV':len(count_SV),'count_intersect_GRCH37_SNP':len(count_intersect_GRCH37_SNP),'count_intersect_GRCH37_SV':len(count_intersect_GRCH37_SV),'count_intersect_SNP_SV':len(count_intersect_SNP_SV),'count_ex_GRCH37_SNP':len(count_ex_GRCH37_SNP),'count_ex_GRCH37_SV':len(count_ex_GRCH37_SV),'count_ex_SNP_GRCH37':len(count_ex_SNP_GRCH37),'count_ex_SV_GRCH37':len(count_ex_SV_GRCH37),'count_ex_SNP_SV':len(count_ex_SNP_SV),'count_ex_SV_SNP':len(count_ex_SV_SNP)}
        RV.append(rv)
        pass
    #dump results
    RV = pandas.DataFrame(RV)
    RV.to_pickle(os.path.join(out_dir,'summary.pickl'))
