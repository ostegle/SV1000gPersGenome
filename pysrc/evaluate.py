import pandas 
import scipy as sp 
import os 
import sys 
import pdb
import cPickle


if __name__ == '__main__':
    exons = ('exons' in sys.argv)

    gff_file = './../reference/b37_g1k_phase2/personal_transcriptome/gencode.v19.annotation_trunc.gtf'
    gff_cache = gff_file+'.pickle'
    variant_info = cPickle.load(open('/homes/stegle/research/projects/1000GenomesRNASeq/data/personalized_genome/reference/variants/NA12878.wgs.mergedSV.v5.20130502.svs.genotypes.redun.auto.SVdefined.sorted.pass.vcf.pickle','rb'))
    variant_info_start = sp.array(variant_info[:,2],dtype='int')
    variant_info_end = variant_info_start+sp.array(variant_info[:,3],dtype='int')
        
    #filter features to things that have SVs?
    filter_features = 'filter' in sys.argv
    ws = 100000
    ws = 1000
    #ws = 000

    [Rgene,Rexon] = cPickle.load(open(gff_cache,'rb'))
    if exons:
        gff = Rexon
    else:
        gff = Rgene
   
    if exons:
        S = pandas.load('./out/exons/summary.pickl')
        if filter_features:
            out_file = './out/exons/result_filtered_%d.csv' % ws
        else:             
            out_file = './out/exons/result.csv'
    else:
        S = pandas.load('./out/genes/summary.pickl')
        if filter_features:
            out_file = './out/genes/result_filtered_%d.csv' % ws
        else:             
            out_file = './out/genes/result.csv'

    I = sp.ones([S.shape[0]],dtype='bool')
    if filter_features:
        for i in xrange(S.shape[0]):
            element_id = S['element_id'][i]
            element_anno = gff[element_id]
            Ich = variant_info[:,1]==element_anno[1]
            Ioverlap_start = (element_anno[2]>=variant_info_start[Ich]-ws) & (element_anno[2]<variant_info_end[Ich]+ws)
            Ioverlap_end = (element_anno[3]>=variant_info_start[Ich]-ws) & (element_anno[3]<variant_info_end[Ich]+ws)
            I[i] = Ioverlap_start.any() | Ioverlap_end.any()
            pass
    pass

    #subset
    S=S.iloc[I]


    M = []
    columns = ['GRCh37','NA12878 (SNPS only)', 'NA12878 (SNPS only) gain','NA12878 (SNPS only) loss','NA12878 (SVs)','NA12878 (SVs) gain','NA12878(SVs) loss','NA12878 SVs/SNPs gain','NA12878 SVs/SNPs loss','NA12878 (SNPS only) delta','NA12878 (SVs) delta','NA12878 SVs/SNPs delta']
    rows = []

    #elements to consider:
    m0 = [S['count_ref'],S['count_SNP'],S['count_ex_SNP_GRCH37'],S['count_ex_GRCH37_SNP'],S['count_SV'],S['count_ex_SV_GRCH37'],S['count_ex_GRCH37_SV'],S['count_ex_SV_SNP'],S['count_ex_SNP_SV']]

    #1. total read count
    rows.append('total read count')
    M.append([f for f in m0])
    #2. features with >10 
    rows.append('>10 reads')
    M.append([(f>10) for f in m0])
    #2. features with >100 
    rows.append('>100 reads')
    M.append([(f>100) for f in m0])
    #2. features with >1000 
    rows.append('>1000 reads')
    M.append([(f>1000) for f in m0])
    #3. featuers > 10 reads and 1 fold change
    rows.append('>10 reads, 1 fold change')
    m = [sp.array([False])]*9
    Ichange = (S['count_ex_SNP_GRCH37']/(1.0*S['count_ref'])>1)*(S['count_ex_SNP_GRCH37']>10)
    m[2] = Ichange
    Ichange = (S['count_ex_GRCH37_SNP']/(1.0*S['count_SNP'])>1)*(S['count_ex_GRCH37_SNP']>10)
    m[3] = Ichange
 
    Ichange = (S['count_ex_SV_GRCH37']/(1.0*S['count_ref'])>1)*(S['count_ex_SV_GRCH37']>10)
    m[5] = Ichange
    Ichange = (S['count_ex_GRCH37_SV']/(1.0*S['count_SV'])>1)*(S['count_ex_GRCH37_SV']>10)
    m[6] = Ichange
    Ichange = (S['count_ex_SV_SNP']/(1.0*S['count_SNP'])>1)*(S['count_ex_SV_SNP']>10)
    m[7] = Ichange
    Ichange = (S['count_ex_SNP_SV']/(1.0*S['count_SV'])>1)*(S['count_ex_SNP_SV']>10)
    m[8] = Ichange
    M.append(m)
    
    for row in M:
        row.append(sp.array(row[2]) | sp.array(row[3]))
        row.append(sp.array(row[5]) | sp.array(row[6]))
        row.append(sp.array(row[7]) | sp.array(row[8]))
    #combine gain/loss columns into delta columns       
    pdb.set_trace()
    for i in xrange(len(M)):
        for j in xrange(len(M[0])):
            M[i][j] = M[i][j].sum()


    
    M=pandas.DataFrame(M,columns=columns,index=rows)
    #reoder columns
    columns=sp.array(M.columns.tolist())
    columns=columns[[0,1,9,2,3,4,10,5,6,11,7,8]]
    M = M[columns]

    M.to_csv(out_file)
