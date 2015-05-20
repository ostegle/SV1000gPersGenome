"""merge single cell RNA-Seq fastq files"""

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

if __name__ == '__main__':
    if 0:
        in_file = '/homes/stegle/research/users/stegle/hipsci/single-cell/data/pilot3/iPS/alignments/sc_16.bam'
        gff_file = '/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/reference/GRCH37/gencode.v19.annotation.gtf'
        out_dir = '/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/rna-seq/counts/GRCH37'
    base_chr = ''
    in_file = sys.argv[1]
    gff_file = sys.argv[2]
    out_dir = sys.argv[3]

    gff_cache = gff_file+'.pickle'
    if (not os.path.exists(gff_cache)) or 'recalc' in sys.argv:
        Rgene ={}
        Rexon ={}
        re_exon_number = re.compile('.*exon_number (\d*).*')
        re_gene_id     = re.compile('.*gene_id "(.*?)".*')
        #M = SP.loadtxt(gff_file,dtype='str',delimiter='\t')
        with open(gff_file,'r') as f:
            for line in f:
                m = line.split('\t')
                if m[2]=='gene':
                    atype='gene'
                elif m[2]=='exon':
                    atype='exon'
                else:
                    continue
                chrom = m[0]
                start = int(m[3])
                end   = int(m[4])
                meta = m[8]
                if len(meta)==0:
                    continue
                gene_id = re_gene_id.match(meta).group(1)
                if atype=='gene':
                    Rgene[gene_id] = [gene_id,chrom,start,end]
                elif atype=='exon':
                    exon_number = int(re_exon_number.match(meta).group(1))
                    gene_id+='_%d' % exon_number
                    Rexon[gene_id] = [gene_id,chrom,start,end]
        cPickle.dump([Rgene,Rexon],open(gff_cache,'wb'),-1) 
    else:
        [Rgene,Rexon] = cPickle.load(open(gff_cache,'rb'))

    #process exons or genes?
    if 'exon' in sys.argv:
        gff = Rexon
    else:
        gff = Rgene


    #open file
    f = pysam.Samfile( in_file, "rb" )

    gene_ids = gff.keys()

    for gene_id in gene_ids: 
        gene_chrom = gff[gene_id][1]
        gene_start = gff[gene_id][2]
        gene_end   = gff[gene_id][3]

        #fetch=f.fetch('%s' % gene_chrom,gene_start,gene_end)
        fetch=f.fetch('%s%s' % (base_chr,gene_chrom),gene_start,gene_end)
        read_ids = []
        for read in fetch:
            #only consider primary alignments:
            if (read.is_qcfail or read.is_duplicate or (not read.is_proper_pair) or read.is_secondary):
                continue
            qname = read.qname
            if read.is_read1:
                qname+='_1'
            else:
                qname+='_2'
            read_ids.append(qname)
            pass
        read_ids = SP.array(read_ids)
        #print "%s:%d" % (gene_id,len(read_ids))
        _out_dir = os.path.join(out_dir,'chr%s' % gene_chrom)
        if not os.path.exists(_out_dir):
            os.makedirs(_out_dir)
        out_file = os.path.join(_out_dir,gene_id+'.pickle')
        cPickle.dump(read_ids,open(out_file,'wb'))
