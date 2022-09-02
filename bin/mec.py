#!/usr/bin/env python3

import os,sys
import pysam
import re
import collections
import numpy as np

def get_read_size_distribution(all_reads,minmapq):

    """
    Function to calculate global insert size distribution across the whole assembly
    Return a frequency table of insert sizes as a  dictionary with key = insert size, value = frequency

    """

    frq = collections.defaultdict(int) # dictionary of insert sizes
    found = {}
    for read in all_reads:
        # accept read based on mapq, contig alignemnt and insert size
        if (read.mapq > minmapq) and (read.rnext == read.tid):
            if read.qname in found and found[read.qname][0]==read.tid:
                mate = found[read.qname]
                isize = abs(max( mate[1]+mate[2]-read.pos,read.pos+read.rlen-mate[1]))
                frq[isize] += 1
            else:
                found[read.qname] = (read.tid,read.pos,read.rlen)
    return frq

def MAD(frq):

    """

    Function to calculate median and median absolute deviation (MAD) from a dictionary of frequencies.

    Return a tuple of (mean, MAD).

    Arguments:
    frq: a dictionary of frequencies: key = insert size, value = frequency
     
    """

    all_lengths = []
    for k, v in frq.iteritems():
        new_vals = [k] * int(v)
        all_lengths.extend(new_vals)
    all_lengths = np.array(sorted(all_lengths))
    mid = len(all_lengths)/2
    median = all_lengths[mid]
    residuals = sorted(abs(all_lengths - median)) # difference between val and median
    MAD = residuals[mid] # median of residuals
    isize_sd = int(1.4826 * MAD)
    return median, isize_sd

def GCProportion(line):
    c = line.count("g") + line.count("c")+line.count("G") + line.count("C")
    t = len(line)*1.0
    if t==0: return 0
    return c/t

def get_all_interval(samfile,s,d,minq,alpha,beta,gamma,seq,GCavgRate):

    
    """
        Function to get all of the interval containing no misjoin errors among all of contigs.

        Arguments:
        filename: the bam file path or name
     """
    # getting reads from bamfile
    references = samfile.references  #type:list, item type:str
    lengths = samfile.lengths

    if s==0 and d==0:
        all_reads = samfile.fetch()
        frq = get_read_size_distribution(all_reads,minq)
        s,d = MAD(frq)
    interval = {}
    interval_valid = []

    for k,v in zip(references,lengths):
        reads = samfile.fetch(k, 0, v)
        interval_valid = get_valid_interval(k,v,reads,s,d,minq,alpha,beta,gamma,seq,GCavgRate)
        interval[k] = interval_valid

    return interval

def get_valid_interval(cid,length,reads,mu,sigma,minq,alpha,beta,gamma,seq,GCavgRate):

    left_discordant = []
    right_discordant = []
    frag_coverage = [0]*length
    
    for read in reads:
        if read.cigarstring!=None:
            if read.rnext == read.tid:
                if read.rlen<100:
                    if mu<1000:
                        if read.is_proper_pair and read.next_reference_start < read.reference_start:
                            if read.next_reference_start+read.rlen < read.reference_start:
                                for i in xrange(read.next_reference_start+read.rlen,read.reference_start):
                                    frag_coverage[i]+=1
                            else:
                                for i in xrange(read.next_reference_start,read.reference_end):
                                    frag_coverage[i]+=1
                                    
                    else:
                        if mu-3*sigma<= abs(read.isize) <= mu+3*sigma and read.next_reference_start < read.reference_start:
                            if read.next_reference_start+read.rlen < read.reference_start:
                                for i in xrange(read.next_reference_start+read.rlen,read.reference_start):
                                    frag_coverage[i]+=1
                            else:
                                for i in xrange(read.next_reference_start,read.reference_end):
                                    frag_coverage[i]+=1
                else:
                    if mu<1000:
                        if read.is_proper_pair and read.next_reference_start < read.reference_start:
                            for i in xrange(read.next_reference_start,read.reference_end):
                                frag_coverage[i]+=1
                    else:
                        if mu-3*sigma<= abs(read.isize) <= mu+3*sigma and read.next_reference_start < read.reference_start:
                            for i in xrange(read.next_reference_start,read.reference_end):
                                frag_coverage[i]+=1
            else:
                if read.is_reverse:
                    right_discordant.append((read.next_reference_name,read.pnext,read.reference_name,read.reference_start,read.reference_end,read.isize))
                else:
                    left_discordant.append((read.reference_name,read.reference_start,read.reference_end,read.next_reference_name,read.pnext,read.isize))
                                 
    frag_threshold = sum(frag_coverage)*1.0/length
    frag_pos = [i for i in xrange(length) if frag_coverage[i]>alpha*frag_threshold]
    frag_interval = pos_to_interval(frag_pos,100)

    interval = []
    if len(frag_interval)>=2:
        interval_len = len(frag_interval)
        p0 = 0
        for i in xrange(interval_len):
            if i==interval_len-1:
                p1 = length-1
                interval.append((p0,p1))
            else:
                error_pos = [frag_coverage[j] for j in xrange(frag_interval[i][1]+1,frag_interval[i+1][0])]
                break_pos01 = error_pos.index(min(error_pos))
                p1 = frag_interval[i][1]+break_pos01
                interval.append((p0,p1))
                p0 = p1+1
                p1 = frag_interval[i+1][1]
    else:
        interval = [(0,length-1)]
        
    interval_valid=[]
    before_discordant=[]
    after_discordant=[]
    discordant=[]
    discordant_list = []
    accordante_rate=[0]*(len(interval)-1)
    if len(interval)>=2:
        p0 = interval[0][0]
        p1 = interval[0][1]
        interval_len = len(interval)
        for i in xrange(interval_len-1):

            before_ctg =  []
            after_ctg = []
            
            start = max(interval[i][1]-mu,interval[i][0])
            end = min(interval[i+1][0]+mu,interval[i+1][1])
            
            after_discordant.append([item for item in left_discordant if (item[1]>=start and item[2]<=interval[i+1][0]) ])
            before_discordant.append([item for item in right_discordant if (item[3]>=interval[i+1][0] and item[4]<=end) ]) 

            for item in left_discordant:
                if (item[1]>=start and item[2]<=interval[i+1][0]):
                    after_ctg.append(item[3])

            for item in right_discordant:
                if (item[3]>=interval[i+1][0] and item[4]<=end):
                    before_ctg.append(item[0])

            couter01 = collections.Counter(after_ctg)
            couter02 = collections.Counter(before_ctg)
            
            n1=n2=0
            if len(couter01.most_common(1))!=0:
                n1 = couter01.most_common(1)[0][1]
            if len(couter02.most_common(1))!=0:
                n2 = couter02.most_common(1)[0][1]

            isize_discordant = 0
            if n1>=n2:
                isize_discordant += n1
                discordant.append(after_discordant[i])
            else:
                isize_discordant += n2
                discordant.append(before_discordant[i])
            
            isize_accordant = frag_coverage[interval[i+1][0]]
            accordante_rate[i] = isize_accordant*1.0/max((isize_accordant+isize_discordant),0.01)
                
            if accordante_rate[i] >= beta or isize_discordant<=3: #or n1<=3 or n2<=3 or discordante_rate[i] <= 0.65
                pos = i+1
                p1 = interval[pos][1]
                if pos == interval_len-1:
                    interval_valid.append((p0,p1))
            else:
                discordant_list.append((accordante_rate[i],interval[i+1][0],isize_accordant,isize_discordant))
                discordant_list.append((len(discordant[i]),discordant[i]))
                interval_valid.append((p0,p1))
                pos= i+1
                p0 = interval[pos][0]
                p1 = interval[pos][1]
                if pos == interval_len-1:
                    interval_valid.append((p0,p1))       
    else:
        interval_valid = interval
    
    GC_rate=[]
    intervals=[]
    if len(interval_valid)>=2:
        interval_len = len(interval_valid)
        p0 = interval_valid[0][0]
        p1 = interval_valid[0][1]
        for i in xrange(interval_len-1):
            start = max(interval_valid[i][1]-mu/2,0)
            end = min(interval_valid[i+1][0]+mu/2,length-1)

            GC_rate.append(GCProportion(seq[cid][start:interval_valid[i][1]]))
            GC_rate.append(GCProportion(seq[cid][interval_valid[i+1][0]:end]))
            GC_rate.append(GCProportion(seq[cid][start:end]))
            
            p_gc = GCProportion(seq[cid][start:end])
            l_gc = GCProportion(seq[cid][start:interval_valid[i][1]])
            r_gc = GCProportion(seq[cid][interval_valid[i+1][0]:end])
            
            min_gc = GCavgRate-gamma
            max_gc = GCavgRate+gamma
            if (min(p_gc,l_gc,r_gc)<=min_gc or max(p_gc,l_gc,r_gc)>=max_gc):
                pos = i+1
                p1 = interval_valid[pos][1]
                if pos == interval_len-1:
                    intervals.append((p0,p1))
            else:
                intervals.append((p0,p1))
                pos= i+1
                p0 = interval_valid[pos][0]
                p1 = interval_valid[pos][1]
                if pos == interval_len-1:
                    intervals.append((p0,p1))
    else:
        intervals=interval_valid
                   
    return intervals
                
def pos_to_interval(pos_list,minlen):
    
    """
        Function to invert pos list to the interval if pos is contiguous.
        return: the interval that does not contain misassembly error.

        Arguments:
        pos_list: the list of pos

     """
    interval = []
    if len(pos_list)>=2:
        start = pos_list[0]
        end = pos_list[0]
        for i in pos_list[1:]:
            if i-end == 1:
                end = i
            else:
                if end - start >minlen:
                    interval.append((start, end))
                start = i
                end = i
        if end - start >minlen:
            interval.append((start, end))
    else:
        interval = []
    return interval



def NProportion(line):
    c = line.count("n") + line.count("N")
    t = len(line)*1.0
    if t==0: return 0
    return c/t

def GCProportion(line):
    c = line.count("g") + line.count("c")+line.count("G") + line.count("C")
    t = len(line)*1.0
    if t==0: return 0
    return c/t

def parse_contig(contigPath):

    """
    Function to parse the input contig file.
    Return a dict of (contigid:sequence).

    Arguments:
    contigPath: the filename or the path  of contigfile.
    flag: the choice to return seq or contigid.

    """
    contigfile = open(contigPath)
    seq={}
    contigid=""
    for line in contigfile:
        if not line: break
        if re.findall(r'^>',line):
            contigid = line.split(" ")[0][1:-1]
            seq[contigid] = ""
        else:
            seq[contigid] += line.strip()
    contigfile.close()
    return seq

def output_contig(seq,dic_interval,contigid,outfile):
    
    """
    Function to split the contigs at positions identified as assembly errors and write a new fasta file containing all contigs.
    
    Arguments:
    outfile: name of the new fasta file (including filepath).
    interval: list of range that have no misassemblies errors.
    contigid: list of id of the all contig.
    seq: the dictionary of the all contigs, key=id of the contig, value=the seq of each contig.
    """  
    out = open(outfile, "w")
    split_num = 0
    
    for cid in contigid:
        if len(dic_interval[cid])<1:
            print (cid)
        assert len(dic_interval[cid])>=1
        if len(dic_interval[cid])==1:
            p0=dic_interval[cid][0][0]
            p1=dic_interval[cid][0][1]
            if NProportion(seq[cid][p0:p1+1])>0.7: continue
            out.write(">"+cid+"\n"+seq[cid][p0:p1+1]+"\n")
        else:
            intervals = dic_interval[cid]
            cnt = 0
            for item in intervals:  
                if NProportion(seq[cid][item[0]:item[1]+1])>0.7: continue
                out.write(">"+cid+"_"+str(cnt)+"\n"+seq[cid][item[0]:item[1]+1]+"\n")
                cnt = cnt+1
            split_num += len(intervals)-1
    out.close()
    return split_num
                
def main():
    # read command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Tool to identify and correct misassemblies in de novo assemblies')
    parser.add_argument('-i', metavar='infasta', type=str, help='input assembly fasta file')
    parser.add_argument('-o', metavar='outfasta', type=str, help='output the correct file name')
    parser.add_argument('-bam', metavar='bam', type=str, default=None,help='index bam file for alignment')
    parser.add_argument('-m', metavar='mu', type=int,default=0,help='one read pair for pair end reads')
    parser.add_argument('-s', metavar='sigma', type=int,default=0,help='another read pair for pair end reads')
    parser.add_argument('-q', metavar ='minmapq', type=int,default=40, help='Minimum mapping quality value. Default value: 40')
    parser.add_argument('-a', metavar ='alpha', type=float, default=0.4, help='The percentage of the average of the fragment coverage. Default value: 0.4')
    parser.add_argument('-b', metavar ='beta', type=float, default=0.5, help='One cutoff for removing false misassemblies. Default value: 0.5')
    parser.add_argument('-g', metavar ='gamma', type=float, default=0.2, help='One parameter for determining high or low GC content. Default value: 0.2.')
    args = parser.parse_args()

    print ("start parse contig file")
    seq = parse_contig(args.i)

    bamfile = args.bam
    samfile = pysam.Samfile(bamfile, "rb" )
    contigid = samfile.references
    GCRate = 0.0
    for cid in contigid:
        GCRate += GCProportion(seq[cid])
    GCavgRate = GCRate/len(seq)
    interval = get_all_interval(samfile,args.m,args.s,args.q,args.a,args.b,args.g,seq,GCavgRate)

    print ("output all the intervals for each contig")
    interval_file01 = "./intervals.txt"
    f01 = open(interval_file01,"w")
    for cid in contigid:
        f01.write(">"+cid+"\n")
        for item in interval[cid]:
            f01.write(str(item[0])+" "+str(item[1])+"\n")
    f01.close()
 
    print ("output contig file")
    split_num = output_contig(seq,interval,contigid,args.o)
    
    print ("split_num",split_num)
    print ("\n")
    
if __name__ == "__main__":
    main()