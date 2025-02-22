#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 18:05:49 2023

@author: Sachithra Yaddehige
"""

import argparse
from datetime import datetime
import re
import csv
import math
from collections import Counter
startTime = datetime.now()

class GBRAP():
    def __init__(self, input_path, output_path):
        self.input_path = input_path
        self.output_path = output_path
        
    def splitter (self,location):
        ''' splits the start & end cordinates'''
        split=location.split(':')
        if len(split)==1: #for cases where only one base position given in gbff
            start= int(split[0])
            value=(start, start)
        else:
            start,end=map(int,split)
            value=(start, end)
        return value
    
    def split_n_append_all_list (self,i,isoform_list):
        ''' splits the join parts of isoform elements'''
        templist=[]
        if 'join' in i:  
            if 'complement' in i:  
                strand='neg'
                loc=re.search(r'complement\(join\((.*?)\)',i)
                if loc: #to avoid errors when strange formats are met (elements with strange format will be skipped)
                    location=loc.group(1)
                    locatn=re.findall(r'([0-9\:]+)', location)
                    for val in locatn:
                        value=self.splitter(val)  
                        templist.append(value)
                else:
                    print(f'the element with coordinates, {i} was skipped due to being strangly formatted')
                
            else:
                strand='pos'
                loc=re.search(r'join\((.*?)\)',i)   
                if loc: #to avoid errors when strange formats are met (elements with strange format will be skipped)                         
                    location=loc.group(1)
                    locatn=re.findall(r'([0-9\:]+)', location)
                    for val in locatn:
                        value=self.splitter(val)
                        templist.append(value)
                else:
                    print(f'the element with coordinates, {i} was skipped due to being strangly formatted')
                    
        else:
            if 'complement' in i:  
                strand='neg'
                loc=re.search(r'complement\((.*?)\)',i)
                if loc: #to avoid errors when strange formats are met (elements with strange format will be skipped)
                    location=loc.group(1)
                    locatn=re.findall(r'([0-9\:]+)', location)
                    for val in locatn:
                        value=self.splitter(val)  
                        templist.append(value)
                else:
                    print(f'the element with coordinates, {i} was skipped due to being strangly formatted')
                    
            else:                
                strand='pos'
                locatn=re.findall(r'([0-9\:]+)', i)   
                if locatn:
                    for val in locatn:
                        value=self.splitter(val)
                        templist.append(value)
                else:
                    print(f'the element with coordinates, {i} was skipped due to being strangly formatted')
        
        isoform_list.append(templist)
        return (strand)
    
    def split_n_append_all (self,i,element_dict):
        ''' splits the join parts of other elements'''
        if 'join' in i:  
            if 'complement' in i:  
                strand='neg'
                loc=re.search(r'complement\(join\((.*?)\)',i)
                if loc: #to avoid errors when strange formats are met (elements with strange format will be skipped)
                    location=loc.group(1)
                    locatn=re.findall(r'([0-9\:]+)', location)
                    for val in locatn:
                        value=self.splitter(val)  
                        element_dict.append(value)
                else:
                    print(f'the element with coordinates, {i} was skipped due to being strangly formatted')
                    
            else:
                strand='pos'
                loc=re.search(r'join\((.*?)\)',i) 
                if loc: #to avoid errors when strange formats are met (elements with strange format will be skipped)                           
                    location=loc.group(1)
                    locatn=re.findall(r'([0-9\:]+)', location)
                    for val in locatn:
                        value=self.splitter(val)
                        element_dict.append(value)
                else:
                    print(f'the element with coordinates, {i} was skipped due to being strangly formatted')
        else:
            if 'complement' in i:  
                strand='neg'
                loc=re.search(r'complement\((.*?)\)',i)
                if loc: #to avoid errors when strange formats are met (elements with strange format will be skipped)
                    location=loc.group(1)
                    locatn=re.findall(r'([0-9\:]+)', location)
                    for val in locatn:
                        value=self.splitter(val)  
                        element_dict.append(value)
                else:
                    print(f'the element with coordinates, {i} was skipped due to being strangly formatted')
            else:                
                strand='pos'
                locatn=re.findall(r'([0-9\:]+)', i) 
                if locatn:
                    for val in locatn:
                        value=self.splitter(val)
                        element_dict.append(value)
                else:
                    print(f'the element with coordinates, {i} was skipped due to being strangly formatted')
        
        return (strand)
    
    
    
    def introns (self,cds_list,i,list_dict,strand):
        '''creates the intron cordinates'''
        intron=(cds_list[i][1])+1,(cds_list[i+1][0])-1
        if strand=='pos':
            list_dict['positives'].append(intron)
        elif strand=='neg':
            list_dict['negatives'].append(intron)
    
    def cordinates (self,i,list_dict,strand):
        ''' creates the element cordinates'''
        if strand=='pos':
            list_dict['positives'].append(i)
        elif strand=='neg':
            list_dict['negatives'].append(i)
            
    def cds_selecter (self,seq,list_of_joins,strand):
        ''' This selects one CDS isoform which satisfy the conditions '''
        ''' Conditions:Having only triplets, starting with ATG|GTG,Being the longest '''
        chosen_dict={} #to save the real CDS (having 3x and ATG|GTG)
        chosen=None #to set the variable empty in case no isoform is selected
        if strand == 'pos':
            isoform_count=0 
            for i in list_of_joins:
                if len(i)!=0:
                    start_codon=seq[(i[0][0])-1:(i[0][0])+2] #the first three bases
                    if start_codon == 'atg' or start_codon =='gtg':
                        count=0
                        tot_len=0 #to save the lenght of the isoform
                        for j in i:
                            short_len=j[1]-(j[0]-1) #the length of one exon
                            tot_len+=short_len
                            count=count+1
                        if tot_len%3==0:
                            chosen_dict[str(isoform_count)]=tot_len
                    isoform_count=isoform_count+1 #update the isoform number
            if chosen_dict: #to make sure that the dict is not empty
                max_key = int(max(chosen_dict, key=lambda k: chosen_dict[k])) #select the key with the highest value and get the int of it
                chosen=list_of_joins[max_key] #chosen isoform 
                        
        elif strand == 'neg':
            isoform_count=0
            for n in list_of_joins:
                if len(n)!=0:
                    start_codon=(seq[(n[-1][1])-3]+seq[(n[-1][1])-2]+seq[(n[-1][1])-1]) #the first three bases
                    if start_codon == 'cat' or start_codon =='cac':
                        count=0
                        tot_len=0 #to save the lenght of the isoform
                        for j in n:
                            short_len=j[1]-(j[0]-1) #the length of one exon
                            tot_len+=short_len
                            count=count+1
                        if tot_len%3==0:
                            chosen_dict[str(isoform_count)]=tot_len
                    isoform_count=isoform_count+1 #update the isoform number
            if chosen_dict: #to make sure that the dict is not empty
                max_key = int(max(chosen_dict, key=lambda k: chosen_dict[k])) #select the key with the highest value and get the int of it
                chosen=list_of_joins[max_key] #chosen isoform
        return chosen
            
    def nc_selecter (self,seq,list_of_joins):
        ''' This selects one non coding rna isoform which satisfy the conditions '''
        ''' Conditions:Being the longest '''
        chosen_dict={} #to save the real CDS (having 3x and ATG|GTG)
        chosen=None #to set the variable empty in case no isoform is selected
        isoform_count=0 
        for i in list_of_joins:
            if len(i)!=0:
                count=0
                tot_len=0 #to save the length of the isoform
                for j in i:
                    short_len=j[1]-(j[0]-1) #the length of one exon
                    tot_len+=short_len
                    count=count+1
                chosen_dict[str(isoform_count)]=tot_len
                isoform_count=isoform_count+1 #update the isoform number
        if chosen_dict: #to make sure that the dict is not empty
            max_key = int(max(chosen_dict, key=lambda k: chosen_dict[k])) #select the key with the highest value and get the int of it
            chosen=list_of_joins[max_key] #chosen isoform 
                    
        return chosen
    
    def element_counter (self,list_dict,strand):
        if strand == 'pos':
            list_dict['counters']['pos_number']+=1
        elif strand == 'neg':
            list_dict['counters']['neg_number']+=1
            
    def count_cols (self,count_dict):
        col_pos=count_dict['counters']['pos_number']
        col_neg=count_dict['counters']['neg_number']
        col_tot=col_pos+col_neg
        
        return [col_pos,col_neg,col_tot]
    
    def des_sorter (self,list_dict):
        '''sorts in descending order'''
        pos_set = set(list_dict['positives'])
        neg_set = set(list_dict['negatives'])
        
        pos_sorted = sorted(pos_set, key=lambda x: (x[0],x[1]),reverse=True) #orders the list based on the starting position 
        neg_sorted = sorted(neg_set, key=lambda x: (x[0],x[1]),reverse=True)
        
        return [pos_sorted,neg_sorted]
    
    def asc_sorter (self,list_dict):
        '''sorts in ascending order'''
        pos_set = set(list_dict['positives'])
        neg_set = set(list_dict['negatives'])
        
        pos_sorted = sorted(pos_set, key=lambda x: (x[0],x[1])) #orders the list based on the starting position 
        neg_sorted = sorted(neg_set, key=lambda x: (x[0],x[1]))
        
        return [pos_sorted,neg_sorted]
    
    def smallest (self,i,sorted_list,del_list,prev_val,prev_val_id): 
        '''gets the shortest piece in ''gene ins gene'' situation'''
        if prev_val is not None and prev_val[0]<=sorted_list[i][0]  and prev_val[1]>=sorted_list[i][1]:
            del_list.append(prev_val_id)
            prev_val=sorted_list[i]
            prev_val_id=i
        elif prev_val is not None and prev_val[0]>=sorted_list[i][0] and prev_val[1]<=sorted_list[i][1]:
            del_list.append(i)        
        else:
            prev_val = sorted_list[i]
            prev_val_id=i 
        return prev_val,prev_val_id

    def biggest (self,i,sorted_list,del_list,prev_val,prev_val_id):
        '''gets the biggest piece in ''gene in gene'' situation'''
        if prev_val is not None and prev_val[0]>=sorted_list[i][0]  and prev_val[1]<=sorted_list[i][1]:
            del_list.append(prev_val_id)
            prev_val=sorted_list[i]
            prev_val_id=i
        elif prev_val is not None and prev_val[0]<=sorted_list[i][0] and prev_val[1]>=sorted_list[i][1]:
            del_list.append(i)        
        else:
            prev_val = sorted_list[i]
            prev_val_id=i 
        return prev_val,prev_val_id
            
    def overlaps (self,i,sorted_list,list_dict):
        '''gets the overlaps'''
        g1_start, g1_end = sorted_list[i]
        g2_start, g2_end = sorted_list[i+1]
        if g1_start <= g2_end and g2_start <= g1_end:
            list_dict.append((g2_start,g1_end))
            
    def seq (self,strand,i,seque):
        '''converts the sequence in capital or in complement from the original sequence'''
        cdsseq=seque[(i[0]-1):i[1]]
        if strand == 'pos':
            sequence=cdsseq.upper()
        elif strand=='neg':
            reverse=cdsseq[::-1]
            sequence=reverse.replace('a','T').replace('t','A').replace('g','C').replace('c','G').replace('n','N') #reverse complement                                                  
        return sequence        
            
    def countbases(self,seq): #class method to count the bases in each element sequence
        '''counts the bases'''
        A=seq.count('A')
        C=seq.count('C')
        G=seq.count('G')
        T=seq.count('T')
        N=seq.count('N')
        total=A+C+G+T+N
        return [total,A,C,G,T,N]
            
    def freq_bases(self,seq,base_counts): 
        ''' calculates the frequency of bases'''
        #nt=self.countbases(seq)
        nt=base_counts
        if nt[0]==0:
            return ''
        else:
            pA = float(nt[1]/nt[0])
            pC = float(nt[2]/nt[0])
            pG = float(nt[3]/nt[0])
            pT = float(nt[4]/nt[0])
            pN = float(nt[5]/nt[0])
            GC = float ((nt[2]+nt[3])/nt[0])
            total=pA+pC+pG+pT+pN
            return [total,pA,pC,pG,pT,pN,GC]          
    
    def chargaff_ct(self,sequence):
        ''' perform chargaff score CT '''
        counts = Counter(sequence)

        def safe_ratio(x, y):
            if x == 0 or y == 0:
                return 0
            return min(x, y) / max(x, y)

        a_t_ratio = safe_ratio(counts['A'], counts['T'])
        c_g_ratio = safe_ratio(counts['C'], counts['G'])

        CT = (a_t_ratio + c_g_ratio) / 2
        return CT
    
    def chargaff_pf(self,sequence):
        '''chargaff pf'''
        counts = Counter(sequence)
        
        def safe_pf(x,y):
            if x==0 and y==0:
                return 0
            return abs((x-y)/(x+y))
        
        a_t=safe_pf(counts['A'],counts['T'])
        c_g=safe_pf(counts['C'],counts['G'])
        
        PF=a_t+c_g
        
        return PF
     
    def shannon (self,seq,base_counts):
        ''' calculates shanon'''
        #s=seq.replace("N","") #Removes the N in the sequence
        #nt=self.countbases(s)
        nt=base_counts
        if nt[0]==0:
            return ''
        else:
            pA = float(nt[1]/nt[0])
            pC = float(nt[2]/nt[0])
            pG = float(nt[3]/nt[0])
            pT = float(nt[4]/nt[0])
            
            if pA == 0:
                pA=1
            if pC == 0:
                pC=1
            if pG == 0:
                pG=1
            if pT == 0:
                pT=1
            
            return -((pA*math.log2(pA))+(pT*math.log2(pT))+(pC*math.log2(pC))+(pG*math.log2(pG)))
           
    def topology (self,seq): 
        '''calculates topological entropy'''
        s=seq.replace("N","") ##Removes the N in the sequence
        n = len(set(s)) # number of bases, this is usually equal to 4
        length = len(s)
        if len(s) != 0 and n > 1: # if len of s is not 0 and number of bases is more than 1
            logg = math.floor(math.log(length,n))
            neww = (s[:(n**logg + logg)]).lower()
            result_mr = math.log(len(set([neww[i:(i+logg)] for i in range(1, n**logg+1)])),n)/logg #the topological entropy score
        else:
            result_mr = ''
        return result_mr 
    
    def codon_count (self,sequence,codon_dict):
        '''calculates codon counts'''
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]  # Get a triplet of characters
            if codon in codon_dict: 
                codon_dict[codon] += 1 #saves the count of each codon in the dict
    
    def calculations (self,sequence,g_counters):
        '''calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary'''
        nt=self.countbases(sequence)
        freq=self.freq_bases(sequence,nt)
        charg_ct=self.chargaff_ct(sequence)
        charg_pf=self.chargaff_pf(sequence)
        shan=self.shannon(sequence,nt)
        topo=self.topology(sequence)
        ###CALCULATIONS###
        #total bases
        g_counters['totA']=g_counters['totA']+nt[1]
        g_counters['totC']=g_counters['totC']+nt[2]
        g_counters['totG']=g_counters['totG']+nt[3]
        g_counters['totT']=g_counters['totT']+nt[4]
        g_counters['totN']=g_counters['totN']+nt[5]
        g_counters['totBase']=g_counters['totBase']+nt[0]
       
        #base frequencies
        if type(freq) is str:
            pass
        else:   
            g_counters['freqA']=g_counters['freqA']+freq[1]
            g_counters['freqC']=g_counters['freqC']+freq[2]
            g_counters['freqG']=g_counters['freqG']+freq[3]
            g_counters['freqT']=g_counters['freqT']+freq[4]
            g_counters['freqN']=g_counters['freqN']+freq[5]
            g_counters['freqGC']=g_counters['freqGC']+freq[6]
            g_counters['freq_counter']=g_counters['freq_counter']+1
        
        #Chargaff scores
        g_counters['tot_PF_charg']=g_counters['tot_PF_charg']+charg_pf
        g_counters['tot_CT_charg']=g_counters['tot_CT_charg']+charg_ct
        g_counters['chargaff_counter']=g_counters['chargaff_counter']+1
        
        #Shannon
        if type(shan) is str:
            pass
        else:
            g_counters['tot_shan']=g_counters['tot_shan']+shan
            g_counters['shan_counter']=g_counters['shan_counter']+1
        
        if type(topo) is str:
            pass
        else:
            g_counters['tot_topo']=g_counters['tot_topo']+topo
            g_counters['topo_counter']=g_counters['topo_counter']+1
            
        return nt[5]
        
            
    def outputs (self,g_counters):
        '''saves the output in relevent column name'''
        column1=g_counters['totA']
        column2=g_counters['totT']
        column3=g_counters['totC']
        column4=g_counters['totG']
        column5=g_counters['totN']
        column6=g_counters['totBase']
            
        try:    #Mean frequencies of each base in all the genes 
            column7='{:.2f}'.format(g_counters['freqA']/g_counters['freq_counter'])
            column8='{:.2f}'.format(g_counters['freqT']/g_counters['freq_counter'])
            column9='{:.2f}'.format(g_counters['freqC']/g_counters['freq_counter'])
            column10='{:.2f}'.format(g_counters['freqG']/g_counters['freq_counter'])
            column11='{:.2f}'.format(g_counters['freqN']/g_counters['freq_counter'])
            column12='{:.2f}'.format(g_counters['freqGC']/g_counters['freq_counter'])
        except ZeroDivisionError:    
            column7=''
            column8=''
            column9=''
            column10=''
            column11=''
            column12=''
        
        try:    #Mean Topological score all the genes in the locus
            column13='{:.4f}'.format(g_counters['tot_topo']/g_counters['topo_counter'])
        except ZeroDivisionError:    
            column13=''
            
        try:    #Mean Chargaff score all the genes in the locus
            column14='{:.4f}'.format(g_counters['tot_PF_charg']/g_counters['chargaff_counter'])
            column15='{:.4f}'.format(g_counters['tot_CT_charg']/g_counters['chargaff_counter'])
        except ZeroDivisionError:    
            column14=''
            column15=''
        
        try:    #Mean Shanon score all the genes in the locus
            column16='{:.4f}'.format(g_counters['tot_shan']/g_counters['shan_counter'])
        except ZeroDivisionError:    
            column16=''
            
        return [column1,column2,column3,column4,column5,column6,column7,column8,column9,column10,column11,column12,column13,column14,column15,column16]
    
    def codon_outputs (self,cod_counter):
        '''saves the output of codon counts in relevent column name'''
        column1=cod_counter["ATG"]
        column2=cod_counter["AAG"]
        column3=cod_counter["GTA"]
        column4=cod_counter["ACT"]
        column5=cod_counter["GCA"]
        column6=cod_counter["GAG"]
        column7=cod_counter["GCT"]
        column8=cod_counter["ATT"]
        column9=cod_counter["TCC"]
        column10=cod_counter["TGG"]
        column11=cod_counter["AAT"]
        column12=cod_counter["GAA"]
        column13=cod_counter["TCA"]
        column14=cod_counter["ACG"]
        column15=cod_counter["AGT"]
        column16=cod_counter["AAC"]
        column17=cod_counter["TCT"]
        column18=cod_counter["GTG"]
        column19=cod_counter["TTC"]
        column20=cod_counter["TTT"]
        column21=cod_counter["CTG"]
        column22=cod_counter["GGT"]
        column23=cod_counter["CTC"]
        column24=cod_counter["GAT"]
        column25=cod_counter["CAG"]
        column26=cod_counter["ACC"]
        column27=cod_counter["CTA"]
        column28=cod_counter["TTG"]
        column29=cod_counter["TAT"]
        column30=cod_counter["GGA"]
        column31=cod_counter["ATC"]
        column32=cod_counter["CTT"]
        column33=cod_counter["GTC"]
        column34=cod_counter["ATA"]
        column35=cod_counter["ACA"]
        column36=cod_counter["GAC"]
        column37=cod_counter["CAC"]
        column38=cod_counter["CCC"]
        column39=cod_counter["TAC"]
        column40=cod_counter["GCC"]
        column41=cod_counter["AGC"]
        column42=cod_counter["CGC"]
        column43=cod_counter["AAA"]
        column44=cod_counter["GGC"]
        column45=cod_counter["TGC"]
        column46=cod_counter["GTT"]
        column47=cod_counter["GGG"]
        column48=cod_counter["AGA"]
        column49=cod_counter["TGT"]
        column50=cod_counter["CAT"]
        column51=cod_counter["TCG"]
        column52=cod_counter["GCG"]
        column53=cod_counter["TTA"]
        column54=cod_counter["CCT"]
        column55=cod_counter["AGG"]
        column56=cod_counter["CCA"]
        column57=cod_counter["TAG"]
        column58=cod_counter["CGG"]
        column59=cod_counter["CCG"]
        column60=cod_counter["CGT"]
        column61=cod_counter["CAA"]
        column62=cod_counter["TAA"]
        column63=cod_counter["TGA"]
        column64=cod_counter["CGA"]
        
        return [column1, column2, column3, column4, column5, column6, column7, column8, column9, column10, column11, column12, column13, column14, column15, column16, column17, column18, column19, column20, column21, column22, column23, column24, column25, column26, column27, column28, column29, column30, column31, column32, column33, column34, column35, column36, column37, column38, column39, column40, column41, column42, column43, column44, column45, column46, column47, column48, column49, column50, column51, column52, column53, column54, column55, column56, column57, column58, column59, column60, column61, column62, column63, column64]
 
    
    def run(self,analysis):
        gbff_file_path = open(self.input_path)        
        
        with open(self.output_path,'w') as outputfile:
            outputfile=csv.writer(outputfile)
            outputfile.writerow(['Class','Organism','Taxon','Assembly','Locus_ID','Version','Definition',
                                 'bp_chromo_A','bp_chromo_T','bp_chromo_C','bp_chromo_G','bp_chromo_N','bp_chromo_tot','fr_chromo_A','fr_chromo_T','fr_chromo_C','fr_chromo_G',
                                 'fr_chromo_N','GC_chromo','topo_entropy_chromo','chargaff_pf_chromo','chargaff_ct_chromo','shannon_chromo',
                                 'n_gene_pos','n_gene_neg','n_gene_tot','bp_gene_A','bp_gene_T','bp_gene_C','bp_gene_G','bp_gene_N','bp_gene_tot','fr_gene_A','fr_gene_T','fr_gene_C',
                                 'fr_gene_G','fr_gene_N','GC_gene','topo_entropy_gene','chargaff_pf_gene','chargaff_ct_gene','shannon_gene','bp_gene_overlap_tot',
                                 'n_cds_pos','n_cds_neg','n_cds_tot','bp_cds_A','bp_cds_T','bp_cds_C','bp_cds_G','bp_cds_N','bp_cds_tot','fr_cds_A','fr_cds_T','fr_cds_C','fr_cds_G',
                                 'fr_cds_N','GC_cds','topo_entropy_cds','chargaff_pf_cds','chargaff_ct_cds','shannon_cds','bp_cds_overlap_tot',
                                 'bp_cds_intron_A','bp_cds_intron_T','bp_cds_intron_C','bp_cds_intron_G','bp_cds_intron_N','bp_cds_intron_tot','fr_cds_intron_A','fr_cds_intron_T',
                                 'fr_cds_intron_C','fr_cds_intron_G','fr_cds_intron_N','GC_cds_intron','topo_entropy_cds_intron','chargaff_pf_cds_intron','chargaff_ct_cds_intron','shannon_cds_intron','bp_cds_intron_overlap_tot',
                                 'n_ncRNA_pos','n_ncRNA_neg','n_ncRNA_tot','bp_ncRNA_A','bp_ncRNA_T','bp_ncRNA_C','bp_ncRNA_G','bp_ncRNA_N','bp_ncRNA_tot','fr_ncRNA_A','fr_ncRNA_T',
                                 'fr_ncRNA_C','fr_ncRNA_G','fr_ncRNA_N','GC_ncRNA','topo_entropy_ncRNA','chargaff_pf_ncRNA','chargaff_ct_ncRNA','shannon_ncRNA','bp_ncRNA_overlap_tot',
                                 'bp_nc_intron_A','bp_nc_intron_T','bp_nc_intron_C','bp_nc_intron_G','bp_nc_intron_N','bp_nc_intron_tot','fr_nc_intron_A','fr_nc_intron_T','fr_nc_intron_C',
                                 'fr_nc_intron_G','fr_nc_intron_N','GC_nc_intron','topo_entropy_nc_intron','chargaff_pf_nc_intron','chargaff_ct_nc_intron','shannon_nc_intron','bp_nc_intron_overlap_tot',
                                 'n_tRNA_pos','n_tRNA_neg','n_tRNA_tot','bp_tRNA_A','bp_tRNA_T','bp_tRNA_C','bp_tRNA_G','bp_tRNA_N','bp_tRNA_tot','fr_tRNA_A','fr_tRNA_T','fr_tRNA_C','fr_tRNA_G','fr_tRNA_N',
                                 'GC_tRNA','topo_entropy_tRNA','chargaff_pf_tRNA','chargaff_ct_tRNA','shannon_tRNA','bp_tRNA_overlap_tot',
                                 'n_rRNA_pos','n_rRNA_neg','n_rRNA_tot','bp_rRNA_A','bp_rRNA_T','bp_rRNA_C','bp_rRNA_G','bp_rRNA_N','bp_rRNA_tot','fr_rRNA_A','fr_rRNA_T','fr_rRNA_C',
                                 'fr_rRNA_G','fr_rRNA_N','GC_rRNA','topo_entropy_rRNA','chargaff_pf_rRNA','chargaff_ct_rRNA','shannon_rRNA','bp_rRNA_overlap_tot',
                                 'ATG','AAG','GTA','ACT','GCA','GAG',
                                 'GCT','ATT','TCC','TGG','AAT','GAA','TCA','ACG','AGT','AAC','TCT','GTG','TTC','TTT','CTG','GGT','CTC','GAT','CAG','ACC','CTA','TTG','TAT','GGA',
                                 'ATC','CTT','GTC','ATA','ACA','GAC','CAC','CCC','TAC','GCC','AGC','CGC','AAA','GGC','TGC','GTT','GGG','AGA','TGT','CAT','TCG','GCG','TTA','CCT','AGG',
                                 'CCA','TAG','CGG','CCG','CGT','CAA','TAA','TGA','CGA'])
            
            #gbrp_obj=GBRAP()
            
            #to save the counts of whole genome
            genome={'chromodicts':{'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}},
            'genedicts':{'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                        'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}},
            'cdsdicts':{'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                      'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                      'nucl_codon_counter':{'ATG': 0, 'AAG': 0, 'GTA': 0, 'ACT': 0, 'GCA': 0, 'GAG': 0, 'GCT': 0, 'ATT': 0, 'TCC': 0, 'TGG': 0, 'AAT': 0, 'GAA': 0, 'TCA': 0, 'ACG': 0, 
                                         'AGT': 0, 'AAC': 0, 'TCT': 0, 'GTG': 0, 'TTC': 0, 'TTT': 0, 'CTG': 0, 'GGT': 0, 'CTC': 0, 'GAT': 0, 'CAG': 0, 'ACC': 0, 'CTA': 0, 'TTG': 0, 
                                         'TAT': 0, 'GGA': 0, 'ATC': 0, 'CTT': 0, 'GTC': 0, 'ATA': 0, 'ACA': 0, 'GAC': 0, 'CAC': 0, 'CCC': 0, 'TAC': 0, 'GCC': 0, 'AGC': 0, 'CGC': 0, 
                                         'AAA': 0, 'GGC': 0, 'TGC': 0, 'GTT': 0, 'GGG': 0, 'AGA': 0, 'TGT': 0, 'CAT': 0, 'TCG': 0, 'GCG': 0, 'TTA': 0, 'CCT': 0, 'AGG': 0, 'CCA': 0, 
                                         'TAG': 0, 'CGG': 0, 'CCG': 0, 'CGT': 0, 'CAA': 0, 'TAA': 0, 'TGA': 0, 'CGA': 0}},
            'cds_introndicts':{'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                      'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}},
            'ncRNAdicts':{'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                      'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}},
            'nc_introndicts':{'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                      'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}},
            'tRNAdicts':{'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                      'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}},
            'rRNAdicts':{'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                      'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                              'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}}}
        

            locus_counter=1
            skipped=0
            analysed_chromo=0
            for result in re.findall('LOCUS(.*?)//\n',gbff_file_path.read(), re.S): #takes all the info in each locus
                if locus_counter==1: #takes the domain name of the organism
                    org=re.search(r'ORGANISM\s{2}(.*)\n\s*([a-zA-Z;\s]+).',result) #takes the organism row
                    taxa=org.group(2).replace(' ','').replace('\n','') #takes the taxonomy
                    organism=org.group(1).replace(' ','_') #takes the organism name
                    print(f'\nGBRAP Analysis started for {organism}\n--')
                    
                    category_dict={'Mammal':['Mammalia'],'Amphibia':['Amphibia'],'Reptile':['Lepidosauria','Crocodylia','Testudines'],'Bird':['Aves,'],
                                 'Fish':['Dipnomorpha','Actinopterygii','Chondrichthyes','Cyclostomata','Coelacanthiformes'],'Plant':['Viridiplantae','Rhodophyta'],
                                 'Protozoa':['Sar','Discoba','Amoebozoa','Metamonada','Cryptophyceae'],'Fungi':['Fungi'],'Bacteria':['Bacteria'],'Archaea':['Archaea'],'Virus':['Viruses'],
                                 'Invertebrate':['Arthropoda','Cephalochordata','Tunicata','Nematoda','Cnidaria','Mollusca','Echinodermata','Hemichordata','Platyhelminthes','Ctenophora','Porifera','Brachiopoda','Choanoflagellata','Scalidophora','Xenacoelomorpha']}
                   
                    category = 'Not_classified'
                    for key, values in category_dict.items():
                        if any(value in taxa.split(';') for value in values):  #taxa.split to exactly match values with the complete words in taxa (ex in viral the taxa level, Mammalianorthoreovirus contains Mammalian) 
                            category=key
                            break
                        
                # if locus_counter==2: #to stop the code after the desired number of chromosomes
                #     break   
                
                definition=re.search(r'DEFINITION(.*?)ACCESSION',result,re.S) #takes the definition of the locus
                col2=definition.group(1).strip().replace('  ','').replace('\n',' ').replace(',',';') #prints the definition
                #Any sequence which contains either of the followings words in the DEFINITION will be removed from the analysis.
                skiplocus=['scaffold','unlocalized','contig','unplaced','patch','unknown'] 
                if any (word in definition.group(1).lower() for word in skiplocus):
                    skipped=skipped+1
                
                else:
                    analysed_chromo+=1
                    
                    locusline=re.search(r'\s{7}([a-zA-Z_0-9]+)\s+[0-9]+\sbp', result) #takes the locus name 
                    col1=locusline.group(1) #prints the locus name
                    
                    ver=re.search(r'VERSION\s{5}(.*?)\n',result) 
                    if ver is not None:
                        version=ver.group(1)
                    else:
                        version='NA'
                        
                    ass=re.search(r'\s{12}Assembly:\s(.*?)\n',result)
                    if ass is not None:
                        assembly=ass.group(1)
                    else:
                        assembly='NA'
                    
                    fullseq=re.findall('ORIGIN\s*\n\s*(.*)',result,re.S) #ignore this warning. it is only due to new python version
                    se=[g.replace(' ','').replace('\n','').replace('\t','') for g in fullseq] #replaces the \s \n and \t
                    seq=re.sub("[^a-zA-Z]","",str(se).lower())   #the whole sequene of the locus ready to use removes the numbers
                    features=re.search(r'FEATURES\s*Location/Qualifiers(.*?)ORIGIN',result,re.S) #takes all the details of the annotations
                    chromodicts={'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                        'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}}
                    genedicts={'lists':{'positives':[],'negatives':[],'pos_overlaps':[],'neg_overlaps':[]},
                                'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                        'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                                'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                        'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}}
                    cdsdicts={'lists':{'positives':[],'negatives':[],'pos_overlaps':[],'neg_overlaps':[]},
                              'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                              'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                              'nucl_codon_counter':{'ATG': 0, 'AAG': 0, 'GTA': 0, 'ACT': 0, 'GCA': 0, 'GAG': 0, 'GCT': 0, 'ATT': 0, 'TCC': 0, 'TGG': 0, 'AAT': 0, 'GAA': 0, 'TCA': 0, 'ACG': 0, 
                                                 'AGT': 0, 'AAC': 0, 'TCT': 0, 'GTG': 0, 'TTC': 0, 'TTT': 0, 'CTG': 0, 'GGT': 0, 'CTC': 0, 'GAT': 0, 'CAG': 0, 'ACC': 0, 'CTA': 0, 'TTG': 0, 
                                                 'TAT': 0, 'GGA': 0, 'ATC': 0, 'CTT': 0, 'GTC': 0, 'ATA': 0, 'ACA': 0, 'GAC': 0, 'CAC': 0, 'CCC': 0, 'TAC': 0, 'GCC': 0, 'AGC': 0, 'CGC': 0, 
                                                 'AAA': 0, 'GGC': 0, 'TGC': 0, 'GTT': 0, 'GGG': 0, 'AGA': 0, 'TGT': 0, 'CAT': 0, 'TCG': 0, 'GCG': 0, 'TTA': 0, 'CCT': 0, 'AGG': 0, 'CCA': 0, 
                                                 'TAG': 0, 'CGG': 0, 'CCG': 0, 'CGT': 0, 'CAA': 0, 'TAA': 0, 'TGA': 0, 'CGA': 0},
                              'overlap_nucl_codon_counter':{'ATG': 0, 'AAG': 0, 'GTA': 0, 'ACT': 0, 'GCA': 0, 'GAG': 0, 'GCT': 0, 'ATT': 0, 'TCC': 0, 'TGG': 0, 'AAT': 0, 'GAA': 0, 'TCA': 0, 'ACG': 0, 
                                                 'AGT': 0, 'AAC': 0, 'TCT': 0, 'GTG': 0, 'TTC': 0, 'TTT': 0, 'CTG': 0, 'GGT': 0, 'CTC': 0, 'GAT': 0, 'CAG': 0, 'ACC': 0, 'CTA': 0, 'TTG': 0, 
                                                 'TAT': 0, 'GGA': 0, 'ATC': 0, 'CTT': 0, 'GTC': 0, 'ATA': 0, 'ACA': 0, 'GAC': 0, 'CAC': 0, 'CCC': 0, 'TAC': 0, 'GCC': 0, 'AGC': 0, 'CGC': 0, 
                                                 'AAA': 0, 'GGC': 0, 'TGC': 0, 'GTT': 0, 'GGG': 0, 'AGA': 0, 'TGT': 0, 'CAT': 0, 'TCG': 0, 'GCG': 0, 'TTA': 0, 'CCT': 0, 'AGG': 0, 'CCA': 0, 
                                                 'TAG': 0, 'CGG': 0, 'CCG': 0, 'CGT': 0, 'CAA': 0, 'TAA': 0, 'TGA': 0, 'CGA': 0}}
                    cds_introndicts={'lists':{'positives':[],'negatives':[],'pos_overlaps':[],'neg_overlaps':[]},
                              'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                              'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}}
                    ncRNAdicts={'lists':{'positives':[],'negatives':[],'pos_overlaps':[],'neg_overlaps':[]},
                              'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                              'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}}
                    nc_introndicts={'lists':{'positives':[],'negatives':[],'pos_overlaps':[],'neg_overlaps':[]},
                              'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                              'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}}
                    tRNAdicts={'lists':{'positives':[],'negatives':[],'pos_overlaps':[],'neg_overlaps':[]},
                              'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                              'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}}
                    rRNAdicts={'lists':{'positives':[],'negatives':[],'pos_overlaps':[],'neg_overlaps':[]},
                              'counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0},
                              'overlap_counters':{'pos_number':0,'neg_number':0,'totA':0,'totC':0,'totG':0,'totT':0,'totN':0,'totBase':0,'freqA':0,'freqC':0,'freqG':0,'freqT':0,'freqN':0,'freqGC':0,'freq_counter':0,
                                      'tot_PF_charg':0,'tot_CT_charg':0,'chargaff_counter':0,'tot_shan':0,'shan_counter':0,'tot_topo':0,'topo_counter':0}}
                    for each in re.findall(r'\s{5}gene\s{12}.*?(?=\s{5}gene\s{12}|$)',features.group(1),re.S): #for each gene
                        if '/pseudo' in each: ###maybe have to use a different word (ex 'LOC[0-9]+')####
                            pass
                        else:
                            note=re.search(r'\s{21}/note="([^"]+)"', each) #checks for the note 
                            if note: 
                                prot_info=note.group(1) #saves the info in note
                            else:
                                prot_info='' #to avoid erros if note does not exist
                                
                            gene=re.findall(r'\s{5}gene\s{12}(.*?)\/',each,re.S)#a list of all the positions of the desired element
                            clean_gene =[g.replace(' ','').replace('\n','').replace('\t','').replace('..',':')
                                          .replace('<','').replace('>','') for g in gene] #the list of the posit
                            genedict=[]
                            for i in clean_gene:
                                strand=self.split_n_append_all(i,genedict) #splits the joined parts in one isoform
                            if genedict != None:
                                self.element_counter(genedicts,strand)
                            for v in genedict:   
                                self.cordinates(v,genedicts['lists'],strand)
                            
                            cds=re.findall(r'\s{5}CDS\s{13}(.*?)\/',each,re.S)#a list of all the positions of the desired element
                            clean_cds =[g.replace(' ','').replace('\n','').replace('\t','').replace('..',':')
                                          .replace('<','').replace('>','') for g in cds] #the list of the posit
                            if 'ribosomal' in prot_info:
                                cds_type="Ribo"
                            else:
                                cds_type="Nucl"
                            #exonlist=[] #a dict is needed when taking introns of each isoform
                            #tempdict={}
                            isoforms=[]
                            for i in clean_cds:
                                strand=self.split_n_append_all_list(i,isoforms) #splits the joined parts in one isoform    
                            chosen_CDS=self.cds_selecter(seq,isoforms, strand)
                            if chosen_CDS != None:
                                self.element_counter(cdsdicts,strand)
                                for v in chosen_CDS: 
                                    self.cordinates(v,cdsdicts['lists'],strand)
                                for i in range (len(chosen_CDS)-1): #introns from CDS
                                    self.introns(chosen_CDS, i, cds_introndicts['lists'],strand)
                               
             
                            ncrna=re.findall(r'\s{5}ncRNA\s{11}(.*?)\/',each,re.S)#a list of all the positions of the desired element
                            clean_ncrna =[g.replace(' ','').replace('\n','').replace('\t','').replace('..',':')
                                          .replace('<','').replace('>','') for g in ncrna] #the list of the posit
                    
                            
                            ncdict=[] #a dict is needed when taking introns of each isoform
                            for i in clean_ncrna:
                                strand=self.split_n_append_all_list(i,ncdict) #splits the joined parts in one isoform
                            
                            chosen_nc=self.nc_selecter(seq,ncdict)
                            if chosen_nc != None:
                                self.element_counter(ncRNAdicts,strand)
                                for v in chosen_nc: 
                                    self.cordinates(v,ncRNAdicts['lists'],strand)
                                for i in range (len(chosen_nc)-1): #introns from CDS
                                    self.introns(chosen_nc, i, nc_introndicts['lists'],strand)
                            
                            
                            trna=re.findall(r'\s{5}tRNA\s{12}(.*?)\/',each,re.S)#a list of all the positions of the desired element
                            clean_trna =[g.replace(' ','').replace('\n','').replace('\t','').replace('..',':')
                                          .replace('<','').replace('>','') for g in trna] #the list of the posit
                            trnadict=[] #a dict is needed when taking introns of each isoform
                            for i in clean_trna:
                                strand=self.split_n_append_all_list(i,trnadict) #splits the joined parts in one isoform
                            chosen_trna=self.nc_selecter(seq,trnadict)
                            if chosen_trna != None:
                                self.element_counter(tRNAdicts,strand)
                                for v in chosen_trna:                                                           
                                    self.cordinates(v,tRNAdicts['lists'],strand)    
                                    
                            rrna=re.findall(r'\s{5}rRNA\s{12}(.*?)\/',each,re.S)#a list of all the positions of the desired element
                            clean_rrna =[g.replace(' ','').replace('\n','').replace('\t','').replace('..',':')
                                          .replace('<','').replace('>','') for g in rrna] #the list of the posit
                            rrnadict=[] #a dict is needed when taking introns of each isoform
                            for i in clean_rrna:
                                strand=self.split_n_append_all_list(i,rrnadict) #splits the joined parts in one isoform
                            chosen_rrna=self.nc_selecter(seq,rrnadict)
                            if chosen_rrna != None:
                                self.element_counter(rRNAdicts,strand)
                                for v in chosen_rrna:                                                           
                                    self.cordinates(v,rRNAdicts['lists'],strand)    
                            
                    #element count outputs
                    gen_counts=self.count_cols(genedicts)
                    cds_counts=self.count_cols(cdsdicts)
                    nc_counts=self.count_cols(ncRNAdicts)
                    tr_counts=self.count_cols(tRNAdicts)
                    rr_counts=self.count_cols(rRNAdicts)
                    
                    #Calculations and Outputs starts from here#
                    ####Whole chromosome####
                    sequence=seq.upper()
                    self.calculations (sequence,chromodicts['counters'])
                    chromo_cols=self.outputs(chromodicts['counters'])
                    
                    #####GENE####
                    gene_sorts=self.asc_sorter(genedicts['lists'])#sorts the list and removes duplicates     
                    #pos gene filtering
                    gen_pos_to_del = []
                    gen_pos_prev_value = None
                    gen_pos_prev_value_index=None
                    for i in range(len(gene_sorts[0])): #deletes the smaller ones, keeping the biggest gene
                        gen_pos_prev_value,gen_pos_prev_value_index=self.biggest(i, gene_sorts[0], gen_pos_to_del,gen_pos_prev_value,gen_pos_prev_value_index)                 
                    gen_pos_to_del=list(set(gen_pos_to_del))
                    gen_pos_to_del.sort(reverse=True) 
                    for i in gen_pos_to_del:     
                        del(gene_sorts[0][i])
                    #pos overlaps
                    for i in range(len(gene_sorts[0])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,gene_sorts[0],genedicts['lists']['pos_overlaps'])
                                
                    #Neg gene filtering
                    gen_neg_to_del = []
                    gen_neg_prev_value = None
                    gen_neg_prev_value_index=None
                    for i in range(len(gene_sorts[1])): #deletes the bigger ones, keeping the smallest intron    
                        gen_neg_prev_value,gen_neg_prev_value_index=self.biggest(i, gene_sorts[1], gen_neg_to_del,gen_neg_prev_value,gen_neg_prev_value_index)  
                    gen_neg_to_del=list(set(gen_neg_to_del))
                    gen_neg_to_del.sort(reverse=True)
                    for i in gen_neg_to_del:     
                        del(gene_sorts[1][i])   
                    #neg overlaps
                    for i in range(len(gene_sorts[1])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,gene_sorts[1],genedicts['lists']['neg_overlaps'])
                    
                    #Calculations and outputs
                    for i in gene_sorts[0]: #for positive genes  
                            sequence=self.seq('pos',i,seq)              
                            self.calculations (sequence,genedicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in gene_sorts[1]: #for negative genes
                            sequence=self.seq('neg',i,seq)              
                            self.calculations (sequence,genedicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    gen_cols=self.outputs(genedicts['counters'])
                    
                    for i in genedicts['lists']['pos_overlaps']: #for positive gene overlaps
                            sequence=self.seq('pos',i,seq)                               
                            self.calculations (sequence,genedicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in genedicts['lists']['neg_overlaps']: #for negative gene overlaps  
                            sequence=self.seq('neg',i,seq)                             
                            self.calculations (sequence,genedicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    gen_over_cols=self.outputs(genedicts['overlap_counters'])
                    
                    #####CDS####
                    cds_sorts=self.asc_sorter(cdsdicts['lists'])#sorts the list and removes duplicates     
                    #pos exons filtering
                    pos_cds_to_del = []
                    pos_prev_value = None
                    pos_prev_value_index=None
                    for i in range(len(cds_sorts[0])): #deletes the bigger ones, keeping the smallest intron
                        pos_prev_value,pos_prev_value_index=self.biggest(i, cds_sorts[0], pos_cds_to_del,pos_prev_value,pos_prev_value_index)                 
                    pos_cds_to_del=list(set(pos_cds_to_del))
                    pos_cds_to_del.sort(reverse=True) 
                    for i in pos_cds_to_del:     
                        del(cds_sorts[0][i])
                    #pos overlaps
                    for i in range(len(cds_sorts[0])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,cds_sorts[0],cdsdicts['lists']['pos_overlaps'])
                                
                    #Neg exons filtering
                    neg_cds_to_del = []
                    neg_prev_value = None
                    neg_prev_value_index=None
                    for i in range(len(cds_sorts[1])): #deletes the bigger ones, keeping the smallest intron    
                        neg_prev_value,neg_prev_value_index=self.biggest(i, cds_sorts[1], neg_cds_to_del,neg_prev_value,neg_prev_value_index)  
                    neg_cds_to_del=list(set(neg_cds_to_del))
                    neg_cds_to_del.sort(reverse=True)
                    for i in neg_cds_to_del:     
                        del(cds_sorts[1][i])   
                    #neg overlaps
                    for i in range(len(cds_sorts[1])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,cds_sorts[1],cdsdicts['lists']['neg_overlaps'])
                    
                    #Calculations and outputs
                    for i in cds_sorts[0]: #for positive exons  
                            sequence=self.seq('pos',i,seq) 
                            self.calculations (sequence,cdsdicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                            self.codon_count(sequence, cdsdicts['nucl_codon_counter']) #calculates the codon counts and saves them in the dict
                    for i in cds_sorts[1]: #for negative exons
                            sequence=self.seq('neg',i,seq)              
                            self.calculations (sequence,cdsdicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                            self.codon_count(sequence, cdsdicts['nucl_codon_counter'])
                    cds_cols=self.outputs(cdsdicts['counters'])
                    codon_cols=self.codon_outputs(cdsdicts['nucl_codon_counter'])
                    #overlaps
                    for i in cdsdicts['lists']['pos_overlaps']: #for positive exon overlaps
                            sequence=self.seq('pos',i,seq)                               
                            self.calculations (sequence,cdsdicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                            self.codon_count(sequence, cdsdicts['overlap_nucl_codon_counter'])
                    for i in cdsdicts['lists']['neg_overlaps']: #for negative exon overlaps  
                            sequence=self.seq('neg',i,seq)                             
                            self.calculations (sequence,cdsdicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                            self.codon_count(sequence, cdsdicts['overlap_nucl_codon_counter'])
                    cds_over_cols=self.outputs(cdsdicts['overlap_counters'])
                    cds_over_codon_cols=self.codon_outputs(cdsdicts['overlap_nucl_codon_counter'])
                    
                    
                    #####CDS INTRONS####
                    cdsintron_sorts=self.des_sorter(cds_introndicts['lists']) #sorts the list and removes duplicates     
                    #deletes the bigger ones, keeping the smallest intron
                    pos_intron_to_del = []
                    pos_prev_value = None
                    pos_prev_value_index=None
                    for i in range(len(cdsintron_sorts[0])): 
                        pos_prev_value,pos_prev_value_index=self.smallest(i, cdsintron_sorts[0], pos_intron_to_del,pos_prev_value,pos_prev_value_index)                 
                    pos_intron_to_del=list(set(pos_intron_to_del))
                    pos_intron_to_del.sort(reverse=True) 
                    for i in pos_intron_to_del:     
                        del(cdsintron_sorts[0][i])
                    #overlaps#
                    for i in range(len(cdsintron_sorts[0])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,cdsintron_sorts[0],cds_introndicts['lists']['pos_overlaps'])
                            
                    #deletes the bigger ones, keeping the smallest intron
                    neg_intron_to_del = []
                    neg_prev_value = None
                    neg_prev_value_index=None
                    for i in range(len(cdsintron_sorts[1])): 
                        neg_prev_value,neg_prev_value_index=self.smallest(i, cdsintron_sorts[1], neg_intron_to_del,neg_prev_value,neg_prev_value_index)                 
                        
                    neg_intron_to_del=list(set(neg_intron_to_del))
                    neg_intron_to_del.sort(reverse=True)
                    for i in neg_intron_to_del:     
                        del(cdsintron_sorts[1][i])   
                    #overlaps
                    for i in range(len(cdsintron_sorts[1])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,cdsintron_sorts[1],cds_introndicts['lists']['neg_overlaps'])
                            
                                
                    for i in cdsintron_sorts[0]: #for positive introns   
                            sequence=self.seq('pos',i,seq)              
                            self.calculations (sequence,cds_introndicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in cdsintron_sorts[1]: #for negative introns
                            sequence=self.seq('neg',i,seq)              
                            self.calculations (sequence,cds_introndicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    CDSint_cols=self.outputs(cds_introndicts['counters'])
                    
                    for i in cds_introndicts['lists']['pos_overlaps']: #for positive intron overlaps   
                            sequence=self.seq('pos',i,seq)                            
                            self.calculations (sequence,cds_introndicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in cds_introndicts['lists']['neg_overlaps']: #for negative intron overlaps     
                            sequence=self.seq('neg',i,seq)                          
                            self.calculations (sequence,cds_introndicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    
                    CDSint_over_cols=self.outputs(cds_introndicts['overlap_counters'])
                    
                    #####ncRNA####
                    nc_sorts=self.asc_sorter(ncRNAdicts['lists'])#sorts the list and removes duplicates     
                    #pos ncRNA filtering
                    nc_pos_to_del = []
                    nc_pos_prev_value = None
                    nc_pos_prev_value_index=None
                    for i in range(len(nc_sorts[0])): #deletes the smaller ones, keeping the biggest 
                        nc_pos_prev_value,nc_pos_prev_value_index=self.biggest(i, nc_sorts[0], nc_pos_to_del,nc_pos_prev_value,nc_pos_prev_value_index)                 
                    nc_pos_to_del=list(set(nc_pos_to_del))
                    nc_pos_to_del.sort(reverse=True) 
                    for i in nc_pos_to_del:     
                        del(nc_sorts[0][i])
                    #pos overlaps
                    for i in range(len(nc_sorts[0])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,nc_sorts[0],ncRNAdicts['lists']['pos_overlaps'])
                                
                    #Neg ncRNA filtering
                    nc_neg_to_del = []
                    nc_neg_prev_value = None
                    nc_neg_prev_value_index=None
                    for i in range(len(nc_sorts[1])): #deletes the smaller ones, keeping the biggest   
                        nc_neg_prev_value,nc_neg_prev_value_index=self.biggest(i, nc_sorts[1], nc_neg_to_del,nc_neg_prev_value,nc_neg_prev_value_index)  
                    nc_neg_to_del=list(set(nc_neg_to_del))
                    nc_neg_to_del.sort(reverse=True)
                    for i in nc_neg_to_del:     
                        del(nc_sorts[1][i])   
                    #neg overlaps
                    for i in range(len(nc_sorts[1])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,nc_sorts[1],ncRNAdicts['lists']['neg_overlaps'])
                    
                    #Calculations and outputs
                    for i in nc_sorts[0]: #for positive ncRNA  
                            sequence=self.seq('pos',i,seq)              
                            self.calculations (sequence,ncRNAdicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in nc_sorts[1]: #for negative ncRNA
                            sequence=self.seq('neg',i,seq)              
                            self.calculations (sequence,ncRNAdicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    nc_cols=self.outputs(ncRNAdicts['counters'])
                    
                    for i in ncRNAdicts['lists']['pos_overlaps']: #for positive ncRNA overlaps
                            sequence=self.seq('pos',i,seq)                               
                            self.calculations (sequence,ncRNAdicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in ncRNAdicts['lists']['neg_overlaps']: #for negative ncRNA overlaps  
                            sequence=self.seq('neg',i,seq)                             
                            self.calculations (sequence,ncRNAdicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    nc_over_cols=self.outputs(ncRNAdicts['overlap_counters'])
                    
                    #####ncRNA INTRONS####
                    ncintron_sorts=self.des_sorter(nc_introndicts['lists']) #sorts the list and removes duplicates    
                    #pos ncRNA intron filtering
                    ncIntro_pos_to_del = []
                    ncIntro_pos_prev_value = None
                    ncIntro_pos_prev_value_index=None
                    for i in range(len(ncintron_sorts[0])): #deletes the bigger ones, keeping the smallest intron
                        ncIntro_pos_prev_value,ncIntro_pos_prev_value_index=self.smallest(i, ncintron_sorts[0], ncIntro_pos_to_del,ncIntro_pos_prev_value,ncIntro_pos_prev_value_index)                 
                    ncIntro_pos_to_del=list(set(ncIntro_pos_to_del))
                    ncIntro_pos_to_del.sort(reverse=True) 
                    for i in ncIntro_pos_to_del:     
                        del(ncintron_sorts[0][i])
                    #overlaps#
                    for i in range(len(ncintron_sorts[0])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,ncintron_sorts[0],nc_introndicts['lists']['pos_overlaps'])
                    
                    #neg ncRNA intron filtering
                    ncIntro_neg_to_del = []
                    ncIntro_neg_prev_value = None
                    ncIntro_neg_prev_value_index=None
                    for i in range(len(ncintron_sorts[1])): #deletes the bigger ones, keeping the smallest intron
                        ncIntro_neg_prev_value,ncIntro_neg_prev_value_index=self.smallest(i, ncintron_sorts[1], ncIntro_neg_to_del,ncIntro_neg_prev_value,ncIntro_neg_prev_value_index)                 
                    ncIntro_neg_to_del=list(set(ncIntro_neg_to_del))
                    ncIntro_neg_to_del.sort(reverse=True)
                    for i in ncIntro_neg_to_del:     
                        del(ncintron_sorts[1][i])   
                    #overlaps
                    for i in range(len(ncintron_sorts[1])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,ncintron_sorts[1],nc_introndicts['lists']['neg_overlaps'])
                            
                    for i in ncintron_sorts[0]: #for positive introns   
                            sequence=self.seq('pos',i,seq)              
                            self.calculations (sequence,nc_introndicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in ncintron_sorts[1]: #for negative introns
                            sequence=self.seq('neg',i,seq)              
                            self.calculations (sequence,nc_introndicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    ncint_cols=self.outputs(nc_introndicts['counters'])
                    
                    for i in nc_introndicts['lists']['pos_overlaps']: #for positive intron overlaps   
                            sequence=self.seq('pos',i,seq)                            
                            self.calculations (sequence,nc_introndicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in nc_introndicts['lists']['neg_overlaps']: #for negative intron overlaps     
                            sequence=self.seq('neg',i,seq)                          
                            self.calculations (sequence,nc_introndicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    
                    ncint_over_cols=self.outputs(nc_introndicts['overlap_counters'])
                    
                    #####tRNA####
                    trna_sorts=self.asc_sorter(tRNAdicts['lists'])#sorts the list and removes duplicates     
                    #pos tRNA filtering
                    tr_pos_to_del = []
                    tr_pos_prev_value = None
                    tr_pos_prev_value_index=None
                    for i in range(len(trna_sorts[0])): #deletes the smaller ones, keeping the biggest 
                        tr_pos_prev_value,tr_pos_prev_value_index=self.biggest(i, trna_sorts[0], tr_pos_to_del,tr_pos_prev_value,tr_pos_prev_value_index)                 
                    tr_pos_to_del=list(set(tr_pos_to_del))
                    tr_pos_to_del.sort(reverse=True) 
                    for i in tr_pos_to_del:     
                        del(trna_sorts[0][i])
                    #pos overlaps
                    for i in range(len(trna_sorts[0])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,trna_sorts[0],tRNAdicts['lists']['pos_overlaps'])
                                
                    #Neg tRNA filtering
                    tr_neg_to_del = []
                    tr_neg_prev_value = None
                    tr_neg_prev_value_index=None
                    for i in range(len(trna_sorts[1])): #deletes the smaller ones, keeping the biggest   
                        tr_neg_prev_value,tr_neg_prev_value_index=self.biggest(i, trna_sorts[1], tr_neg_to_del,tr_neg_prev_value,tr_neg_prev_value_index)  
                    tr_neg_to_del=list(set(tr_neg_to_del))
                    tr_neg_to_del.sort(reverse=True)
                    for i in tr_neg_to_del:     
                        del(trna_sorts[1][i])   
                    #neg overlaps
                    for i in range(len(trna_sorts[1])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,trna_sorts[1],tRNAdicts['lists']['neg_overlaps'])
                    
                    #Calculations and outputs
                    for i in trna_sorts[0]: #for positive tRNA  
                            sequence=self.seq('pos',i,seq)              
                            self.calculations (sequence,tRNAdicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in trna_sorts[1]: #for negative tRNA
                            sequence=self.seq('neg',i,seq)              
                            self.calculations (sequence,tRNAdicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    tr_cols=self.outputs(tRNAdicts['counters'])
                    
                    for i in tRNAdicts['lists']['pos_overlaps']: #for positive tRNA overlaps
                            sequence=self.seq('pos',i,seq)                               
                            self.calculations (sequence,tRNAdicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in tRNAdicts['lists']['neg_overlaps']: #for negative tRNA overlaps  
                            sequence=self.seq('neg',i,seq)                             
                            self.calculations (sequence,tRNAdicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    tr_over_cols=self.outputs(tRNAdicts['overlap_counters'])
                    
                    #####rRNA####
                    rrna_sorts=self.asc_sorter(rRNAdicts['lists'])#sorts the list and removes duplicates     
                    #pos rRNA filtering
                    rr_pos_to_del = []
                    rr_pos_prev_value = None
                    rr_pos_prev_value_index=None
                    for i in range(len(rrna_sorts[0])): #deletes the smaller ones, keeping the biggest 
                        rr_pos_prev_value,rr_pos_prev_value_index=self.biggest(i, rrna_sorts[0], rr_pos_to_del,rr_pos_prev_value,rr_pos_prev_value_index)                 
                    rr_pos_to_del=list(set(rr_pos_to_del))
                    rr_pos_to_del.sort(reverse=True) 
                    for i in rr_pos_to_del:     
                        del(rrna_sorts[0][i])
                    #pos overlaps
                    for i in range(len(rrna_sorts[0])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,rrna_sorts[0],rRNAdicts['lists']['pos_overlaps'])
                                
                    #Neg rRNA filtering
                    rr_neg_to_del = []
                    rr_neg_prev_value = None
                    rr_neg_prev_value_index=None
                    for i in range(len(rrna_sorts[1])): #deletes the smaller ones, keeping the biggest   
                        rr_neg_prev_value,rr_neg_prev_value_index=self.biggest(i, rrna_sorts[1], rr_neg_to_del,rr_neg_prev_value,rr_neg_prev_value_index)  
                    rr_neg_to_del=list(set(rr_neg_to_del))
                    rr_neg_to_del.sort(reverse=True)
                    for i in rr_neg_to_del:     
                        del(rrna_sorts[1][i])   
                    #neg overlaps
                    for i in range(len(rrna_sorts[1])-1): #fills the overlap list with the regions that are overlapping
                        self.overlaps(i,rrna_sorts[1],rRNAdicts['lists']['neg_overlaps'])
                    
                    #Calculations and outputs
                    for i in rrna_sorts[0]: #for positive rRNA  
                            sequence=self.seq('pos',i,seq)              
                            self.calculations (sequence,rRNAdicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in rrna_sorts[1]: #for negative rRNA
                            sequence=self.seq('neg',i,seq)              
                            self.calculations (sequence,rRNAdicts['counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    rr_cols=self.outputs(rRNAdicts['counters'])
                    
                    for i in rRNAdicts['lists']['pos_overlaps']: #for positive rRNA overlaps
                            sequence=self.seq('pos',i,seq)                               
                            self.calculations (sequence,rRNAdicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    for i in rRNAdicts['lists']['neg_overlaps']: #for negative rRNA overlaps  
                            sequence=self.seq('neg',i,seq)                             
                            self.calculations (sequence,rRNAdicts['overlap_counters']) #calculates all the output scores and saves the numbers in the relevent key,value in the counter dictionary
                    rr_over_cols=self.outputs(rRNAdicts['overlap_counters'])
                    
                    
                    #genomewide counts
                    for key in genome.keys():
                        element_dict = locals()[key]['counters']
                        
                        if key=='chromodicts':
                            for small_key in element_dict.keys():
                                genome[key]['counters'][small_key]+=element_dict[small_key]
                            
                        elif key=='cdsdicts':
                            codon_dict=cdsdicts['nucl_codon_counter']
                            for codon_key in codon_dict.keys():
                                genome[key]['nucl_codon_counter'][codon_key]+=codon_dict[codon_key]
                            
                            overlap_element_dict = locals()[key]['overlap_counters']
                            for small_key in element_dict.keys():
                                genome[key]['counters'][small_key]+=element_dict[small_key]
                                genome[key]['overlap_counters'][small_key]+=overlap_element_dict[small_key]
                                
                            
                        else:      
                            overlap_element_dict = locals()[key]['overlap_counters']
                            for small_key in element_dict.keys():
                                genome[key]['counters'][small_key]+=element_dict[small_key]
                                genome[key]['overlap_counters'][small_key]+=overlap_element_dict[small_key]
                                
        
                    #generate chromosomal level outputs
                    if analysis=='c':
                        outputfile.writerow([category,organism,taxa,assembly,col1,version,col2,
                                             chromo_cols[0],chromo_cols[1],chromo_cols[2],chromo_cols[3],chromo_cols[4],chromo_cols[5],chromo_cols[6],chromo_cols[7],chromo_cols[8],chromo_cols[9],
                                             chromo_cols[10],chromo_cols[11],chromo_cols[12],chromo_cols[13],chromo_cols[14],chromo_cols[15],gen_counts[0],gen_counts[1],gen_counts[2],gen_cols[0],gen_cols[1],gen_cols[2],gen_cols[3],gen_cols[4],gen_cols[5],gen_cols[6],gen_cols[7],
                                             gen_cols[8],gen_cols[9],gen_cols[10],gen_cols[11],gen_cols[12],gen_cols[13],gen_cols[14],gen_cols[15],gen_over_cols[5],
                                             cds_counts[0],cds_counts[1],cds_counts[2],cds_cols[0],cds_cols[1],cds_cols[2],cds_cols[3],cds_cols[4],cds_cols[5],cds_cols[6],cds_cols[7],
                                             cds_cols[8],cds_cols[9],cds_cols[10],cds_cols[11],cds_cols[12],cds_cols[13],cds_cols[14],cds_cols[15],cds_over_cols[5],
                                             CDSint_cols[0],CDSint_cols[1],CDSint_cols[2],CDSint_cols[3],CDSint_cols[4],CDSint_cols[5],CDSint_cols[6],CDSint_cols[7],CDSint_cols[8],CDSint_cols[9],
                                             CDSint_cols[10],CDSint_cols[11],CDSint_cols[12],CDSint_cols[13],CDSint_cols[14],CDSint_cols[15],CDSint_over_cols[5],
                                             nc_counts[0],nc_counts[1],nc_counts[2],nc_cols[0],nc_cols[1],nc_cols[2],nc_cols[3],nc_cols[4],nc_cols[5],nc_cols[6],nc_cols[7],nc_cols[8],nc_cols[9],
                                             nc_cols[10],nc_cols[11],nc_cols[12],nc_cols[13],nc_cols[14],nc_cols[15],nc_over_cols[5],
                                             ncint_cols[0],ncint_cols[1],ncint_cols[2],ncint_cols[3],ncint_cols[4],ncint_cols[5],ncint_cols[6],ncint_cols[7],ncint_cols[8],ncint_cols[9],ncint_cols[10],
                                             ncint_cols[11],ncint_cols[12],ncint_cols[13],ncint_cols[14],ncint_cols[15],ncint_over_cols[5],
                                             tr_counts[0],tr_counts[1],tr_counts[2],tr_cols[0],tr_cols[1],tr_cols[2],tr_cols[3],tr_cols[4],tr_cols[5],tr_cols[6],tr_cols[7],tr_cols[8],tr_cols[9],
                                             tr_cols[10],tr_cols[11],tr_cols[12],tr_cols[13],tr_cols[14],tr_cols[15],tr_over_cols[5],
                                             rr_counts[0],rr_counts[1],rr_counts[2],rr_cols[0],rr_cols[1],rr_cols[2],rr_cols[3],rr_cols[4],rr_cols[5],rr_cols[6],rr_cols[7],rr_cols[8],rr_cols[9],
                                             rr_cols[10],rr_cols[11],rr_cols[12],rr_cols[13],rr_cols[14],rr_cols[15],rr_over_cols[5],
                                             codon_cols[0], codon_cols[1], codon_cols[2],codon_cols[3], codon_cols[4], codon_cols[5],codon_cols[6], codon_cols[7],codon_cols[8], codon_cols[9],
                                             codon_cols[10], codon_cols[11], codon_cols[12], codon_cols[13],codon_cols[14], codon_cols[15], codon_cols[16],codon_cols[17], codon_cols[18], 
                                             codon_cols[19], codon_cols[20], codon_cols[21],codon_cols[22], codon_cols[23],codon_cols[24], codon_cols[25], codon_cols[26], codon_cols[27], 
                                             codon_cols[28], codon_cols[29],codon_cols[30],codon_cols[31], codon_cols[32], codon_cols[33], codon_cols[34], codon_cols[35], codon_cols[36], 
                                             codon_cols[37],codon_cols[38], codon_cols[39], codon_cols[40], codon_cols[41], codon_cols[42], codon_cols[43], codon_cols[44],codon_cols[45],
                                             codon_cols[46], codon_cols[47], codon_cols[48], codon_cols[49], codon_cols[50], codon_cols[51],codon_cols[52], codon_cols[53],codon_cols[54], 
                                             codon_cols[55], codon_cols[56], codon_cols[57], codon_cols[58],codon_cols[59], codon_cols[60], codon_cols[61],codon_cols[62], codon_cols[63]]) 
                                            
                        print(f'Completed:{col2}')
                
                locus_counter=locus_counter+1
                
            #generate genome level outputs
            if analysis=='g' and analysed_chromo!=0:
                genome_gen_counts=self.count_cols(genome['genedicts'])
                genome_cds_counts=self.count_cols(genome['cdsdicts'])
                genome_nc_counts=self.count_cols(genome['ncRNAdicts'])
                genome_tr_counts=self.count_cols(genome['tRNAdicts'])
                genome_rr_counts=self.count_cols(genome['rRNAdicts'])
            
                genome_chromo_cols=self.outputs(genome['chromodicts']['counters'])
                genome_gen_cols=self.outputs(genome['genedicts']['counters'])
                genome_cds_cols=self.outputs(genome['cdsdicts']['counters'])
                genome_CDSint_cols=self.outputs(genome['cds_introndicts']['counters'])
                genome_nc_cols=self.outputs(genome['ncRNAdicts']['counters'])
                genome_ncint_cols=self.outputs(genome['nc_introndicts']['counters'])
                genome_tr_cols=self.outputs(genome['tRNAdicts']['counters'])
                genome_rr_cols=self.outputs(genome['rRNAdicts']['counters'])
                
                genome_gen_over_cols=self.outputs(genome['genedicts']['overlap_counters'])
                genome_cds_over_cols=self.outputs(genome['cdsdicts']['overlap_counters'])
                genome_CDSint_over_cols=self.outputs(genome['cds_introndicts']['overlap_counters'])
                genome_nc_over_cols=self.outputs(genome['ncRNAdicts']['overlap_counters'])
                genome_ncint_over_cols=self.outputs(genome['nc_introndicts']['overlap_counters'])
                genome_tr_over_cols=self.outputs(genome['tRNAdicts']['overlap_counters'])
                genome_rr_over_cols=self.outputs(genome['rRNAdicts']['overlap_counters'])
                
                genome_codon_cols=self.codon_outputs(genome['cdsdicts']['nucl_codon_counter'])
            
                outputfile.writerow([category,organism,taxa,assembly,'NA','NA','Whole Genome',
                                     genome_chromo_cols[0],genome_chromo_cols[1],genome_chromo_cols[2],genome_chromo_cols[3],genome_chromo_cols[4],
                                     genome_chromo_cols[5],genome_chromo_cols[6],genome_chromo_cols[7],genome_chromo_cols[8],genome_chromo_cols[9],genome_chromo_cols[10],genome_chromo_cols[11],
                                     genome_chromo_cols[12],genome_chromo_cols[13],genome_chromo_cols[14],genome_chromo_cols[15],genome_gen_counts[0],genome_gen_counts[1],genome_gen_counts[2],
                                     genome_gen_cols[0],genome_gen_cols[1],genome_gen_cols[2],genome_gen_cols[3],genome_gen_cols[4],genome_gen_cols[5],genome_gen_cols[6],genome_gen_cols[7],
                                      genome_gen_cols[8],genome_gen_cols[9],genome_gen_cols[10],genome_gen_cols[11],genome_gen_cols[12],genome_gen_cols[13],genome_gen_cols[14],genome_gen_cols[15],
                                      genome_gen_over_cols[5],genome_cds_counts[0],genome_cds_counts[1],genome_cds_counts[2],genome_cds_cols[0],genome_cds_cols[1],genome_cds_cols[2],genome_cds_cols[3],genome_cds_cols[4],genome_cds_cols[5],
                                      genome_cds_cols[6],genome_cds_cols[7],genome_cds_cols[8],genome_cds_cols[9],genome_cds_cols[10],genome_cds_cols[11],genome_cds_cols[12],genome_cds_cols[13],
                                      genome_cds_cols[14],genome_cds_cols[15],genome_cds_over_cols[5],genome_CDSint_cols[0],genome_CDSint_cols[1],genome_CDSint_cols[2],
                                      genome_CDSint_cols[3],genome_CDSint_cols[4],genome_CDSint_cols[5],genome_CDSint_cols[6],genome_CDSint_cols[7],genome_CDSint_cols[8],genome_CDSint_cols[9],
                                      genome_CDSint_cols[10],genome_CDSint_cols[11],genome_CDSint_cols[12],genome_CDSint_cols[13],genome_CDSint_cols[14],genome_CDSint_cols[15],
                                      genome_CDSint_over_cols[5],genome_nc_counts[0],genome_nc_counts[1],genome_nc_counts[2],genome_nc_cols[0],genome_nc_cols[1],genome_nc_cols[2],genome_nc_cols[3],genome_nc_cols[4],genome_nc_cols[5],genome_nc_cols[6],
                                      genome_nc_cols[7],genome_nc_cols[8],genome_nc_cols[9],genome_nc_cols[10],genome_nc_cols[11],genome_nc_cols[12],genome_nc_cols[13],genome_nc_cols[14],genome_nc_cols[15],
                                      genome_nc_over_cols[5],genome_ncint_cols[0],genome_ncint_cols[1],genome_ncint_cols[2],genome_ncint_cols[3],genome_ncint_cols[4],genome_ncint_cols[5],
                                      genome_ncint_cols[6],genome_ncint_cols[7],genome_ncint_cols[8],genome_ncint_cols[9],genome_ncint_cols[10],genome_ncint_cols[11],genome_ncint_cols[12],
                                      genome_ncint_cols[13],genome_ncint_cols[14],genome_ncint_cols[15],genome_ncint_over_cols[5],genome_tr_counts[0],genome_tr_counts[1],genome_tr_counts[2],genome_tr_cols[0],genome_tr_cols[1],genome_tr_cols[2],
                                      genome_tr_cols[3],genome_tr_cols[4],genome_tr_cols[5],genome_tr_cols[6],genome_tr_cols[7],genome_tr_cols[8],genome_tr_cols[9],genome_tr_cols[10],genome_tr_cols[11],genome_tr_cols[12],
                                      genome_tr_cols[13],genome_tr_cols[14],genome_tr_cols[15],genome_tr_over_cols[5],genome_rr_counts[0],genome_rr_counts[1],genome_rr_counts[2],genome_rr_cols[0],genome_rr_cols[1],genome_rr_cols[2],genome_rr_cols[3],genome_rr_cols[4],
                                      genome_rr_cols[5],genome_rr_cols[6],genome_rr_cols[7],genome_rr_cols[8],genome_rr_cols[9],genome_rr_cols[10],genome_rr_cols[11],genome_rr_cols[12],genome_rr_cols[13],genome_rr_cols[14],
                                      genome_rr_cols[15],genome_rr_over_cols[5],genome_codon_cols[0],genome_codon_cols[1],genome_codon_cols[2], 
                                      genome_codon_cols[3], genome_codon_cols[4],genome_codon_cols[5],genome_codon_cols[6], genome_codon_cols[7], genome_codon_cols[8], genome_codon_cols[9], 
                                      genome_codon_cols[10], genome_codon_cols[11], genome_codon_cols[12], genome_codon_cols[13],genome_codon_cols[14], genome_codon_cols[15], genome_codon_cols[16], 
                                      genome_codon_cols[17], genome_codon_cols[18], genome_codon_cols[19], genome_codon_cols[20], genome_codon_cols[21],genome_codon_cols[22], genome_codon_cols[23], 
                                      genome_codon_cols[24], genome_codon_cols[25], genome_codon_cols[26], genome_codon_cols[27], genome_codon_cols[28], genome_codon_cols[29],genome_codon_cols[30], 
                                      genome_codon_cols[31], genome_codon_cols[32], genome_codon_cols[33], genome_codon_cols[34], genome_codon_cols[35], genome_codon_cols[36], genome_codon_cols[37], 
                                      genome_codon_cols[38], genome_codon_cols[39], genome_codon_cols[40], genome_codon_cols[41], genome_codon_cols[42], genome_codon_cols[43], genome_codon_cols[44], 
                                      genome_codon_cols[45],genome_codon_cols[46], genome_codon_cols[47], genome_codon_cols[48], genome_codon_cols[49], genome_codon_cols[50], genome_codon_cols[51], 
                                      genome_codon_cols[52], genome_codon_cols[53],genome_codon_cols[54], genome_codon_cols[55], genome_codon_cols[56], genome_codon_cols[57], genome_codon_cols[58], 
                                      genome_codon_cols[59], genome_codon_cols[60], genome_codon_cols[61],genome_codon_cols[62], genome_codon_cols[63]])
                
                print("\nAnalysis completed at whole genome level")
                
            if analysed_chromo==0:
                print("\nUnfortunately the input file does not contain any complete chromosomal/mitochondrial or plasmid sequences to analyse.There won't be any data in the output file\nPlease re try with another version or another organism.")
                
            print ("\nSKIPPED",skipped, "sequences due to being either","/".join(skiplocus) )
              

def main():
    try:
        parser = argparse.ArgumentParser(description='''This script parses GenBank flat file format (gbff) and outputs a CSV file.
    The chromosomes with any of the words 'scaffold','unlocalized','contig','unplaced','patch','unknown' will be skipped entirely
    Each column in the output CSV has the following description:
    - Locus ID: Sequence ID
    - Definition: Sequence description
    - n_<element>_pos: Number of <element>s in the positive strand
    - n_<element>_neg: Number of <element>s in the negative strand
    - n_<element>_tot: Total number of <element>s in the chromosome
    - bp_<element>_A: Number of 'A' bases in total <element>
    - bp_<element>_T: Number of 'T' bases in total <element>
    - bp_<element>_C: Number of 'C' bases in total <element>
    - bp_<element>_G: Number of 'G' bases in total <element>
    - bp_<element>_N: Number of 'N' bases in total <element>
    - bp_<element>_total: Total number of bases in total <element>
    - fr_<element>_A: Frequency of 'A' bases among total <element> bases
    - fr_<element>_T: Frequency of 'T' bases among total <element> bases
    - fr_<element>_C: Frequency of 'C' bases among total <element> bases
    - fr_<element>_G: Frequency of 'G' bases among total <element> bases
    - fr_<element>_N: Frequency of 'N' bases among total <element> bases
    - GC_<element>: Percentage of 'G' and 'C' bases in total <element>
    - <element>_topological_entropy: Topological entropy calculated from total <element> sequences
    - <element>_chargaff_pf: Chargaff's second parity rule score (method1) for total <element>
    - <element>_chargaff_ct: Chargaff's second parity rule score (method2) for total <element>
    - <element>_shannon: Shannon entropy score for total <element>
    - bp_<element>_overlap_total: Number of bases overlapping between <element> on two strands
    - ATG, AAG, GTA, etc.: Count of respective codons in CDS.
                ''',formatter_class=argparse.RawTextHelpFormatter)
        
        parser.add_argument("-in", help="gbff input filename.gbff|.txt with path(optional) ", dest="input", type=str, required=True)
        parser.add_argument("-out", help="csv output filename.csv with path(optional)", dest="output", type=str, required=True)
        
        # Add mutually exclusive group for analysis type
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument("-c", help="Run Chromosome level analysis", action="store_true")
        group.add_argument("-g", help="Run Genome level analysis", action="store_true")
        
        
        args = parser.parse_args()

        cool_tool = GBRAP(args.input, args.output)
        
        #Run the relevent analysis based on the user's choice
        
        if args.c:
            cool_tool.run('c')  # Chromosomal level analysis
        elif args.g:
            cool_tool.run('g')  # Genome level analysis
            
    except Exception as e:
        print(f"Error{e}")

if __name__ == "__main__":
    main()
    
print('\nAnalysis Finished\nRun Time :',datetime.now() - startTime)    
    
    

 
