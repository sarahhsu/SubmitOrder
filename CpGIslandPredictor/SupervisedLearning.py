"""
Author: Caryn Willis
This program takes in a file with sequences on each line
and trains parameters for an HMM for predicting CpG islands
using Maximum Likelyhood Estimation.

"""
import sys
import csv
import math
import numpy as np
import collections
import time
np.set_printoptions(threshold=np.nan)

#learns prior matrix
def LearnPi(sequences,identifiers,hmmprior):
    num_labels=2
    num_nucleotides=4
    pi=np.ones((num_labels,1))
    for i in range (0,len(identifiers)):
        if identifiers[i][0]=='Y':
            pi[0]=pi[0]+1
        else:
            pi[1]=pi[1]+1
    pi_sum=pi[0]+pi[1]
    pi=pi/(pi_sum)
    prior = open(hmmprior, "w")
    for i in range (0,len(pi)):
        j=str(pi[i])
        j=j[1:len(j)-1]
        prior.write("%s\n" %str(j))
    prior.close()
    prior.close()
    return pi

#learns transition matrix
def LearnA(identifiers,hmmtrans):
    A=np.ones((2,2))
    #A is +,- by +,-
    for i in range(0,len(identifiers)):
        for j in range (0,len(identifiers[i])-1):
                if identifiers[i][j]=="+" and identifiers[i][j]=="+":
                    A[0][0]+=1
                elif identifiers[i][j]=="+" and identifiers[i][j]=="-":
                    A[0][1]+=1
                elif identifiers[i][j]=="-" and identifiers[i][j]=="+":
                    A[1][0]+=1
                else:
                    A[1][1]+=1
    for i in range(0,2):
        A_sum=0
        for j in range(0,2):
            A_sum+=A[i][j]
        for j in range(0,2):
            A[i][j]=A[i][j]/A_sum
    trans = open(hmmtrans, "w")
    for i in range (0,len(A)):
        j=str(A[i][0])
        for k in range(1,len(A[i])):
            j=j+" "+str(A[i][k])
        trans.write(j+'\n')
    trans.close()
    return A

#learns emmition matrix
def LearnB(sequences, identifiers, hmmemit):
    nucleotide_dict=dict()
    nucleotide_dict['A']=0
    nucleotide_dict['C']=1
    nucleotide_dict['G']=2
    nucleotide_dict['T']=3
    identifier_dict=dict()
    identifier_dict['+']=0
    identifier_dict['-']=1
    num_nucleotides=4
    num_tags=2
    B=np.ones((num_tags,num_nucleotides))
    #[identifier][nucleotide]
    #A,C,G,T
    for i in range(0,len(identifiers)):
        for j in range (0,len(identifiers[i])):
            B[identifier_dict[identifiers[i][j]]][nucleotide_dict[sequences[i][j]]]+=1
    for i in range(0,num_tags):
        B_sum=0
        for j in range(0,num_nucleotides):
            B_sum+=B[i][j]
        for j in range(0,num_nucleotides):
            B[i][j]=B[i][j]/B_sum
    emit = open(hmmemit, "w")
    for i in range (0,len(B)):
        j=str(B[i][0])
        for k in range(1,len(B[i])):
            j=j+" "+str(B[i][k])
        emit.write(j+'\n')
    emit.close()
    return B

def main():
    train_input=sys.argv[1]
    hmmprior="prior.txt"
    hmmemit="emit.txt"
    hmmtrans="trans.txt"
    ident="trainlabels.txt"
    test="testsequences.txt"
    testid="testlabels.txt"
    in_file = open(train_input, 'rt')
    lines=in_file.readlines()
    sequences=[]
    identifiers=[]
    for i in range(0,35):
        line=lines[i].strip()
        sequences.append(line)
        identifier=50*'-'+(len(line)-100)*'+'+50*'-'
        identifiers.append(identifier)

    identi = open(ident, "w")
    for i in range (0,len(identifiers)):
        identi.write(identifiers[i]+'\n')
    identi.close()

    test15 = open(test, "w")
    for i in range (35,50):
        test15.write(lines[i])

    testid15=open(testid,"w")
    for i in range(35,50):
        line=lines[i].strip()
        identifier=50*'-'+(len(line)-100)*'+'+50*'-'+'\n'
        testid15.write(identifier)
    testid15.close()

    pi=LearnPi(sequences,identifiers,hmmprior)
    A=LearnA(identifiers,hmmtrans)
    B=LearnB(sequences, identifiers, hmmemit)

if __name__ == "__main__":
    main()
