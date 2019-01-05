"""
Author: Caryn Willis
This program takes in files with the prior, emission, and transition matrices
and files with test sequences and labels and make predictions of whether or not it is in a CpG island.
Then it gives the amount correct outside the island, inside the island, and total.
"""
import sys
import csv
import math
import numpy as np
import collections
import time
np.set_printoptions(threshold=np.nan)


#Translates the sequences into lists
def ReadFileMakeLists(file, wordDict):
    in_file = open(file, 'rt')
    words_list=[]
    for line in in_file.readlines():
        line=line.strip()
        array=[]
        for i in range (0,len(line)):
            array.append(line[i])

        words=[]
        #labels=[]
        for i in range (0,len(array)):
            words.append(wordDict[array[i]])
            #labels.append(tagDict[array1[1].strip()])
        words_list.append(words)
        #labels_list.append(labels)
    in_file.close()
    return words_list
#takes the matrix files and turns them into matrixes
def ReadFileMakeMatrix(file):
    in_file = open(file, 'rt')
    matrix=[]
    for line in in_file.readlines():
        line.strip()
        if line !="":
            array=line.split(" ")
            array1=[]
            for i in range (0,len(array)):
                        array[i]=array[i].strip()
                        array1.append(float(array[i]))
            matrix.append(array1)
    in_file.close()
    return matrix
#takes the prior file and turns it into a matrix
def ReadFileMakeMatrixPrior(file):
    in_file = open(file, 'rt')
    matrix=[]
    for line in in_file.readlines():
        line.strip()
        if line !="":
            matrix.append(float(line.strip()))

    in_file.close()
    return matrix
#hmmForwardBackword algorithm
def ForwardBackword(words_test, num_words,num_tags,pi, A,B):
    predictLabels=[]
    B=np.transpose(B)
    bigSum=0.0
    for s in range (0,len(words_test)):
        T=len(words_test[s])
        alpha=np.zeros((T,num_tags))
        alpha[0]=np.transpose(B[words_test[s][0]])*pi
        for i in range(1,T):
            alpha[i]=np.transpose(B[words_test[s][i]])*np.dot(np.transpose(A),np.transpose(alpha[i-1]))
        beta=np.ones((T,num_tags))
        i=T-2
        while i >=0:
            beta[i]=np.dot(A,(B[words_test[s][i+1]]*np.transpose(beta[i+1])))
            i-=1
        predict=np.zeros((T,num_tags))
        for i in range (0,T):
            predict[i]=alpha[i]*beta[i]
        predicted=[]
        for i in range (0,T):
            max=0
            maxindex=0
            for j in range (0,len(predict[i])):
                if predict[i][j]>max:
                    max=predict[i][j]
                    maxindex=j
            predicted.append(maxindex)
        likelihood=0

        sum=0.0
        for j in range (0,len(alpha[0])):
            sum=sum+alpha[len(alpha)-1][j]
        if sum==0:
            likelihood=0
        else:
            likelihood=math.log(sum)
        bigSum=bigSum+likelihood
        predictLabels.append(predicted)
    #print("Likelihood",bigSum)
    #print("AverageLikelihood",bigSum/float(len(words_test)))
    return predictLabels

#reverses a dictionary
def ReverseDict(wordDict):
    reverse=dict()
    for key,value in wordDict.items():
        reverse[value]=key
    return reverse

#prints the labels for each sequence
def PrintLabels(words_test, tags, reverseWord, reverseTag, predicted):
    labeled = open(predicted, "w")
    for i in range (0,len(words_test)):
        line=reverseTag[tags[i][0]]
        for j in range (1,len(words_test[i])):
            line=line+reverseTag[tags[i][j]]
        labeled.write(line+'\n')
        labeled.write('\n')
    labeled.close()

def CalculateSuccess(tags,true_labels):
    letter_count=0
    correct_count=0
    incorrect_count=0
    false_pos=0
    false_neg=0
    sequence_rate=[]
    in_island=0
    for i in range (0,len(tags)):
        seq=0
        for j in range(0,len(tags[i])):
            if tags[i][j]==true_labels[i][j]:
                correct_count+=1
                seq+=1
            else:
                incorrect_count+=1
                if tags[i][j]==0:
                    false_pos+=1
                else:
                    false_neg+=1
            if true_labels[i][j]==0:
                in_island+=1
            letter_count+=1
        sequence_rate.append(seq/len(tags[i]))
    false_neg_rate=false_neg/letter_count
    false_pos_rate=false_pos/letter_count
    success_rate=correct_count/letter_count
    island_percent=in_island/letter_count
    out_island=letter_count-in_island
    false_positive=false_pos/out_island
    false_negative=false_neg/in_island
    print("Out of Island Correct:", 1-false_positive)
    print("In Island Correct:", 1-false_negative)
    print("Overall Correct:", success_rate)



def main ():
    test_input=sys.argv[1]
    true_labels_test=sys.argv[2]
    hmmprior="prior.txt"
    hmmemit="emit.txt"
    hmmtrans="trans.txt"
    predicted="train_results.txt"

    nucleotide_dict=dict()
    nucleotide_dict['A']=0
    nucleotide_dict['C']=1
    nucleotide_dict['G']=2
    nucleotide_dict['T']=3
    reverse_nuc_dict=dict()
    reverse_nuc_dict[0]='A'
    reverse_nuc_dict[1]='C'
    reverse_nuc_dict[2]='G'
    reverse_nuc_dict[3]='T'
    reverse_iden_dict=dict()
    reverse_iden_dict[0]='+'
    reverse_iden_dict[1]='-'
    identifier_dict=dict()
    identifier_dict['+']=0
    identifier_dict['-']=1
    num_nucleotides=4
    num_tags=2

    words_test=ReadFileMakeLists(test_input,nucleotide_dict)
    true_labels_test=ReadFileMakeLists(labels,identifier_dict)
    pi=ReadFileMakeMatrixPrior(hmmprior)
    A=ReadFileMakeMatrix(hmmtrans)
    B=ReadFileMakeMatrix(hmmemit)
    tags=ForwardBackword(words_test, num_nucleotides,num_tags,pi, A,B)
    PrintLabels(words_test, tags, reverse_nuc_dict, reverse_iden_dict, predicted)
    CalculateSuccess(tags,true_labels_test)

if __name__ == "__main__":
    main()
