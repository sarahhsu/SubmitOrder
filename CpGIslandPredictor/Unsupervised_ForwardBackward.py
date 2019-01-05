"""
Authors: Sarah Hsu and Xuecong Fu
Adapted from Caryn Willis
This program takes a files with the prior, emission, and transition matrices
and make predictions of whether or not it is in a CpG island.
Then it gives the amount correct outside the island, inside the island, and total.
"""
import math
import numpy as np
np.set_printoptions(threshold=np.nan)

def ReadFileMakeLists(file, nucleotideDict):
    in_file = open(file, 'rt')
    observation_list = []
    for line in in_file.readlines():
        sequence = line.strip()
        observations = []
        for i in range (0,len(sequence)):
            observations.append(nucleotideDict[sequence[i]])
        observation_list.append(observations)
    in_file.close()
    return observation_list

def ReadFileMakeMatrix(file):
    in_file = open(file, 'rt')
    matrix = []
    for line in in_file.readlines():
        line = line.strip()
        if line != "":
            array = line.split(" ")
            for i in range (0, len(array)):
                        array[i] = float(array[i].strip())
            matrix.append(array)
    in_file.close()
    return matrix

def ReadFileMakeMatrixPrior(file):
    in_file = open(file, 'rt')
    matrix = []
    for line in in_file.readlines():
        line = line.strip()
        if line != "":
            array = line.split(" ")
            for i in range (0, len(array)):
                        matrix.append(float(array[i].strip()))
    in_file.close()
    return matrix

def CalcAlpha(T, num_tags, B, words_test, pi, s):
    alpha = np.zeros((T,num_tags))
    alpha[0] = np.transpose(B[words_test[s][0]])  *  pi
    for i in range(1,T):
        alpha[i] = np.transpose(B[words_test[s][i]]) * np.dot(np.transpose(A), np.transpose(alpha[i-1]))
    return alpha

def CalcBeta(T, num_tags, A, B, words_test, s):
    beta = np.ones((T, num_tags))
    for i in range(T-2, -1, -1):
        beta[i] = np.dot(A, (B[words_test[s][i+1]] * np.transpose(beta[i+1])))
    return beta

def Predict(alpha, beta, T, num_tags):
    predict = np.zeros((T, num_tags))
    for i in range (0,T):
        predict[i] = alpha[i] * beta[i]
    return predict

def GetMaxIndexes(T, predict):
    predicted = []
    for i in range (0,T):
        max, maxIndex = 0, 0
        for j in range (0, len(predict[i])):
            if predict[i][j] > max:
                max = predict[i][j]
                maxIndex = j
        predicted.append(maxIndex)
    return predicted

def CalcLikelihood(alpha):
    sum = 0.0
    for j in range (0,len(alpha[0])):
        sum += alpha[len(alpha)-1][j]
    if sum == 0:
        return 0
    else:
        return math.log(sum)

def ForwardBackward(words_test, num_words, num_tags, pi, A, B):
    predictLabels=[]
    total = 0.0
    for s in range (0,len(words_test)):
        T = len(words_test[s])
        alpha = CalcAlpha(T, num_tags, B, words_test, pi, s)
        beta = CalcBeta(T, num_tags, A, B, words_test, s)
        predict = Predict(alpha, beta, T, num_tags)
        predicted = GetMaxIndexes(T, predict)
        total += CalcLikelihood(alpha)
        predictLabels.append(predicted)
    print("Likelihood", total)
    print("AverageLikelihood", total/float(len(words_test)))
    return predictLabels


def PrintLabels(words_test, tags, reverseWord, reverseTag, predicted):
    labeled = open(predicted, "w")
    for i in range (0,len(words_test)):
        line = reverseTag[tags[i][0]]
        for j in range (1,len(words_test[i])):
            line += reverseTag[tags[i][j]]
        labeled.write(line + '\n')
    labeled.close()


def InitializeDictionaries():
    nucleotide_dict = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3,
    }
    reverse_nuc_dict = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T',
    }
    identifier_dict = {
    'A+': 0,
    'C+': 1,
    'G+': 2,
    'T+': 3,
    'A-': 4,
    'C-': 5,
    'G-': 6,
    'T-': 7,
    }
    reverse_iden_dict = {
    0: 'A+',
    1: 'C+',
    2: 'G+',
    3: 'T+',
    4: 'A-',
    5: 'C-',
    6: 'G-',
    7: 'T-',
    }
    return nucleotide_dict, reverse_nuc_dict, identifier_dict, reverse_iden_dict


# Gives the amount of predictions that were correct in decimal form
def CalculateCorrect(predicted):
    file = open(predicted, 'rt')
    correct, incorrect = 0, 0
    correct_in, incorrect_in = 0, 0
    correct_out, incorrect_out = 0, 0
    for line in file.readlines():
        line = line.strip()
        for i in range(1, len(line), 2):
            if i < 100 or i > (len(line)-100):
                if line[i] == '-':
                    correct += 1
                    correct_out += 1
                elif line[i] == '+':
                    incorrect += 1
                    incorrect_out += 1
            else:
                if line[i] == '+':
                    correct += 1
                    correct_in += 1
                elif line[i] == '-':
                    incorrect += 1
                    incorrect_in += 1
    return float(correct)/(correct + incorrect), float(correct_in/(correct_in+incorrect_in)), float(correct_out/(correct_out+incorrect_out))


if __name__ == "__main__":
    ##Change these parameters
    num_sequences = 30
    num_iterations = 2000
    test_input = "flydata50.txt"
    hmmprior = "prior.txt"
    hmmemit = "emit.txt"
    hmmtrans = "gradient_descent_matrices/trans_2000_both25.txt"
    predicted = "results_25seq_2000it.txt"

    nucleotide_dict, reverse_nuc_dict, identifier_dict, reverse_iden_dict = InitializeDictionaries()
    num_nucleotides, num_tags = 4, 8
    words_test = ReadFileMakeLists(test_input,nucleotide_dict)
    pi = ReadFileMakeMatrixPrior(hmmprior)
    A = ReadFileMakeMatrix(hmmtrans)
    B = ReadFileMakeMatrix(hmmemit)
    words_test = words_test[35:50]
    tags = ForwardBackward(words_test, num_nucleotides,num_tags,pi, A,B)
    PrintLabels(words_test, tags, reverse_nuc_dict, reverse_iden_dict, predicted)
    accuracy, recall_in, recall_out = CalculateCorrect(predicted)
    print("testing on", test_input, ":", num_sequences, "sequences and", num_iterations, "iterations: overall accuracy", accuracy, "in island", recall_in, "out island", recall_out)
