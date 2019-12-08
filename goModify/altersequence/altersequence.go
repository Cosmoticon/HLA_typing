package altersequence

import (
	"sort"
	"math/rand"
)

func RandomSubsequence(sequence [][]byte, lengthSubsequence int) [][]byte {
	// Given a DNA sequence, it outputs the maximum number of non-overlapping subsequences of length "lengthSubsequence" 
	// within that DNA sequence, and modifies the sequence quality information accordingly
	lengthSequence := len(sequence[1])
	numSubsequences := lengthSequence / lengthSubsequence
	outputSubsequences := make([][]byte, numSubsequences * 4)
	indices := make([]int, lengthSequence - (lengthSubsequence-1) * numSubsequences)
	for index := range indices {
		indices[index] = index
	}
	subsequenceNum := 0
	offset := 0
	listAllIndices := rand.Perm(lengthSequence - (lengthSubsequence-1) * numSubsequences)
	listIndices := listAllIndices[0:numSubsequences]
	sort.Ints(listIndices)
	for _, randomIndex := range listIndices {
		outputSubsequences[subsequenceNum*4] = sequence[0]
		outputSubsequences[subsequenceNum*4 + 1] = sequence[1][randomIndex+offset:randomIndex+offset+lengthSubsequence]
		outputSubsequences[subsequenceNum*4 + 2] = sequence[2]
		outputSubsequences[subsequenceNum*4 + 3] = sequence[3][randomIndex+offset:randomIndex+offset+lengthSubsequence]
		offset += lengthSubsequence - 1
		subsequenceNum += 1
	}
	return outputSubsequences
}

func RandomDiscard(sequences [][]byte, proportionKeep float64) [][]byte {
	// Randomly discards with probability "proportionKeep" each sequence from a list of sequences
	var outputSequences [][]byte
	for i := 0; i < len(sequences)-1; i += 4 {
		randomInt := rand.Float64()
		if randomInt < proportionKeep{
			outputSequences = append(outputSequences, sequences[i:i+4]...)
		}
	}
	return outputSequences
}