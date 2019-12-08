package main
 
import (
   "os"
   "sync"
   "bytes"
   "strconv"
   "strings"
   "bufio"
   "goModify/altersequence"
   "goModify/concurrentreading"
)

type Command struct {
	// Keeps track of all the user input command information
	InPath string
	OutPath string
	ProportionKeep float64
	LengthSubsequence int
}

func goWriter(writerChannel <-chan [][]byte, done <-chan int, command *Command, numWorkers int, wg *sync.WaitGroup){
	// Receives altered sequences from goWorker through the writerChannel and writes them to output file
	outFile, _ := os.Create(command.OutPath)
	numDoneWorkers := 0
	for {
		select {
		case sequences := <- writerChannel:
			for i := 0; i < len(sequences); i++ {
				outFile.Write(sequences[i])
				outFile.Write([]byte("\n"))
			}
		case doneWorker := <- done:
			numDoneWorkers = numDoneWorkers + doneWorker
			if numDoneWorkers == numWorkers {
				wg.Done()
				return
			}
		}
	}
}

func goWorker(command *Command, readChannel <-chan []byte, writeChannel chan<- [][]byte, doneChannel chan<- int) {
	// Calls on the functions from altersequence to alter the sequences according to user input
	// Reads from readChannel and sends results to writeChannel
	for chunkSequences := range readChannel {
		splitBytes := bytes.SplitAfter(chunkSequences, []byte("\n"))
		for i := 0; i < len(splitBytes) - 1; i += 4 {
			sequenceTruncated := altersequence.RandomSubsequence(splitBytes[i:i+4], command.LengthSubsequence)
			sequenceDiscarded := altersequence.RandomDiscard(sequenceTruncated, command.ProportionKeep)
			writeChannel <- sequenceDiscarded
		}
	}
	doneChannel <- 1
}

func checkChunkSize(fastaFile *os.File) int {
	// Checks the size of each sequence name/sequence/sequence quality pattern
	// This allows for concurrent reading of file chunks (since chunks can then be read independently of one another)
	scanner := bufio.NewScanner(fastaFile)
	chunkSize := 0
	for i := 0; i < 4; i++ {
		if scanner.Scan() {
			chunkSize = chunkSize + len(scanner.Bytes()) + 1
		}
	}
	return chunkSize
}


func runParallel(numWorkerThreads int, numReaderThreads int, command *Command) {
	// Coordinates all the goroutines and channels required for this program to run in parallel
	inFile, _ := os.Open(command.InPath)
	chunkSize := checkChunkSize(inFile)
	reader := concurrentreading.Reader{InFile: inFile, BytesInterval: int64(chunkSize), CurrentOffset: int64(0)}
	writeChannel := make(chan [][]byte)
	doneChannel := make(chan int)
	var wg sync.WaitGroup
	for i := 0; i < numReaderThreads; i++ {
		readChannel := make(chan []byte)
		go reader.ReadConcurrently(readChannel)
		for j := 0; j < numWorkerThreads; j++ {
			go goWorker(command, readChannel, writeChannel, doneChannel)
		}
	}
	wg.Add(1)
	go goWriter(writeChannel, doneChannel, command, numWorkerThreads * numReaderThreads, &wg)
	wg.Wait()
}

func main() {
	// Gets user input and launches program
	proportionKeep, _ := strconv.ParseFloat(strings.Split(os.Args[len(os.Args) - 4], "=")[1], 64)
	lengthSubsequence, _ := strconv.Atoi(strings.Split(os.Args[len(os.Args) - 3], "=")[1])
	inPath := os.Args[len(os.Args) - 2]
	outPath := os.Args[len(os.Args) - 1]
	command := Command{InPath: inPath,
						OutPath: outPath,
						ProportionKeep: proportionKeep,
						LengthSubsequence: lengthSubsequence}
	numWorkerThreads, _ := strconv.Atoi(strings.Split(os.Args[1], "=")[1])
	numReaderWriterThreads, _ := strconv.Atoi(strings.Split(os.Args[2], "=")[1])
	runParallel(numWorkerThreads, numReaderWriterThreads, &command)
}