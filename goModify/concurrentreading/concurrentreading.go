package concurrentreading
 
import (
   "os"
   "sync/atomic"
   "io"
)

type Reader struct {
	// Keeps track of file object address as well as shared pointer within the file
	InFile *os.File
	BytesInterval int64
	CurrentOffset int64
}

func (reader *Reader) ReadConcurrently(readChannel chan<- []byte) {
	// Each reader can read concurrently by upadating the shared reader pointer "CurrentOffset"
	// and reading a chunk of length "BytesInterval", and sends the result to the readChannel
	for {
		currentOffset := atomic.AddInt64(&reader.CurrentOffset, reader.BytesInterval) - reader.BytesInterval
		readBytes := make([]byte, int(reader.BytesInterval))
		_, err := reader.InFile.ReadAt(readBytes, currentOffset)
		readChannel <- readBytes
		if err == io.EOF {
			close(readChannel)
			break
		}
	}
}


