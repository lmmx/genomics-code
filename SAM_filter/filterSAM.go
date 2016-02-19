package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

/*
func FlagBitToBinary(fbit int) int {
	bin_num := strconv.FormatInt(fbit, 2)
	if len(bin_num) > 8 {
		bin_num = bin_num[len(bin_num)-8:]
	}
	for len(bin_num) < 8 {
		bin_num = "0" + bin_num
	}
	return bin_num
	// process binary as switches with bitwise AND (&) to check presence
}
*/

// I only four valid binary flag bits: 83, 99, 147, 163
// but to parse logically by conversion to binary use the above

var l *log.Logger = log.New(os.Stderr, "", 0)

func ValidFlagBit(fbit int, e error) bool {
	if e != nil {
		return false
		// not a valid FLAG bit, not coercible to integer from string
	}
	switch fbit {
	case
		83,
		99,
		147,
		163:
		return true
	}
	return false
}

func ValidMapQ(q int, e error) bool {
	if e != nil {
		return false
		// not a valid MapQ score, not coercible to integer from string
	}
	// use MAPQ of 5, but see biofinysics.blogspot.co.uk/2014/05/how-does-bowtie2-assign-mapq-scores.html
	min_mapq := 5
	if q < min_mapq {
		return false
	}
	return true
}

func main() {
	file, err := os.Open(os.Args[1])
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	// must use pointers since maps are immutable
	// alternative wastes memory copying each map http://stackoverflow.com/a/24221843/2668831
	paired := make(map[string]bool)    // 1:1 map of string (ID) to boolean (seen or not seen)
	reads := make(map[string][]string) // 1:2 map of string (ID) to 2-tuple of strings (read 1 and read 2)
	// if a third read is seen then something is wrong... highly doubt that happens
	gr, err := gzip.NewReader(file)
	if err != nil {
		log.Fatal(err)
	}
	scan := bufio.NewScanner(gr)
	for scan.Scan() {
		line := scan.Bytes()
		if len(line) == 0 {
			continue
		}
		s_str := string(line)
		s := strings.Split(s_str, "\t")
		if line[0] == '@' {
			fmt.Println(s_str)
		} else {
			flag := s[1]
			mapq := s[4]
			if ValidFlagBit(strconv.Atoi(flag)) && ValidMapQ(strconv.Atoi(mapq)) {
				read_id := s[0]
				// if this is the first read added the value of
				// paired is initialised false for the read's ID key
				reads[read_id] = append(reads[read_id], s_str)
				paired[read_id] = len(reads[read_id]) > 1
				// fmt.Println(read_id, "\tpaired: ", paired[read_id])
			}
		}
	}
	for read_id, is_paired := range paired {
		if is_paired {
			for _, read_str := range reads[read_id] {
				fmt.Println(read_str)
			}
		}
	}
	l.Println("\n\nRead file.")
}
