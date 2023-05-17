package main

import (
	"encoding/csv"
	"fmt"
	"github.com/geraldywy/go-cseq/cseq"
	"github.com/geraldywy/go-cseq/model"
	"io"
	"os"
	"strconv"
	"strings"
)

func main() {
	// Open the CSV file for reading
	file, err := os.Open("gaode_poi_2018_value_1000000.csv")
	if err != nil {
		fmt.Println("Error:", err)
		return
	}
	defer file.Close()

	// Create a new CSV reader
	reader := csv.NewReader(file)

	// Set the delimiter to comma
	reader.Comma = ','

	data := make([]*model.DataPoint, 0)
	row, err := reader.Read()
	i := 0
	for row != nil {
		lat, _ := strconv.ParseFloat(row[1], 64)
		lng, _ := strconv.ParseFloat(row[2], 64)
		attVec := make([]bool, 0)
		for i := 4; i < len(row); i++ {
			b, _ := strconv.ParseBool(row[i])
			attVec = append(attVec, b)
		}
		data = append(data, &model.DataPoint{
			Name:            row[0],
			Lat:             lat,
			Lng:             lng,
			Categories:      strings.Split(row[3], ","),
			AttributeVector: attVec,
		})
		i++

		row, err = reader.Read()
	}
	if err != nil && err != io.EOF {
		panic(err)
	}

	cs, err := cseq.Init(data)
	if err != nil {
		panic(err)
	}
	queryRes, err := cs.Query([]int{943103, 837464, 830989}, 5)
	if err != nil {
		panic(err)
	}

	fmt.Printf("query node 1, lat: %f, lng: %f  ", data[943103].Lat, data[943103].Lng)
	fmt.Println()
	fmt.Printf("query node 2, lat: %f, lng: %f  ", data[837464].Lat, data[837464].Lng)
	fmt.Println()
	fmt.Printf("query node 3, lat: %f, lng: %f  ", data[830989].Lat, data[830989].Lng)
	fmt.Println()

	for _, res := range queryRes {
		for _, id := range res.Ids {
			fmt.Printf("id: %d (lat: %f, lng: %f)", id, data[id].Lat, data[id].Lng)
		}
		fmt.Println()
	}
}
