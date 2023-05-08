package cseq

import (
	"container/heap"
	"errors"
	"github.com/geraldywy/go-cseq/model"
	"math"
	"sort"
)

var (
	ErrNullDataPoint                   = errors.New("null data point")
	ErrInconsistentAttributeVectorSize = errors.New("inconsistent attribute vector size")
	ErrIllegalIdInQuery                = errors.New("illegal id in query")
)

func Init(data []*model.DataPoint) (CSEQ, error) {
	if err := verifyInitData(data); err != nil {
		return nil, err
	}

	index := make(map[string]map[int]bool)
	for i, d := range data {
		for _, cat := range d.Categories {
			if _, exist := index[cat]; !exist {
				index[cat] = make(map[int]bool)
			}

			index[cat][i] = true
		}
	}

	return &cseq{
		data:  data,
		index: index,
	}, nil
}

type CSEQ interface {
	// Query accepts a list of object ids as a query, a desired resultSet size. The result size specified is not
	// guaranteed, and the actual result set returned could be smaller if not enough matches are found.
	Query(query []int, resultSet int, cutoffDist float64, cellSplitParam int) ([]*QueueNode, error)
}

var _ CSEQ = (*cseq)(nil)

type cseq struct {
	data  []*model.DataPoint
	index map[string]map[int]bool
}

func (c *cseq) Query(query []int, resultSet int, cutoffDist float64, cellSplitParam int) ([]*QueueNode, error) {
	for _, id := range query {
		if id < 0 || id >= len(c.data) { // we differ from the original library, using 0 indexed object ids
			return nil, ErrIllegalIdInQuery
		}
	}

	return c.getTopK(query, resultSet, cutoffDist, cellSplitParam)
}

// we stick by param naming convention used in the original CSEQ repository
func (c *cseq) getTopK(query []int, K int, r float64, D int) ([]*QueueNode, error) {
	oriList := make([][]*pNode, 0)
	//oriListSpatial := make([][]*pNode, 0)

	for _, id := range query {
		pNodeList := c.filterType(c.data[id], r)

		pNodeListSpatial := make([]*pNode, len(pNodeList))
		for i, p := range pNodeList {
			pNodeListSpatial[i] = &pNode{
				id:       p.id,
				weight:   p.weight,
				distance: p.distance,
			}
		}

		sort.Slice(pNodeList, func(i, j int) bool {
			return pNodeList[i].weight > pNodeList[j].weight
		})
		//sort.Slice(pNodeListSpatial, func(i, j int) bool {
		//	return pNodeListSpatial[i].distance < pNodeListSpatial[j].distance
		//})
		oriList = append(oriList, pNodeList)
		//oriListSpatial = append(oriListSpatial, pNodeListSpatial)
	}

	minLat := make([]float64, len(oriList))
	minLng := make([]float64, len(oriList))
	maxLat := make([]float64, len(oriList))
	maxLng := make([]float64, len(oriList))

	totMinLat := 1e18
	totMinLng := 1e18
	totMaxLat := -1e18
	totMaxLng := -1e18
	for i, pNodeLis := range oriList {
		minLat[i] = 1e18
		minLng[i] = 1e18
		maxLat[i] = -1e18
		maxLng[i] = -1e18

		for _, node := range pNodeLis {
			minLat[i] = math.Min(minLat[i], c.data[node.id].Lat)
			minLng[i] = math.Min(minLng[i], c.data[node.id].Lng)
			maxLat[i] = math.Max(maxLat[i], c.data[node.id].Lat)
			maxLng[i] = math.Max(maxLng[i], c.data[node.id].Lng)
		}

		totMinLat = math.Min(totMinLat, minLat[i])
		totMinLng = math.Min(totMinLng, minLng[i])
		totMaxLat = math.Max(totMaxLat, maxLat[i])
		totMaxLng = math.Max(totMaxLng, maxLng[i])
	}

	spatialQueryVector := c.getSpatialVector(query)
	var lengthLim float64
	for _, wi := range spatialQueryVector {
		lengthLim += wi * wi
	}
	lengthLim = math.Sqrt(lengthLim)

	combination := make([]int, len(oriList))
	que := make(priorityQueue, 0)
	c.splitDFS(oriList, query, lengthLim, D, K, 0, &combination, &que, totMinLat, totMaxLat+1e-8, totMinLng, totMaxLng+1e-8)

	res := make([]*QueueNode, 0)
	for len(que) > 0 {
		res = append(res, heap.Pop(&que).(*QueueNode))
	}
	for i := 0; i < len(res)/2; i++ {
		res[i], res[len(res)-1-i] = res[len(res)-1-i], res[i]
	}

	return res, nil
}

func (c *cseq) splitDFS(oriList [][]*pNode, query []int, lengthLim float64, D, K int, odd int, combination *[]int,
	que *priorityQueue, totMinLat, totMaxLat, totMinLng, totMaxLng float64) {

}

func (c *cseq) getSpatialVector(y []int) []float64 {
	res := make([]float64, 0)
	n := len(y)
	for i := 0; i < n-1; i++ {
		node1 := c.data[i]
		for j := i + 1; j < n; j++ {
			node2 := c.data[j]
			res = append(res, sphericalDist(node1.Lat, node1.Lng, node2.Lat, node2.Lng))
		}
	}

	return res
}

func (c *cseq) filterType(egNode *model.DataPoint, r float64) []*pNode {
	res := make([]*pNode, 0)
	used := make([]bool, len(c.data))

	for _, cat := range egNode.Categories {
		for candNodeId := range c.index[cat] {
			candNode := c.data[candNodeId]
			if used[candNodeId] || !checkDist(r, egNode, candNode) {
				continue
			}

			res = append(res, &pNode{
				id:       candNodeId,
				weight:   sim(candNode, egNode),
				distance: sphericalDist(egNode.Lat, egNode.Lng, candNode.Lat, candNode.Lng),
			})
			used[candNodeId] = true
		}
	}

	return res
}

func sim(candNode, egNode *model.DataPoint) float64 {
	var fenzi, fenmu1, fenmu2 float64

	for i := range candNode.AttributeVector {
		a := float64(bool2Int(candNode.AttributeVector[i]))
		b := float64(bool2Int(egNode.AttributeVector[i]))

		fenzi += a * b
		fenmu1 += a * a
		fenmu2 += b * b
	}

	if fenmu1 == 0 || fenmu2 == 0 {
		return 0
	}

	return fenzi / (math.Sqrt(fenmu1) * math.Sqrt(fenmu2))
}

func bool2Int(b bool) int {
	if b {
		return 1
	}

	return 0
}

func checkDist(r float64, egNode, candNode *model.DataPoint) bool {
	return sphericalDist(egNode.Lat, egNode.Lng, candNode.Lat, candNode.Lng) <= r
}

func sphericalDist(lat1, lng1, lat2, lng2 float64) float64 {
	radLat1 := rad(lat1)
	radLat2 := rad(lat2)
	a := radLat1 - radLat2
	b := rad(lng1) - rad(lng2)

	s := 2 * math.Asin(math.Sqrt(math.Pow(math.Sin(a/2), 2)+math.Cos(radLat1)*math.Cos(radLat2)*math.Pow(math.Sin(b/2), 2)))
	EarthRadius := 6378.137 // in km

	return s * EarthRadius * 1000
}

func rad(d float64) float64 {
	return d * math.Pi / 180.0
}

func verifyInitData(data []*model.DataPoint) error {
	attVectorSize := -1
	for _, d := range data {
		if d == nil {
			return ErrNullDataPoint
		}
		if attVectorSize == -1 {
			attVectorSize = len(d.AttributeVector)
		} else if attVectorSize != len(d.AttributeVector) {
			return ErrInconsistentAttributeVectorSize
		}
	}

	return nil
}
