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

	return &cseq{
		data:       data,
		index:      buildIndex(data),
		rateP:      1.5,
		errorBound: 1e-9,
	}, nil
}

// InitExtras is a overloaded init, accepting extra params to better control cseq execution
func InitExtras(data []*model.DataPoint, rateP, errorBound float64) (CSEQ, error) {
	if err := verifyInitData(data); err != nil {
		return nil, err
	}

	return &cseq{
		data:       data,
		index:      buildIndex(data),
		rateP:      rateP,
		errorBound: errorBound,
	}, nil
}

func buildIndex(data []*model.DataPoint) map[string]map[int]bool {
	index := make(map[string]map[int]bool)
	for i, d := range data {
		for _, cat := range d.Categories {
			if _, exist := index[cat]; !exist {
				index[cat] = make(map[int]bool)
			}

			index[cat][i] = true
		}
	}

	return index
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

	rateP      float64
	errorBound float64
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
	newMinLat := 1e18
	newMinLng := 1e18
	newMaxLat := -1e18
	newMaxLng := -1e18

	minLatList := make([]float64, len(oriList))
	minLngList := make([]float64, len(oriList))
	maxLatList := make([]float64, len(oriList))
	maxLngList := make([]float64, len(oriList))
	for i, pList := range oriList {
		minLatList[i] = 1e18
		minLngList[i] = 1e18
		maxLatList[i] = -1e18
		maxLngList[i] = -1e18

		for _, p := range pList {
			newMinLat = math.Min(newMinLat, c.data[p.id].Lat)
			newMinLng = math.Min(newMinLng, c.data[p.id].Lng)
			newMaxLat = math.Max(newMaxLat, c.data[p.id].Lat)
			newMaxLng = math.Max(newMaxLng, c.data[p.id].Lng)
			minLatList[i] = math.Min(minLatList[i], c.data[p.id].Lat)
			minLngList[i] = math.Min(minLngList[i], c.data[p.id].Lng)
			maxLatList[i] = math.Max(maxLatList[i], c.data[p.id].Lat)
			maxLngList[i] = math.Max(maxLngList[i], c.data[p.id].Lng)
		}
	}

	midLat := (totMinLat + totMaxLat) / 2.0
	midLng := (totMinLng + totMaxLng) / 2.0
	distRan1 := sphericalDist(midLat, midLng, midLat, totMaxLng)
	distRan2 := sphericalDist(midLat, midLng, totMaxLat, midLng)

	if (distRan1 <= lengthLim*c.rateP+c.errorBound) || (distRan2 <= lengthLim*c.rateP+c.errorBound) {
		gridList := make([]map[int][]*pNode, len(query))

		totDist := sphericalDist(newMinLat, newMinLng, newMaxLat, newMaxLng)
		rr := math.Min(2.0, 2.0*math.Sqrt(3.0)*c.rateP*totDist/float64(D)/lengthLim+1.0)

		regionMaxVal := make([]float64, 0)
		for i, pList := range oriList {
			gridList[i] = make(map[int][]*pNode)

			var maxAttVal float64
			for _, p := range pList {
				gridId := c.getGridId(p.id, minLatList[i], minLngList[i], maxLatList[i], maxLngList[i], D)
				if float64(len(gridList[i][gridId])) < float64(K)*rr {
					gridList[i][gridId] = append(gridList[i][gridId], p)
				}

				maxAttVal = math.Max(maxAttVal, p.weight)
			}

			regionMaxVal = append(regionMaxVal, maxAttVal)
		}

		dfs(gridList, 0, query, combination, que, totMinLat, totMaxLat, totMinLng, totMaxLng, regionMaxVal)
	} else {
		leftListSet := make([][]*pNode, 0)
		rightListSet := make([][]*pNode, 0)

		leftFlag := true
		rightFlag := true
		for _, pList := range oriList {
			leftList := make([]*pNode, 0)
			rightList := make([]*pNode, 0)

			for _, p := range pList {
				if odd%2 == 0 {
					dist1 := sphericalDist(midLat, c.data[p.id].Lng, c.data[p.id].Lat, c.data[p.id].Lng)

					if (c.data[p.id].Lat < midLat) || (dist1 < lengthLim*c.rateP) {
						leftList = append(leftList, p)
					}

					if (c.data[p.id].Lat >= midLat) || (dist1 < lengthLim*c.rateP) {
						rightList = append(rightList, p)
					}
				} else {
					dist1 := sphericalDist(c.data[p.id].Lat, midLng, c.data[p.id].Lat, c.data[p.id].Lng)

					if (c.data[p.id].Lng < midLng) || (dist1 < lengthLim*c.rateP) {
						leftList = append(leftList, p)
					}

					if (c.data[p.id].Lng >= midLng) || (dist1 < lengthLim*c.rateP) {
						rightList = append(rightList, p)
					}
				}
			}

			if len(leftList) == 0 {
				leftFlag = false
			}
			if len(rightList) == 0 {
				rightFlag = false
			}

			leftListSet = append(leftListSet, leftList)
			rightListSet = append(rightListSet, rightList)
		}

		if leftFlag {
			if odd%2 == 0 {
				c.splitDFS(leftListSet, query, lengthLim, D, K, odd^1, combination, que, totMinLat, midLat, totMinLng, totMaxLng)
			} else {
				c.splitDFS(leftListSet, query, lengthLim, D, K, odd^1, combination, que, totMinLat, totMaxLat, totMinLng, midLng)
			}
		}
		if rightFlag {
			if odd%2 == 0 {
				c.splitDFS(rightListSet, query, lengthLim, D, K, odd^1, combination, que, midLat, totMaxLat, totMinLng, totMaxLng)
			} else {
				c.splitDFS(rightListSet, query, lengthLim, D, K, odd^1, combination, que, totMinLat, totMaxLat, midLng, totMaxLng)
			}
		}
	}
}

func (c *cseq) getGridId(nodeId int, minLat, minLng, maxLat, maxLng float64, dCnt int) int {
	latLength := (maxLat - minLat + 1e-9) / float64(dCnt)
	lngLength := (maxLng - minLng + 1e-9) / float64(dCnt)

	x := int((c.data[nodeId].Lat - minLat) / latLength)
	y := int((c.data[nodeId].Lng - minLng) / lngLength)

	return x*dCnt + y
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
