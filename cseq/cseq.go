package cseq

import (
	"container/heap"
	"errors"
	"fmt"
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
		alpha:      0.5,
	}, nil
}

// InitExtras is a overloaded init, accepting extra params to better control cseq execution
func InitExtras(data []*model.DataPoint, rateP, errorBound, alpha float64) (CSEQ, error) {
	if err := verifyInitData(data); err != nil {
		return nil, err
	}

	return &cseq{
		data:       data,
		index:      buildIndex(data),
		rateP:      rateP,
		errorBound: errorBound,
		alpha:      alpha,
	}, nil
}

func buildIndex(data []*model.DataPoint) map[string][]int {
	added := make(map[string]map[int]bool)
	index := make(map[string][]int)
	for i, d := range data {
		for _, cat := range d.Categories {
			if _, exist := added[cat]; !exist {
				added[cat] = make(map[int]bool)
			}
			if _, exist := added[cat][i]; exist {
				continue
			}

			index[cat] = append(index[cat], i)
		}
	}

	return index
}

type CSEQ interface {
	// Query accepts a list of object ids as a query, a desired resultSet size. The result size specified is not
	// guaranteed, and the actual result set returned could be smaller if not enough matches are found.
	Query(query []int, resultSetSize int) ([]*QueueNode, error)
	QueryExplicit(query []int, resultSetSize int, cutoffDist float64, cellSplitParam int) ([]*QueueNode, error)
}

var _ CSEQ = (*cseq)(nil)

type cseq struct {
	data  []*model.DataPoint
	index map[string][]int

	rateP      float64
	errorBound float64
	alpha      float64
}

func (c *cseq) Query(query []int, resultSetSize int) ([]*QueueNode, error) {
	for _, id := range query {
		if id < 0 || id >= len(c.data) { // we differ from the original library, using 0 indexed object ids
			return nil, ErrIllegalIdInQuery
		}
	}

	return c.getTopK(query, resultSetSize, 10000000, 5)
}

func (c *cseq) QueryExplicit(query []int, resultSetSize int, cutoffDist float64, cellSplitParam int) ([]*QueueNode, error) {
	for _, id := range query {
		if id < 0 || id >= len(c.data) { // we differ from the original library, using 0 indexed object ids
			return nil, ErrIllegalIdInQuery
		}
	}

	return c.getTopK(query, resultSetSize, cutoffDist, cellSplitParam)
}

// we stick by param naming convention used in the original CSEQ repository
func (c *cseq) getTopK(query []int, K int, r float64, D int) ([]*QueueNode, error) {
	oriList := make([][]*pNode, 0)
	for _, id := range query {
		pNodeList := c.filterType(c.data[id], r)
		sort.Slice(pNodeList, func(i, j int) bool {
			if pNodeList[i].weight == pNodeList[j].weight {
				if pNodeList[i].distance == pNodeList[j].distance {
					return pNodeList[i].id < pNodeList[j].id
				}
				return pNodeList[i].distance < pNodeList[j].distance
			}
			return pNodeList[i].weight > pNodeList[j].weight
		})
		oriList = append(oriList, pNodeList)
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
	c.splitDFS(oriList, query, lengthLim, D, K, 0, &combination, &que, totMinLat, totMaxLat+1e-8, totMinLng, totMaxLng+1e-8, spatialQueryVector)

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
	que *priorityQueue, totMinLat, totMaxLat, totMinLng, totMaxLng float64, spatialQueryVector []float64) {
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

		c.dfs(gridList, 0, query, combination, que, totMinLat, totMaxLat, totMinLng, totMaxLng, regionMaxVal, K, spatialQueryVector)
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
				c.splitDFS(leftListSet, query, lengthLim, D, K, odd^1, combination, que, totMinLat, midLat, totMinLng, totMaxLng, spatialQueryVector)
			} else {
				c.splitDFS(leftListSet, query, lengthLim, D, K, odd^1, combination, que, totMinLat, totMaxLat, totMinLng, midLng, spatialQueryVector)
			}
		}
		if rightFlag {
			if odd%2 == 0 {
				c.splitDFS(rightListSet, query, lengthLim, D, K, odd^1, combination, que, midLat, totMaxLat, totMinLng, totMaxLng, spatialQueryVector)
			} else {
				c.splitDFS(rightListSet, query, lengthLim, D, K, odd^1, combination, que, totMinLat, totMaxLat, midLng, totMaxLng, spatialQueryVector)
			}
		}
	}
}

func (c *cseq) dfs(gridList []map[int][]*pNode, num int, query []int, combination *[]int, que *priorityQueue,
	totMinLat, totMaxLat, totMinLng, totMaxLng float64, regionMaxVal []float64, K int, spatialQueryVector []float64) {
	if num == len(query) {
		pq := make(itemPriorityQueue, 0)
		var sum float64
		ps := make(map[string]bool) // hash of []int as key
		sumVec := make([]int, 0)

		for i := 0; i < num; i++ {
			sum += gridList[i][(*combination)[i]][0].weight
			sumVec = append(sumVec, 0)
		}

		heap.Push(&pq, &item{
			priority: sum,
			value:    sumVec,
		})

		for i := 0; i < K; i++ {
			if len(pq) == 0 {
				break
			}

			wc := heap.Pop(&pq).(*item)
			candiVector := make([]int, 0)
			inFlag := false

			for j := 0; j < num; j++ {
				candId := gridList[j][(*combination)[j]][wc.value[j]].id
				candiVector = append(candiVector, candId)
				if j == 0 {
					lat1 := c.data[candId].Lat
					lng1 := c.data[candId].Lng

					if lat1 >= totMinLat && lat1 < totMaxLat && lng1 >= totMinLng && lng1 < totMaxLng {
						inFlag = true
					}
				}
			}

			similarity := c.sim(query, candiVector, wc.priority, spatialQueryVector)
			rate := c.distanceRate(query, candiVector)

			if len(*que) >= K && ((1-c.alpha)*wc.priority/float64(len(query))+c.alpha*1.0) < (*que)[0].Weight+c.errorBound {
				return
			}
			if rate > c.rateP || rate < 1.0/c.rateP || !inFlag {
				if rate > (1+c.rateP)*c.rateP {
					return
				}

				i--
			} else {
				if len(*que) < K {
					*que = append(*que, &QueueNode{
						Ids:     candiVector,
						Weight:  similarity,
						AvgDist: c.getAvgDist(candiVector),
					})
				} else if similarity > (*que)[0].Weight {
					heap.Pop(que)
					heap.Push(que, &QueueNode{
						Ids:     candiVector,
						Weight:  similarity,
						AvgDist: c.getAvgDist(candiVector),
					})
				}
			}

			num1 := 0
			if !inFlag {
				num1 = 1
			} else {
				num1 = num
			}

			for j := 0; j < num1; j++ {
				if wc.value[j]+1 < len(gridList[j][(*combination)[j]]) {
					sumTmp := wc.priority - gridList[j][(*combination)[j]][wc.value[j]].weight + gridList[j][(*combination)[j]][wc.value[j]+1].weight
					wc.value[j]++
					if _, exist := ps[fmt.Sprint(wc.value)]; !exist {
						ps[fmt.Sprint(wc.value)] = true
						valCopy := make([]int, len(wc.value))
						copy(valCopy, wc.value)
						heap.Push(&pq, &item{
							value:    valCopy,
							priority: sumTmp,
						})
					}
					wc.value[j]--
				}
			}
		}

		return
	}

	for key := range gridList[num] {
		(*combination)[num] = key
		if len(*que) >= K {
			var UB, cw float64
			for j := 0; j <= num; j++ {
				cw += gridList[j][(*combination)[j]][0].weight
			}
			for j := num + 1; j < len(query); j++ {
				cw += regionMaxVal[j]
			}
			UB = (1-c.alpha)*cw/float64(len(query)) + c.alpha*1.0
			if UB < math.Min(c.errorBound+(*que)[0].Weight, 1.0) {
				continue
			}
		}

		c.dfs(gridList, num+1, query, combination, que, totMinLat, totMaxLat, totMinLng, totMaxLng, regionMaxVal, K, spatialQueryVector)
	}
}

func (c *cseq) distanceRate(query []int, combination []int) float64 {
	v1 := c.getSpatialVector(query)
	v2 := c.getSpatialVector(combination)

	var X, Y float64
	for j := 0; j < len(v1); j++ {
		X += v1[j] * v1[j]
		Y += v2[j] * v2[j]
	}

	X = math.Sqrt(X)
	Y = math.Sqrt(Y)

	return X / Y
}

func (c *cseq) getAvgDist(combination []int) float64 {
	var ans float64
	n := len(combination)
	for i := 0; i < n-1; i++ {
		node1 := c.data[combination[i]]
		for j := i + 1; j < n; j++ {
			node2 := c.data[combination[j]]
			ans += sphericalDist(node1.Lat, node1.Lng, node2.Lat, node2.Lng)
		}
	}

	return 2 * ans / float64(n*(n-1))
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
		node1 := c.data[y[i]]
		for j := i + 1; j < n; j++ {
			node2 := c.data[y[j]]
			res = append(res, sphericalDist(node1.Lat, node1.Lng, node2.Lat, node2.Lng))
		}
	}

	return res
}

func (c *cseq) filterType(egNode *model.DataPoint, r float64) []*pNode {
	res := make([]*pNode, 0)
	used := make([]bool, len(c.data))

	for _, cat := range egNode.Categories {
		for _, candNodeId := range c.index[cat] {
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

func (c *cseq) sim(query []int, combination []int, currentWeight float64, spatialQueryVector []float64) float64 {
	ans := c.alpha * c.getSpatialSim(combination, spatialQueryVector)

	return (1-c.alpha)*currentWeight/float64(len(query)) + ans
}

func (c *cseq) getSpatialSim(combination []int, spatialQueryVector []float64) float64 {
	v2 := c.getSpatialVector(combination)
	var ans float64

	for i := 0; i < len(spatialQueryVector); i++ {
		ans += spatialQueryVector[i] * v2[i]
	}

	var X, Y float64
	for i := 0; i < len(spatialQueryVector); i++ {
		X += spatialQueryVector[i] * spatialQueryVector[i]
		Y += v2[i] * v2[i]
	}

	ans /= math.Sqrt(Y)
	ans = math.Sqrt(X)

	return ans
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
