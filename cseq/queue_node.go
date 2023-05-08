package cseq

import (
	"container/heap"
)

type QueueNode struct {
	Ids     []int
	Weight  float64
	AvgDist float64
	index   int
}

type priorityQueue []*QueueNode

func (pq priorityQueue) Len() int { return len(pq) }

func (pq priorityQueue) Less(i, j int) bool {
	if pq[i].Weight == pq[j].Weight {
		return pq[i].AvgDist < pq[j].AvgDist
	}

	return pq[i].Weight > pq[j].Weight
}

func (pq priorityQueue) Swap(i, j int) {
	pq[i], pq[j] = pq[j], pq[i]
	pq[i].index = i
	pq[j].index = j
}

func (pq *priorityQueue) Push(x any) {
	n := len(*pq)
	item := x.(*QueueNode)
	item.index = n
	*pq = append(*pq, item)
}

func (pq *priorityQueue) Pop() any {
	old := *pq
	n := len(old)
	item := old[n-1]
	old[n-1] = nil  // avoid memory leak
	item.index = -1 // for safety
	*pq = old[0 : n-1]
	return item
}

// update modifies the priority and value of an Item in the queue.
func (pq *priorityQueue) update(item *QueueNode, ids []int, weight float64, avgDist float64) {
	item.Ids = ids
	item.Weight = weight
	item.AvgDist = avgDist
	heap.Fix(pq, item.index)
}
