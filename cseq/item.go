package cseq

import "container/heap"

// https://pkg.go.dev/container/heap

type item struct {
	value    []int
	priority float64
	index    int
}

type itemPriorityQueue []*item

func (pq itemPriorityQueue) Len() int { return len(pq) }

func (pq itemPriorityQueue) Less(i, j int) bool {
	return pq[i].priority > pq[j].priority
}

func (pq itemPriorityQueue) Swap(i, j int) {
	pq[i], pq[j] = pq[j], pq[i]
	pq[i].index = i
	pq[j].index = j
}

func (pq *itemPriorityQueue) Push(x any) {
	n := len(*pq)
	item := x.(*item)
	item.index = n
	*pq = append(*pq, item)
}

func (pq *itemPriorityQueue) Pop() any {
	old := *pq
	n := len(old)
	item := old[n-1]
	old[n-1] = nil  // avoid memory leak
	item.index = -1 // for safety
	*pq = old[0 : n-1]
	return item
}

func (pq *itemPriorityQueue) update(item *item, value []int, priority float64) {
	item.value = value
	item.priority = priority
	heap.Fix(pq, item.index)
}
