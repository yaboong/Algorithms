/**	Algorithms Analysis Homework3 Knapsack Problem
*
*	Student ID: 20800399
*	Author: 유예본.
*
*	(참고) 
*	채점하시는 분의 편의를 위해 Brute Force의 제한 시간을 바로 아래 BRUTE_FORCE_TIME_LIMIT_SEC 에 지정해 뒀습니다.
*	초 단위로 작동하며, 현재는 주어진 조건에 맞게 20분으로 설정하기 위해 1200으로 되어있습니다.
*
*	N=30 부터는 5분이상이 걸리기 때문에 빠른 확인을 위해서는 BRUTE_FORCE_TIME_LIMIT_SEC 를 1로 설정하시면 됩니다.
*
*   N=x 에서 Time Over 가 발생했다면 x 이후의 input에 대해서는 실행하지 않고 time over를 출력하도록 했습니다.
*/

#define BRUTE_FORCE_TIME_LIMIT_SEC 1200;





#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#pragma warning(disable:4996)

typedef enum { false, true } bool;
int timeOver = false;

//branch and bound 에서 사용 될 Priority Queue와 Node 정의
typedef struct PQ *pq;
typedef struct PQ {
	struct Node* node;
	int capacity;
	int N;
} PQ;

typedef struct Node {
	int benefit;
	int weight;
	float bound;
	int level;
	struct Node* left;
	struct Node* right;
}Node;

//주어진 input 개수에 따라 만들어진 item을 담을 구조체
typedef struct {
	int benefit;
	int weight;
	float benefitPerWeight;
} Item;



//주어진 input 개수에 따라 무작위로 benefit, weight를 정의하고 
//각각의 경우에 따라 DP, Greedy, BruteForce, BranchBound 계산 method호출
void benefit_weight_generator(int N, FILE* fp);


//주어진 4가지 접근방법에 따라 Max Benefit을 계산하는 method
void dynamicProgramming(Item* item, int N, int W, FILE* fp);
void greedy(Item* item, int N, int W, FILE* fp);
void bruteForce(Item* item, int N, int W, FILE* fp);
void branchBound(Item* item, int N, int W, FILE* fp);


/********** Greedy 관련 method **********/
//Greedy에서 benefitPerWeight에 따라 오름차순으로 item을 정렬하기 위해서 heapsort 사용
//(Quicksort, Mergesort 사용시 Stack Overflow 발생하므로 Heapsort사용)
//Descending order heapsort를 구현하기 위한 Min Heap관련 method
void buildMinHeap(Item* arr, int size);
void minHeapify(Item* arr, int size, int k);
void heapSort(Item* arr, int size);
void exch(Item* a, Item* b);


/********** Branch & Bound 관련 method**********/
float boundCalculator(Node* node, Item* item, int N, int MAX, int level, int curBenefit, int curWeight);

Node copyNode(Node* node); //priority queue에 node를 넣을 때 값들만 copy한 새로운 Node를 만들어서 enqueue하기 위한 copy method

//Priority Queue
pq newPQ(int capacity);
int size(pq p);
int resize(pq p, int size);
void enqueue(pq p, Node node);
void swap(pq p, int i, int j);
void swim(pq p, int k);
int less(pq p, int i, int j);
void sink(pq p, int k);
Node dequeue(pq p);
int isEmpty(pq p);
int isFull(pq p);


int main(){
	srand(time(0));

	printf("Start!\n\n\n");
	FILE *fp;
	if (fopen_s(&fp, "output.txt", "wt") != NULL)
	{
		return;
	}
	fprintf(fp, "(Processing time in seconds / Maximum benefit value)\n\n");
	fprintf(fp, "n\tBruteForce\t\t  DP\t\t\t    Greedy\t\t       BranchBound\n\n");
	benefit_weight_generator(10, fp);
	benefit_weight_generator(20, fp);
	benefit_weight_generator(30, fp);
	benefit_weight_generator(40, fp);
	benefit_weight_generator(50, fp);
	benefit_weight_generator(100, fp);
	benefit_weight_generator(500, fp);
	benefit_weight_generator(1000, fp);
	benefit_weight_generator(5000, fp);
	benefit_weight_generator(10000, fp);

	printf("Done! Plese check the output.txt file.\n\n\n");
	return 0;
}

void benefit_weight_generator(int N, FILE* fp){
	Item* item;
	int MAX = N * 25;

	item = (Item*)malloc(sizeof(Item)*(N + 1));

	for (int i = 0; i <= N; i++){
		item[i].benefit = rand() % 500 + 1;
	}

	for (int i = 0; i <= N; i++){
		item[i].weight = rand() % 100 + 1;
	}

	bruteForce(item, N, MAX, fp);

	dynamicProgramming(item, N, MAX, fp);

	greedy(item, N, MAX, fp);

	branchBound(item, N, MAX, fp);

	free(item);
}

void bruteForce(Item* item, int N, int MAX, FILE* fp){ 
	if (timeOver == true){
		printf("N=%d, BruteForce Time Over\n", N);
		fprintf(fp, "%d\t  t_over\t\t  ", N);
		return;
	}

	time_t startTime = 0, endTime = 0, checkPoint = 0, loopStartTime = 0;
	float elapsedTime, loopRunningTime;
	startTime = clock();
	int timeLimit = BRUTE_FORCE_TIME_LIMIT_SEC;

	int j, n = N, W = MAX, curBenefit, curWeight, maxBenefit = 0;
	unsigned long long int i;

	//item배열의 index 1 ~ index N 까지 체크하던 것을
	//안정적인 비트연산을 위해 index 0 ~ index N-1 까지로 복사하여 실행
	Item* copy_item = (Item*)malloc(sizeof(Item)*N);
	for (int i = 0; i < N; i++){
		copy_item[i] = item[i + 1];
	}
	
	//ULL을 붙이지 않으면 32비트 연산까지만 가능하므로 64비트 연산을 위해서 ULL 사용	
	//2^n을 의미하는 upperbound 만큼 loop를 돌게 됨
	unsigned long long upperbound = 1ULL << n;
	
	loopStartTime = clock();
	for (i = 0; i < upperbound; i++){
		checkPoint = clock();
		loopRunningTime = (float)(checkPoint - loopStartTime) / (CLOCKS_PER_SEC);
		if (loopRunningTime > timeLimit){ //loop의 동작시간이 BRUTE_FORCE_TIME_LIMIT_SEC 보다 크면 종료.
			timeOver = true; //이후에 들어오는 N에 대해 brute force를 수행하지 않도록 하기 위함.
			break;
		}
		curBenefit = 0;
		curWeight = 0;
		for (j = 0; j < n; j++){
			if (i & (1 << j)){	// i의 j번째 비트가 1인지 아닌지를 의미한다. 1이면 넣고 아니면 넣지 않는다.
				if (copy_item[j].weight <= W - curWeight){ //이번 item을 넣을 수 있는 경우
					curWeight += copy_item[j].weight;
					curBenefit += copy_item[j].benefit;
				}
			}
		}
		if (curBenefit > maxBenefit)
			maxBenefit = curBenefit;
	}
	free(copy_item);

	endTime = clock();
	elapsedTime = (float)(endTime - startTime) / (CLOCKS_PER_SEC);
	if (timeOver == false){
		printf("N=%d, BruteForce Complete\n", n);
		fprintf(fp, "%-8d%.3f/%-20d", n, elapsedTime, maxBenefit);
	}
	else{
		printf("N=%d, BruteForce Time Over\n", n);
		fprintf(fp, "%d\t  t_over\t\t  ", N);
	}

}

void dynamicProgramming(Item* item, int N, int MAX, FILE* fp){
	time_t startTime = 0, endTime = 0;
	float elapsedTime;
	startTime = clock();

	int n = N, W = MAX, maxBenefit = 0;

	int* pre = (int*)calloc((W + 1), sizeof(int)); //이전 item까지만 고려했을때, 주어진w(1~250)에따라 가질수있는 Benefit
	for (int i = 1; i <= N; i++){
		int* cur = (int*)calloc((W + 1), sizeof(int)); //현재 넣으려는 item에 대한 benefit
		for (int w = 1; w <= W; w++){
			if (item[i].weight <= w){ // 이번에 넣으려는 item이 주어진 w 보다 작거나 같은경우 (넣을수있는경우)
				int leftSpace = w - item[i].weight;
				if (item[i].benefit + pre[leftSpace] > pre[w]) // 이번 item을 넣는게 benefit이 더 큰 경우
					cur[w] = item[i].benefit + pre[leftSpace];
				else
					cur[w] = pre[w];
			}
			else
				cur[w] = pre[w];

			if (cur[w] > maxBenefit)
				maxBenefit = cur[w];
		}
		// 현재까지 계산 한 결과를 pre에 담고 cur는 free해 줌으로써 총 2개의 1D array만 사용하여 계산하기 위함
		for (int j = 0; j <= W; j++) pre[j] = cur[j];
		free(cur);
	}
	free(pre);

	endTime = clock();
	elapsedTime = (float)(endTime - startTime) / (CLOCKS_PER_SEC);
	printf("N=%d, Dynamic Programming Complete\n", n);
	fprintf(fp, "%-3.3f/%-20d", elapsedTime, maxBenefit);
}

void greedy(Item* item, int N, int MAX, FILE* fp){
	time_t startTime = 0, endTime = 0;
	float elapsedTime;

	startTime = clock();

	int n = N;
	int W = MAX;
	float maxBenefit = 0.0;
	int weightSum = 0;
	int i;
	float fraction;

	for (int i = 1; i <= N; i++){
		item[i].benefitPerWeight = (float)item[i].benefit / (float)item[i].weight;
	}

	//NlogN의 time complexity를 가지는 Quicksort, Mergesort, Heapsort 중
	//Quicksort, Mergesort는 재귀호출로 인해 N이 커지면 Stack Overflow가 발생하므로 heapsort사용
	//단위무게당 benefit이 큰 item을 기준으로 내림차순으로 정렬한다
	heapSort(item, N);

	i = 1;
	while (i <= n && weightSum < W){
		if (weightSum + item[i].weight > W){ // item을 쪼개지 않고는 더이상 배낭에 넣을 수 없는 경우
			fraction = (float)(W - weightSum) / (float)item[i].weight;
			weightSum = weightSum + item[i].weight * fraction;
			maxBenefit = (float)maxBenefit + item[i].benefit * fraction;
		}
		else{ // item을 쪼개지 않고 넣을 수 있는 경우
			weightSum = weightSum + item[i].weight;
			maxBenefit = maxBenefit + item[i].benefit;
		}
		i++;
	}

	endTime = clock();
	elapsedTime = (float)(endTime - startTime) / (CLOCKS_PER_SEC);
	printf("N=%d, Greedy Complete\n", n);
	fprintf(fp, "%-3.3f/%-20.3f ", elapsedTime, maxBenefit);
}

void branchBound(Item* item, int n, int W, FILE* fp){
	time_t startTime = 0, endTime = 0;
	float elapsedTime;
	startTime = clock();

	pq p = newPQ(1);
	int N = n;
	int MAX = W;
	int MaxBenefit = 0;
	Node* node = (Node*)malloc(sizeof(Node)*(1));

	// set root
	node->level = 0;
	node->benefit = 0;
	node->weight = 0;
	node->bound = boundCalculator(node, item, N, MAX, node->level, node->benefit, node->weight);
	node->left = NULL;
	node->right = NULL;

	//enqueue하기위해서
	Node copiedNode = copyNode(node);

	//root item을 Queue에 넣고 시작
	enqueue(p, copiedNode);

	//Queue가 빌 때까지 반복
	while (!isEmpty(p)){
		Node testNode = dequeue(p); //Queue에 있는 값들 중 가장 bound값이 큰 값을 dequeue하며 시작

		//Promising하면 if문 실행
		if (testNode.weight < MAX && testNode.bound > testNode.benefit && testNode.bound > MaxBenefit){
			//현재 Priority Queue내에서 가장 큰 bound 값을 가진 노드의 bound가 MaxBenefit 보다 작거나 (지금 Queue에 있는 애들은 뭘 해도 현재 Max보다 커질 수 없다)
			//현재 Priorioty Queue에서 뽑아낸 값의 benefit이 bound보다 크면 종료 (더 진행해도 benefit이 커질 수 없다)
			if (MaxBenefit > (p->node[1]).bound && testNode.benefit > testNode.bound){
				MaxBenefit = testNode.benefit;
				break;
			}


			Node* parent = &testNode;
			Node* newLeft = (Node*)malloc(sizeof(Node));
			Node* newRight = (Node*)malloc(sizeof(Node));

			//set left child
			newLeft->level = parent->level + 1;
			newLeft->benefit = parent->benefit + item[newLeft->level].benefit;
			newLeft->weight = parent->weight + item[newLeft->level].weight;
			newLeft->bound = boundCalculator(node, item, N, MAX, newLeft->level, newLeft->benefit, newLeft->weight);
			if (newLeft->benefit > MaxBenefit && newLeft->weight <= MAX)
				MaxBenefit = newLeft->benefit;
			parent->left = newLeft;
			Node copy = copyNode(newLeft);
			enqueue(p, copy);

			//set right child
			newRight->level = parent->level + 1;
			newRight->benefit = parent->benefit;
			newRight->weight = parent->weight;
			newRight->bound = boundCalculator(node, item, N, MAX, newRight->level, newRight->benefit, newRight->weight);
			if (newRight->benefit > MaxBenefit && newRight->weight <= MAX)
				MaxBenefit = newRight->benefit;
			parent->right = newRight;
			enqueue(p, copyNode(newRight));
		}
	}

	free(node);
	free(p);
	endTime = clock();
	elapsedTime = (float)(endTime - startTime) / (CLOCKS_PER_SEC);

	printf("N=%d, Branch and Bound Complete\n\n", n);
	fprintf(fp, "%.3f/%d\n\n", elapsedTime, MaxBenefit);
}

float boundCalculator(Node* node, Item* item, int N, int MAX, int level, int curBenefit, int curWeight){
	float fraction, bound;
	int flag = false;
	float maxBenefit = 0.0;
	int next = level + 1;
	while (next <= N && curWeight < MAX){
		flag = true;
		if (curWeight + item[next].weight > MAX){
			fraction = (float)(MAX - curWeight) / (float)item[next].weight;
			curWeight = curWeight + item[next].weight * fraction;
			maxBenefit = (float)maxBenefit + (float)item[next].benefit * fraction;
		}
		else{
			curWeight = curWeight + item[next].weight;
			maxBenefit = maxBenefit + item[next].benefit;
		}
		next++;
	}
	flag ? (bound = curBenefit + maxBenefit) : (bound = 0.0);
	return bound;
}

Node copyNode(Node* node){
	Node copyNode;
	copyNode.level = node->level;
	copyNode.benefit = node->benefit;
	copyNode.weight = node->weight;
	copyNode.bound = node->bound;
	copyNode.left = node->left;
	copyNode.right = node->right;
	return copyNode;
}

/***************Heapsort, MinHeap 관련 메소드 ********************/

void buildMinHeap(Item* item, int size){
	for (int i = size / 2; i >= 1; i--)
		minHeapify(item, size, i);
}

void minHeapify(Item* item, int size, int k){
	int min = 0;
	int l = 2 * k;
	int r = 2 * k + 1;

	if (l <= size && item[l].benefitPerWeight < item[k].benefitPerWeight)
		min = l;
	else
		min = k;

	if (r <= size && item[r].benefitPerWeight < item[min].benefitPerWeight)
		min = r;

	if (min != k){
		exch(&item[k], &item[min]);
		minHeapify(item, size, min);
	}
}

void heapSort(Item* item, int size){
	buildMinHeap(item, size);

	for (int i = size; i >= 2; i--){
		exch(&item[1], &item[i]);
		size--;
		minHeapify(item, size, 1);
	}
}

void exch(Item* a, Item* b){
	Item tmp = *a;
	*a = *b;
	*b = tmp;
}



/*************** Priority Queue 관련 메소드 ********************/

pq newPQ(int capacity) {
	pq p = (pq)malloc(sizeof(PQ));

	p->N = 0;
	p->capacity = capacity < 2 ? 2 : capacity;
	p->node = (Node *)malloc(sizeof(Node)* p->capacity);
	return p;
}

int size(pq p) {
	return p->N;
}

int resize(pq p, int newCapacity) {
	p->node = (Node*)realloc(p->node, sizeof(Node)*newCapacity); //newCapacity만큼 새로 node를 동적할당 한다
	return 0;
}

void enqueue(pq p, Node node) {
	if (isFull(p)){
		resize(p, (p->capacity) * 2);//p의 capacity가 가득 찼을 경우 resize()를 불러 공간을 재할당 한다. 
		//resize()함수의 newCapacity로 현재 capacity의 2배만큼 계산하여 넘겨준다
		p->capacity *= 2; //capacity의 숫자를 두배로 늘려 재설정 한다
	}

	p->node[++p->N] = node;
	swim(p, p->N);
}

Node dequeue(pq p) {
	Node max = p->node[1];
	swap(p, 1, p->N--);
	sink(p, 1);

	if ((p->N > 0) && (p->N == (p->capacity - 1) / 4)){ // 실제 들어있는 N의 개수가 capacity의 1/4이 되었을 경우 
		resize(p, (p->capacity) / 2);					// resize함수를 통해 newCapacity로 현재 capacity의 1/2만큼 보내어 재할당 한다
		p->capacity /= 2;// capacity의 숫자를 반으로 줄여 재설정 한다
	}
	return max;
}

void swap(pq p, int i, int j) {
	Node t = p->node[i];
	p->node[i] = p->node[j];
	p->node[j] = t;
}

void swim(pq p, int k) {
	while (k > 1 && less(p, k / 2, k)) {
		swap(p, k / 2, k);
		k = k / 2;
	}
}

int less(pq p, int i, int j) {
	return p->node[i].bound < p->node[j].bound;
}

void sink(pq p, int k) {
	while (2 * k <= p->N) {
		int j = 2 * k;
		if (j < p->N && less(p, j, j + 1)) j++;
		if (!less(p, k, j)) break;
		swap(p, k, j);
		k = j;
	}
}

int isEmpty(pq p) {
	return (p->N == 0) ? true : false;
}

int isFull(pq p) {
	return (p->N == p->capacity - 1) ? true : false;
}

void freePQ(pq p) {
	if (p == NULL) return;
	free(p->node);
	free(p);
}

