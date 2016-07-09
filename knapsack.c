/**	Algorithms Analysis Homework3 Knapsack Problem
*
*	Student ID: 20800399
*	Author: ������.
*
*	(����) 
*	ä���Ͻô� ���� ���Ǹ� ���� Brute Force�� ���� �ð��� �ٷ� �Ʒ� BRUTE_FORCE_TIME_LIMIT_SEC �� ������ �׽��ϴ�.
*	�� ������ �۵��ϸ�, ����� �־��� ���ǿ� �°� 20������ �����ϱ� ���� 1200���� �Ǿ��ֽ��ϴ�.
*
*	N=30 ���ʹ� 5���̻��� �ɸ��� ������ ���� Ȯ���� ���ؼ��� BRUTE_FORCE_TIME_LIMIT_SEC �� 1�� �����Ͻø� �˴ϴ�.
*
*   N=x ���� Time Over �� �߻��ߴٸ� x ������ input�� ���ؼ��� �������� �ʰ� time over�� ����ϵ��� �߽��ϴ�.
*/

#define BRUTE_FORCE_TIME_LIMIT_SEC 1200;





#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#pragma warning(disable:4996)

typedef enum { false, true } bool;
int timeOver = false;

//branch and bound ���� ��� �� Priority Queue�� Node ����
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

//�־��� input ������ ���� ������� item�� ���� ����ü
typedef struct {
	int benefit;
	int weight;
	float benefitPerWeight;
} Item;



//�־��� input ������ ���� �������� benefit, weight�� �����ϰ� 
//������ ��쿡 ���� DP, Greedy, BruteForce, BranchBound ��� methodȣ��
void benefit_weight_generator(int N, FILE* fp);


//�־��� 4���� ���ٹ���� ���� Max Benefit�� ����ϴ� method
void dynamicProgramming(Item* item, int N, int W, FILE* fp);
void greedy(Item* item, int N, int W, FILE* fp);
void bruteForce(Item* item, int N, int W, FILE* fp);
void branchBound(Item* item, int N, int W, FILE* fp);


/********** Greedy ���� method **********/
//Greedy���� benefitPerWeight�� ���� ������������ item�� �����ϱ� ���ؼ� heapsort ���
//(Quicksort, Mergesort ���� Stack Overflow �߻��ϹǷ� Heapsort���)
//Descending order heapsort�� �����ϱ� ���� Min Heap���� method
void buildMinHeap(Item* arr, int size);
void minHeapify(Item* arr, int size, int k);
void heapSort(Item* arr, int size);
void exch(Item* a, Item* b);


/********** Branch & Bound ���� method**********/
float boundCalculator(Node* node, Item* item, int N, int MAX, int level, int curBenefit, int curWeight);

Node copyNode(Node* node); //priority queue�� node�� ���� �� ���鸸 copy�� ���ο� Node�� ���� enqueue�ϱ� ���� copy method

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

	//item�迭�� index 1 ~ index N ���� üũ�ϴ� ����
	//�������� ��Ʈ������ ���� index 0 ~ index N-1 ������ �����Ͽ� ����
	Item* copy_item = (Item*)malloc(sizeof(Item)*N);
	for (int i = 0; i < N; i++){
		copy_item[i] = item[i + 1];
	}
	
	//ULL�� ������ ������ 32��Ʈ ��������� �����ϹǷ� 64��Ʈ ������ ���ؼ� ULL ���	
	//2^n�� �ǹ��ϴ� upperbound ��ŭ loop�� ���� ��
	unsigned long long upperbound = 1ULL << n;
	
	loopStartTime = clock();
	for (i = 0; i < upperbound; i++){
		checkPoint = clock();
		loopRunningTime = (float)(checkPoint - loopStartTime) / (CLOCKS_PER_SEC);
		if (loopRunningTime > timeLimit){ //loop�� ���۽ð��� BRUTE_FORCE_TIME_LIMIT_SEC ���� ũ�� ����.
			timeOver = true; //���Ŀ� ������ N�� ���� brute force�� �������� �ʵ��� �ϱ� ����.
			break;
		}
		curBenefit = 0;
		curWeight = 0;
		for (j = 0; j < n; j++){
			if (i & (1 << j)){	// i�� j��° ��Ʈ�� 1���� �ƴ����� �ǹ��Ѵ�. 1�̸� �ְ� �ƴϸ� ���� �ʴ´�.
				if (copy_item[j].weight <= W - curWeight){ //�̹� item�� ���� �� �ִ� ���
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

	int* pre = (int*)calloc((W + 1), sizeof(int)); //���� item������ ���������, �־���w(1~250)������ �������ִ� Benefit
	for (int i = 1; i <= N; i++){
		int* cur = (int*)calloc((W + 1), sizeof(int)); //���� �������� item�� ���� benefit
		for (int w = 1; w <= W; w++){
			if (item[i].weight <= w){ // �̹��� �������� item�� �־��� w ���� �۰ų� ������� (�������ִ°��)
				int leftSpace = w - item[i].weight;
				if (item[i].benefit + pre[leftSpace] > pre[w]) // �̹� item�� �ִ°� benefit�� �� ū ���
					cur[w] = item[i].benefit + pre[leftSpace];
				else
					cur[w] = pre[w];
			}
			else
				cur[w] = pre[w];

			if (cur[w] > maxBenefit)
				maxBenefit = cur[w];
		}
		// ������� ��� �� ����� pre�� ��� cur�� free�� �����ν� �� 2���� 1D array�� ����Ͽ� ����ϱ� ����
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

	//NlogN�� time complexity�� ������ Quicksort, Mergesort, Heapsort ��
	//Quicksort, Mergesort�� ���ȣ��� ���� N�� Ŀ���� Stack Overflow�� �߻��ϹǷ� heapsort���
	//�������Դ� benefit�� ū item�� �������� ������������ �����Ѵ�
	heapSort(item, N);

	i = 1;
	while (i <= n && weightSum < W){
		if (weightSum + item[i].weight > W){ // item�� �ɰ��� �ʰ�� ���̻� �賶�� ���� �� ���� ���
			fraction = (float)(W - weightSum) / (float)item[i].weight;
			weightSum = weightSum + item[i].weight * fraction;
			maxBenefit = (float)maxBenefit + item[i].benefit * fraction;
		}
		else{ // item�� �ɰ��� �ʰ� ���� �� �ִ� ���
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

	//enqueue�ϱ����ؼ�
	Node copiedNode = copyNode(node);

	//root item�� Queue�� �ְ� ����
	enqueue(p, copiedNode);

	//Queue�� �� ������ �ݺ�
	while (!isEmpty(p)){
		Node testNode = dequeue(p); //Queue�� �ִ� ���� �� ���� bound���� ū ���� dequeue�ϸ� ����

		//Promising�ϸ� if�� ����
		if (testNode.weight < MAX && testNode.bound > testNode.benefit && testNode.bound > MaxBenefit){
			//���� Priority Queue������ ���� ū bound ���� ���� ����� bound�� MaxBenefit ���� �۰ų� (���� Queue�� �ִ� �ֵ��� �� �ص� ���� Max���� Ŀ�� �� ����)
			//���� Priorioty Queue���� �̾Ƴ� ���� benefit�� bound���� ũ�� ���� (�� �����ص� benefit�� Ŀ�� �� ����)
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

/***************Heapsort, MinHeap ���� �޼ҵ� ********************/

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



/*************** Priority Queue ���� �޼ҵ� ********************/

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
	p->node = (Node*)realloc(p->node, sizeof(Node)*newCapacity); //newCapacity��ŭ ���� node�� �����Ҵ� �Ѵ�
	return 0;
}

void enqueue(pq p, Node node) {
	if (isFull(p)){
		resize(p, (p->capacity) * 2);//p�� capacity�� ���� á�� ��� resize()�� �ҷ� ������ ���Ҵ� �Ѵ�. 
		//resize()�Լ��� newCapacity�� ���� capacity�� 2�踸ŭ ����Ͽ� �Ѱ��ش�
		p->capacity *= 2; //capacity�� ���ڸ� �ι�� �÷� �缳�� �Ѵ�
	}

	p->node[++p->N] = node;
	swim(p, p->N);
}

Node dequeue(pq p) {
	Node max = p->node[1];
	swap(p, 1, p->N--);
	sink(p, 1);

	if ((p->N > 0) && (p->N == (p->capacity - 1) / 4)){ // ���� ����ִ� N�� ������ capacity�� 1/4�� �Ǿ��� ��� 
		resize(p, (p->capacity) / 2);					// resize�Լ��� ���� newCapacity�� ���� capacity�� 1/2��ŭ ������ ���Ҵ� �Ѵ�
		p->capacity /= 2;// capacity�� ���ڸ� ������ �ٿ� �缳�� �Ѵ�
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

