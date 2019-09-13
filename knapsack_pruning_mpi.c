/* Knapsack calculation based on that of */
/* https://www.tutorialspoint.com/cplusplus-program-to-solve-knapsack-problem-using-dynamic-programming */

#include <stdio.h>
#include <sys/time.h>
#include <stdint.h>
#include <mpi.h>

long int knapSack(long int C, long int w[], long int v[], int n);

uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(int argc, char *argv[]) {
    long int C;    /* capacity of backpack */
    int n;    /* number of items */
    int i;    /* loop counter */

    MPI_Init (&argc, &argv);
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        scanf ("%ld", &C);
        scanf ("%d", &n);
    }

    long int v[n], w[n];        /* value, weight */

    if (rank == 0) {
        for (i = 0; i < n; i++) {
            scanf ("%ld %ld", &v[i], &w[i]);
        }
    }


    uint64_t start = GetTimeStamp ();
    long int ks = knapSack(C, w, v, n);

    if (rank == 0) {
        printf ("knapsack occupancy %ld\n", ks);
        printf ("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
    }

    MPI_Finalize ();

    return 0;
}

/* PLACE YOUR CHANGES BELOW HERE */
#include <stdlib.h>
#include <stddef.h>
#define true 1
#define false 0

typedef struct node{
    long int *arr; //length n always
    struct node * next;
} node;

typedef struct list_node {
    node * head;
    int length;
} list_node;

void insert_into_list(list_node * list,long int *arr, int n){
    node * x = (node*) malloc(sizeof(node));
    x->arr = arr;
    x->next = list->head;
    list->head = x;
    list->length += 1;
}

long int * remove_from_list(list_node *list){
    if (list->length == 0) return NULL;
    node * temp = list->head;
    list->head = temp->next;
    list->length -= 1;
    long int *arr = temp->arr;
    return arr;
}

int empty_list(list_node * list){
    if (list->length == 0) return 1;
    return 0;
}

typedef enum { END_TAG, PBM_TAG, SOLVE_TAG, IDLE_TAG, BnB_TAG, DONE} flag;

typedef struct item{
    long int value;
    long int weight;
} item;

item * inp;
int n;

int nextIdle(int *busy, int n, int* idle){
    int i;
    for (i = 0; i < n; i++){
        if (!busy[i]) {*idle = *idle - 1; busy[i] = true;return i;}
    }
    return -1;
}

void branch(long int *arr, int n, list_node * list){
    int upto;
    int i = 0;
    while (i<n && arr[i]!=-1)
        i++;
    upto = i;
    if (i != n){
        long int *arr1, *arr2;
        arr1 = (long int*)malloc(((n+1)*sizeof(long int)));
        for (i=0;i<n+1;i++)
            arr1[i] = arr[i];
        arr2 = arr;
        arr1[upto] = 0;
        arr2[upto] = 1;
        insert_into_list(list, arr1, n+1);
        insert_into_list(list, arr2, n+1);
    }
}

long int findmax(long int *productlist, int n,long int C){ //n为商品的数量
    long int i=0,total_val=0,weight_used = 0;
    for (; i<n && productlist[i]!=-1; i++){
        if (productlist[i] == 1){
            total_val += inp[i].value;
            weight_used += inp[i].weight;
        }
    }
    if (weight_used>C) return -1;
    long int left = C - weight_used;
    for (; i<n && weight_used<=C; i++){  //多装进去一个
        total_val += inp[i].value;
        weight_used += inp[i].weight;
    }
    return total_val;
}

long int findmin(long int *productlist, int n, long int C){  //得到理想的小的
    long int i=0,total_val=0,weight_used = 0;
    for (; i<n && productlist[i]!=-1; i++){
        if (productlist[i] == 1){
            total_val += inp[i].value;
            weight_used += inp[i].weight;
        }
    }
    for (; i<n && weight_used<=C; i++){
        total_val += inp[i].value;
        weight_used += inp[i].weight;
    }
    if (weight_used > C){
        total_val -= inp[i-1].value;
    }
    return total_val;
}

int compar(const void *a, const void*b){
    item * a1 = (item*)a;
    item * b1 = (item*)b;
    double r1, r2;
    r1 = (double)a1->value/a1->weight;
    r2 = (double)b1->value/b1->weight;
    int res;
    if (r1>r2) res = 1;
    else res = -1;
    return -res;
}


long int knapSack(long int C, long int w[], long int v[], int n)
{
  int myrank, size,i;
  long int res;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  int lengthofitemvalue[] = { 1, 1 };
  MPI_Aint itemoffest[] = { offsetof( item, value),
      offsetof(item, weight)};
  MPI_Datatype array_of_types[] = { MPI_LONG, MPI_LONG };
  MPI_Datatype tmp, MPI_ITEM;
  MPI_Type_create_struct( 2, lengthofitemvalue, itemoffest,array_of_types, &tmp );
  MPI_Aint lb, extent;
  MPI_Type_get_extent( tmp, &lb, &extent );
  MPI_Type_create_resized( tmp, lb, extent, &MPI_ITEM );
  MPI_Type_commit( &MPI_ITEM );

  if (myrank == 0){
      inp = (item*)malloc(n*sizeof(item));
      for (i=0;i<n;i++)
      {
          inp[i].value=v[i];
          inp[i].weight=w[i];
      }
      long int *visited = (long int*) malloc(sizeof(long int)*(n+1));
      long int *tempvisited = (long int*)malloc(sizeof(long int)*(n+1));
      for (i=0;i<n;i++) tempvisited[i] = visited[i] = -1; //初始值为-1
      qsort(inp, n, sizeof(item), compar); //平均价值高的排序
      long int pair[2];
      pair[0] = n; pair[1] = C;
      MPI_Bcast(pair, 2, MPI_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(inp, n, MPI_ITEM, 0, MPI_COMM_WORLD);
      int idle = size - 1;
      flag curflag;
      int *busy = (int*)malloc(size*sizeof(int));
      for (i=1;i<size;i++)
          busy[i] = false;
      busy[0] = true;   //主线程是BUSY
      int dst = nextIdle(busy, size, &idle);
      long int currentweight = 0;
      for (i=0;i<n && currentweight<=C;i++){
          currentweight += inp[i].weight;
          res += inp[i].value;
          visited[i] = 1;
      } //贪婪的加满
      if (currentweight > C){
          visited[i-1] = 0;
          res -= inp[i-1].value;
      } else if(currentweight == C) {
          return res;
      }
      visited[n] = res; //存放到bestSol最后一个值
      tempvisited[n] = res;
      MPI_Send(tempvisited, n+1, MPI_LONG, dst, PBM_TAG, MPI_COMM_WORLD); //发送下一个空闲的线程
      while (idle != size-1){ //直到idle为最后一个线程符号index时
          MPI_Status status;
          MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//接受其他线程发送的信息
          curflag = status.MPI_TAG;
          int comfrom = status.MPI_SOURCE;
          if (curflag == SOLVE_TAG) {
              MPI_Recv(tempvisited, n+1, MPI_LONG, comfrom, curflag, MPI_COMM_WORLD, &status);
              if (tempvisited[n] >= visited[n]){
                  for (i=0;i<n+1;i++)
                      visited[i] = tempvisited[i]; //把temp里面的副给全局里面
                  res = visited[n];
              }
          }
          if (curflag == IDLE_TAG) {
              MPI_Recv(&curflag, 1, MPI_INT, comfrom, curflag, MPI_COMM_WORLD,  &status);
              idle++;
              busy[status.MPI_SOURCE] = 0;
          }
          if (curflag == BnB_TAG) {
              long int data[2];
              MPI_Recv(data, 2, MPI_LONG, comfrom, curflag, MPI_COMM_WORLD, &status);
              long int high = data[0], nSlaves = data[1]; //获取上限 和 list的长度
              if (high > res) {
                  int total= ((nSlaves <= idle)?nSlaves:idle);
                  long int *data = (long int*)malloc((total+1)*sizeof(long int)); //data[total] contains bestSolval
                  data[total] = res;
                  for (i=0;i<total;i++){
                      data[i] = nextIdle(busy, size, &idle);
                  }
                  MPI_Send(data, total+1, MPI_LONG, comfrom, BnB_TAG, MPI_COMM_WORLD);
              }
              else{
                  curflag = DONE;
                  MPI_Send(&curflag, 1, MPI_INT, comfrom, DONE, MPI_COMM_WORLD);
              }
          }
      }
      for (i=1;i<size;i++){ //给所有线程发送终止信息
          curflag = END_TAG;
          MPI_Bsend(&curflag, 1, MPI_INT, i, END_TAG, MPI_COMM_WORLD);
      }
      return res;

  }
  else {   //子线程循环
      long int pair[2];
      MPI_Status status;
      MPI_Bcast(pair, 2, MPI_LONG, 0, MPI_COMM_WORLD);
      n = pair[0]; C = pair[1];
      inp = (item*)malloc(n*sizeof(item));
      MPI_Bcast(inp, n, MPI_ITEM, 0, MPI_COMM_WORLD);
      list_node list;
      list.head = NULL; list.length = 0;
      while (1){
          MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
          flag curflag = status.MPI_TAG;
          int comfrom = status.MPI_SOURCE;
          if (curflag == END_TAG){  //接受终止信息
              return 0;
          }
          if (curflag == PBM_TAG) {  //为主线程第一次发送的数据 为next里的 将其置为busy了
              long int *axSol = (long int*)malloc(sizeof(long int)*(n+1));  //axSol[n] contains bestSol 这个结果来自0号进程的计算
              MPI_Recv(axSol, n+1, MPI_LONG, comfrom, curflag, MPI_COMM_WORLD, &status);//axsol为tempSol
              res = axSol[n]; //得 到了0号进程的结果
              insert_into_list(&list, axSol, n+1); //是个链表来保存中间的结果
              while (!empty_list(&list)){ //当list为不空时
                  long int *auxSp = remove_from_list(&list);
                  long int high = findmax(auxSp, n,C);
                  if (high >= res){
                      long int low = findmin(auxSp, n, C);
                      if (low >= res)  { //结果更好时才返回
                          auxSp[n] = low;
                          res = low;
                          MPI_Send(auxSp, n+1, MPI_LONG, 0, SOLVE_TAG, MPI_COMM_WORLD); //problem
                          //发回0线程
                      }
                      if (low != high){ //problem is here 如果两个不相等
                          long int data[2];
                          data[0] = high;
                          data[1] = list.length;
                          MPI_Send(data, 2, MPI_LONG, 0, BnB_TAG, MPI_COMM_WORLD);////发回0线程
                          long int *assigned = (long int*) malloc((list.length+1)*sizeof(long int));
                          MPI_Recv(assigned, list.length+1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                          if (status.MPI_TAG == BnB_TAG){
                              int total;
                              MPI_Get_count(&status, MPI_LONG, &total);
                              res = assigned[total-1]; total--;
                              branch(auxSp, n, &list);
                              for (i=0;i<total;i++){
                                  long int *sp = remove_from_list(&list);
                                  sp[n] = res;
                                  MPI_Send(sp, n+1, MPI_LONG, assigned[i], PBM_TAG, MPI_COMM_WORLD);
                              }
                          }
                          free(assigned);
                      }
                  }
              }
              curflag = IDLE_TAG;
              MPI_Send(&curflag, 1, MPI_INT, 0, IDLE_TAG, MPI_COMM_WORLD);
          }
      }
  }
  return 0;
}
