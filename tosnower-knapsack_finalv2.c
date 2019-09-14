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
    MPI_Bcast(&C, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    long int v[n], w[n];        /* value, weight */
    
    if (rank == 0) {
        for (i = 0; i < n; i++) {
            scanf ("%ld %ld", &v[i], &w[i]);
        }
    }
    MPI_Bcast(v, n, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(w, n, MPI_LONG, 0, MPI_COMM_WORLD);

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
//subblock size
#include<string.h>
#include<stdlib.h>
#include <stddef.h>
#define true 1
#define false 0

//long int knapSack(long int C, long int w[], long int v[], int n) {
//    
//}
typedef enum { END, STARTJOB, FINDBETTER, FREEPRO, SENDJOB, DONE} flag;

typedef struct node{
    long int *curnode;
    struct node * next;
} node;
typedef struct list_node {
    node * head;
    int length;
} list_node;

int hasnode(list_node * list){
    return list->length;
}
void append(list_node * list,long int *curnode, int n){
    node * newenode = (node*) malloc(sizeof(node));
    newenode->curnode = curnode;
    newenode->next = list->head;
    list->head = newenode;
    list->length++;
}
long int * pop(list_node *list){
    if (list->length == 0) return NULL;
    node * temp = list->head;
    list->head = temp->next;
    list->length--;
    return temp->curnode;
}

typedef struct item{
    long int value;
    long int weight;
} item;


int findnextfreeprocess(int *working, int n, int* idle){
    int i;
    for (i = 0; i < n; i++){
        if (!working[i])
        {*idle = *idle - 1;
            working[i] = true;
            return i;
        }
    }
    return -1;
}

void breadennew(list_node * list,long int *curnode, int n ){
    int i = 0;
    while (i<n && curnode[i]!=-1) i++;
    if (i != n){
        long int* newnode = (long int*)malloc(((n+1)*sizeof(long int)));
        for (int j=0;j<n+1;j++)
            newnode[j] = curnode[j];
        curnode[i]=1;
        newnode[i]=0;
        append(list, newnode, n+1);
        append(list, curnode, n+1);
    }
}


long int findmax(long int *productlist,item *items,long int C,int n){ //n为商品的数量
    long int i=0,total_val=0,occupied = 0;
    for (; i<n && productlist[i]!=-1; i++){
        if (productlist[i] == 1){
            total_val += items[i].value;
            occupied += items[i].weight;
        }
    }
    if (occupied>C) return -1;
    long int left = C - occupied;
    for (; i<n && occupied<=C; i++){  //多装进去一个
        total_val += items[i].value;
        occupied += items[i].weight;
    }
    return total_val;
}

long int findmin(long int *productlist,item *items,long int C,int n){  //得到理想的小的
    long int i=0,minvalue=0,occupied = 0;
    for (; i<n && productlist[i]!=-1; i++){
        if (productlist[i] == 1){
            minvalue += items[i].value;
            occupied += items[i].weight;
        }
    }
    for (; i<n && occupied<=C; i++){
        minvalue += items[i].value;
        occupied += items[i].weight;
    }
    if (occupied > C){
        minvalue -= items[i-1].value;
    }
    return minvalue;
}

int comitemaveval(const void *itema, const void*itemb){
    item* item1 = (item*)itema,*item2 = (item*)itemb;
    double avval1=(double)item1->value/item1->weight,avval2=(double)item2->value/item2->weight;
    return avval1>avval2?-1:1;
}

long int pruning(long int C, long int w[], long int v[], int n)
{
    int myrank, size,i,root=0;
    long int res;
    item *items;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int lengthofitemvalue[] = { 1, 1 };
    MPI_Aint itemoffest[] = { offsetof( item, value),
        offsetof(item, weight)};
    MPI_Datatype array_of_types[] = {MPI_LONG, MPI_LONG};
    MPI_Datatype tmp, MPI_ITEM;
    MPI_Type_create_struct( 2, lengthofitemvalue, itemoffest,array_of_types, &tmp );
    MPI_Aint lb, extent;
    MPI_Type_get_extent( tmp, &lb, &extent );
    MPI_Type_create_resized( tmp, lb, extent, &MPI_ITEM );
    MPI_Type_commit( &MPI_ITEM );
    
    if (myrank == 0){
        items = (item*)malloc(n*sizeof(item));
        for (i=0;i<n;i++)
        {
            items[i].value=v[i];
            items[i].weight=w[i];
        }
        long int *visited = (long int*) malloc(sizeof(long int)*(n+1));
        long int *tempvisited = (long int*)malloc(sizeof(long int)*(n+1));
        for (i=0;i<n;i++) tempvisited[i] = visited[i] = -1; //初始值为-1
        qsort(items, n, sizeof(item), comitemaveval); //平均价值高的排序
        long int currentweight = 0;
        for (i=0;i<n && currentweight<=C;i++){
            currentweight += items[i].weight;
            res += items[i].value;
            visited[i] = 1;
        } //贪婪的加满
        if (currentweight > C){
            visited[i-1] = 0;
            res -= items[i-1].value;
        } else if(currentweight == C) {
            return res;
        }
        long int pair[2];
        pair[0] = n; pair[1] = C;
        MPI_Bcast(pair, 2, MPI_LONG, root, MPI_COMM_WORLD);
        MPI_Bcast(items, n, MPI_ITEM, root, MPI_COMM_WORLD);
        int idle = size - 1;
        flag curflag;
        visited[n] = res;
        tempvisited[n] = res;
        int *working = (int*)malloc(size*sizeof(int));
        for (i=1;i<size;i++)
            working[i] = false;
        working[0] = true;
        int dst = findnextfreeprocess(working, size, &idle);
        MPI_Send(tempvisited, n+1, MPI_LONG, dst, STARTJOB, MPI_COMM_WORLD); //发送下一个空闲的线程
        while (idle != size-1){
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//接受其他线程发送的信息
            curflag = status.MPI_TAG;
            int comfrom = status.MPI_SOURCE;
            switch(curflag)
            {
                case FINDBETTER: {
                    MPI_Recv(tempvisited, n + 1, MPI_LONG, comfrom, curflag, MPI_COMM_WORLD, &status);
                    if (tempvisited[n] >= visited[n]) {
                        for (i = 0; i < n + 1; i++)
                            visited[i] = tempvisited[i];
                        res = visited[n];
                    }
                    break;
                }
                case FREEPRO:{
                    MPI_Recv(&curflag, 1, MPI_INT, comfrom, curflag, MPI_COMM_WORLD,  &status);
                    idle++;
                    working[status.MPI_SOURCE] = false;
                    break;
                }
                case SENDJOB:{
                    long int data[2];
                    MPI_Recv(data, 2, MPI_LONG, comfrom, curflag, MPI_COMM_WORLD, &status);
                    long int high = data[0], nSlaves = data[1]; //获取上限 和 list的长度
                    if (high > res) {
                        int total= ((nSlaves <= idle)?nSlaves:idle);
                        long int *data = (long int*)malloc((total+1)*sizeof(long int));
                        data[total] = res;
                        for (i=0;i<total;i++){
                            data[i] = findnextfreeprocess(working, size, &idle);
                        }
                        MPI_Send(data, total+1, MPI_LONG, comfrom, SENDJOB, MPI_COMM_WORLD);
                    }
                    else{
                        MPI_Send(&curflag, 1, MPI_INT, comfrom, DONE, MPI_COMM_WORLD);
                    }
                }
            }
        }
        for (i=1;i<size;i++){ //给所有线程发送终止信息
            MPI_Bsend(&curflag, 1, MPI_INT, i, END, MPI_COMM_WORLD);
        }
        return res;
    }
    else {
        long int pair[2];
        MPI_Bcast(pair, 2, MPI_LONG, 0, MPI_COMM_WORLD);
        n = pair[0]; C = pair[1];
        items = (item*)malloc(n*sizeof(item));
        MPI_Bcast(items, n, MPI_ITEM, 0, MPI_COMM_WORLD);
        list_node list={NULL,0};
        MPI_Status status;
        flag curflag=true;
        int comfrom;
        long int *resfromothers=NULL,*currresult=NULL,high;
        while (curflag){
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            curflag = status.MPI_TAG;
            if (STARTJOB ==curflag) {
                comfrom= status.MPI_SOURCE;
                resfromothers = (long int*)malloc(sizeof(long int)*(n+1));  //axSol[n] contains bestSol 这个结果来自0号进程的计算
                MPI_Recv(resfromothers, n+1, MPI_LONG, comfrom, curflag, MPI_COMM_WORLD, &status);//axsol为tempSol
                res = resfromothers[n]; //得 到了0号进程的结果
                append(&list, resfromothers, n+1); //是个链表来保存中间的结果
                while (hasnode(&list)){ //当list为不空时
                    currresult = pop(&list);
                    high = findmax(currresult,items, C,n);
                    if (high >= res){
                        long int low = findmin(currresult,items, C, n);
                        if (low >= res)  { //结果更好时才返回
                            currresult[n] = low;
                            res = low;
                            MPI_Send(currresult, n+1, MPI_LONG, root, FINDBETTER, MPI_COMM_WORLD); //problem
                            //发回0线程
                        }
                        if (low != high){ //如果两个不相等
                            long int data[2];
                            data[0] = high;
                            data[1] = list.length;
                            MPI_Send(data, 2, MPI_LONG, root, SENDJOB, MPI_COMM_WORLD);//发回0线程
                            long int *givenewjob = (long int*) malloc((list.length+1)*sizeof(long int));
                            MPI_Recv(givenewjob, list.length+1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                            if (SENDJOB == status.MPI_TAG){
                                int len;
                                MPI_Get_count(&status, MPI_LONG, &len);
                                res = givenewjob[(len--)-1];
                                breadennew( &list,currresult,n);
                                for (i=0;i<len;i++){
                                    long int *sp = pop(&list);
                                    sp[n] = res;
                                    MPI_Send(sp, n+1, MPI_LONG, givenewjob[i], STARTJOB, MPI_COMM_WORLD);
                                }
                            }
                            free(givenewjob);
                        }
                    }
                }
                MPI_Send(&curflag, 1, MPI_INT, 0, FREEPRO, MPI_COMM_WORLD);
            }
        }
    }
    return -1;
}


long int dp_mpi(long int C, long int w[], long int v[], int n) {
    ++C;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rowSize = n/size+1;
    int colSize = C/size+1;
    int i = 0;
    if (i==0 && rank==0) {
        int blockRow = (rowSize < n-i)?rowSize:n-i;
        long int blockWeight[blockRow];
        long int blockValue[blockRow];
        for (int j = 0; j < blockRow; j++) {
            blockWeight[j] = w[i + j];
            blockValue[j] = v[i + j];
        }
        long int tmp[blockRow][C];
        for (int j = 0; j < C; j += colSize) {
            int blockCol = (colSize<C-j)? colSize:C-j;
            for (int k = j; k < j + blockCol; k++) {
                if (blockWeight[0] > k) {
                    tmp[0][k] = 0;
                } else {
                    tmp[0][k] = blockValue[0];
                }
            }
            for (int x = 1; x < blockRow; x++) {
                for (int k = j; k < j + blockCol; k++) {
                    if (k < blockWeight[x]) {
                        tmp[x][k] = tmp[x - 1][k];
                    } else if ( tmp[x - 1][k - blockWeight[x]] + blockValue[x] > tmp[x - 1][k]){
                        tmp[x][k] = tmp[x - 1][k - blockWeight[x]] + blockValue[x];
                    } else{
                        tmp[x][k] = tmp[x - 1][k];
                    }
                }
            }
            MPI_Send(&tmp[blockRow - 1][j], blockCol, MPI_LONG, (rank + 1) % size, j, MPI_COMM_WORLD);
        }
        i+=rowSize;
    }
    
    
    for (; i < n; i += rowSize) {
        int blockRow = (rowSize < n-i)?rowSize:n-i;
        long int blockWeight[blockRow];
        long int blockValue[blockRow];
        for (int j = 0; j < blockRow; j++) {
            blockWeight[j] = w[i + j];
            blockValue[j] = v[i + j];
        }
        
        if ((i / rowSize) % size == rank) {

            long int tmp[blockRow + 1][C];
            for (int j = 0; j < C; j += colSize) {
                int blockCol = (colSize<C-j)? colSize:C-j;
//                int blockCol = min(colSize, C - j);
                MPI_Recv(&tmp[0][j], blockCol, MPI_LONG, (rank - 1) % size, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int x = 1; x <= blockRow; x++) {
                    for (int k = j; k < j + blockCol; k++) {
                        if (k < blockWeight[x-1]) {
                            tmp[x][k] = tmp[x - 1][k];
                        } else if ( tmp[x - 1][k - blockWeight[x-1]] + blockValue[x-1] > tmp[x - 1][k]){
                            tmp[x][k] = tmp[x - 1][k - blockWeight[x-1]] + blockValue[x-1];
                        } else{
                            tmp[x][k] = tmp[x - 1][k];
                        }
                    }
                }

                if (i + blockRow == n && j + blockCol == C) {
                    MPI_Send(&tmp[blockRow][C - 1], 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
                } else if (i + blockRow != n) {
                    MPI_Send(&tmp[blockRow][j], blockCol, MPI_LONG, (rank + 1) % size, j, MPI_COMM_WORLD);
                }
            }
        }
    }
    long res = 0;
    if (rank == 0) {
        MPI_Status status;
        int tag;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        tag = status.MPI_TAG;
        int src = status.MPI_SOURCE;
        if (tag == 0) {
            MPI_Recv(&res, n + 1, MPI_LONG, src, tag, MPI_COMM_WORLD, &status);
            return res;
        }
    }
    return 0;
}

long int dp(long int C, long int w[], long int v[], int n) {
    ++C;
    long int tmp[n][C];
    for(int i=0;i<C;++i){
        if (w[0] > i) {
            tmp[0][i] = 0;
        } else {
            tmp[0][i] = v[0];
        }
    }
    for(int i=1;i< n;++i)
    {
        for(int j=0;j<C;++j)
        {
            if (j < w[i]) {
                tmp[i][j] = tmp[i - 1][j];
            } else if ( tmp[i - 1][j - w[i]] + v[i] > tmp[i - 1][j]){
                tmp[i][j] = tmp[i - 1][j - w[i]] + v[i];
            } else{
                tmp[i][j] = tmp[i - 1][j];
            }
        }
    }
    return tmp[n-1][C-1];  //最终求得最大值
}
#include <math.h>
long int bruteForce(long int C, long int w[], long int v[], int n) {
    int i;
    unsigned char used [n], solution [n];
    int done = 0;
    long int wt;
    
    bzero (used, sizeof (used));
    long int max_value = 0;
    long int weight_of_max = 0;
    long int weight = 0;
    long int value = 0;
    while (!done) {
        int carry = 1;
        done = 1;
        for (i = 0; i < n; i++)
            if (!used[i]) {
                used[i] = 1;
                weight += w[i];
                value += v[i];
                done= 0;
                break;
            } else {
                used[i] = 0;
                weight -= w[i];
                value -= v[i];
            }
        if (weight <= C && value > max_value) {
            max_value = value;
            weight_of_max = weight;
            bcopy (used, solution, sizeof (used));
        }
    }
    return max_value;
}

long int knapSack(long int C, long int w[], long int v[], int n) {
//    if (n<31) {
//        if ((1<<n) > n*C) {
//            return dp(C,w,v,n);
//        } else {
//            return pruning(C,w,v,n);
//        }
//    }
//    else {
        return bruteForce(C,w,v,n);
//    }
}
/* mpicc tosnower-knapsack.c -lm -o tosnower-knapsack */
