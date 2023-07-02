#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

/* Structs */
typedef struct Point {
    int x; //x-co-ordinate of the point
    int y; //y-co-ordinate of the point
} Point;

typedef struct Rect {
    Point topRight;      // Top-Right point of MBR
    Point bottomLeft;  // Bottom-Left point of MBR
} Rectangle;

typedef struct Node {
    int isLeaf;         // 1 if the node is a leaf node, 0 otherwise
    int numOfRects;     // Number of rectangles in the node
    Rectangle *rects;        // Rectangles representing MBRs of child nodes
    struct Node **child; // Pointers array to child nodes (if non-leaf node)
    struct Node *parent; //parent of the current node
    Rectangle * extra;
} TreeNode;

typedef struct RTree {
    int m;              // Minimum number of entries in a node
    int M;              // Maximum number of entries in a node
    TreeNode *root;         // Pointer to the root node of the R-tree
} RTree;

typedef struct ListNode {
    Rectangle* data;               // Pointer to data stored in this node
    struct ListNode *next;    // Pointer to the next node in the list
} ListNode;

typedef struct List {
    int size;                 // Number of nodes in the list
    ListNode *head;           // Pointer to the head of the list
    ListNode *tail;           // Pointer to the tail of the list
} List;

//function definitions
TreeNode* createTreeNode(int m, int M, int isLeaf);
Rectangle* createRectFromPoint(Point p);
double calcRectEnlargement(Rectangle *r1, Rectangle *r2);
double calcRectArea(Rectangle r);
TreeNode *chooseLeaf(Rectangle r, TreeNode *n);
Rectangle* getMBR(TreeNode *node);
int getChildIndex(TreeNode *node);
void updateMBR(TreeNode *node);
void splitNode(TreeNode *N, TreeNode *NN);
void adjustTreeSplit(TreeNode *N, TreeNode *NN, RTree *T);
void insert(Point p, RTree *T);
void preOrderTraversal(TreeNode *node);
int overlap(Rectangle r1, Rectangle r2);
void addToList(List *resultList, Rectangle* r);
void searchRectangle(TreeNode *node, Rectangle searchRect, List *resultList);
void freeTree(TreeNode* node);
void freeList(ListNode* head);

//function to create a new Node
TreeNode* createTreeNode(int m, int M, int isLeaf) {
    TreeNode* node = (TreeNode*) malloc(sizeof(TreeNode));
    node->isLeaf = isLeaf; //To decide whether the node is a leaf node or internal node
    node->numOfRects = 0;
    node->extra=NULL;
    node->rects = (Rectangle*) malloc(sizeof(Rectangle) * (M+1));
    node->child = (TreeNode**) malloc(sizeof(TreeNode*) * M);
    node->parent = NULL;
   
    for (int i = 0; i < M; i++) {
        node->child[i] = NULL;
    }
   
    return node;
}


//Considering one point in the data.txt as the center, building a rectangle(square) as the MBR encapsulating the point.
Rectangle* createRectFromPoint(Point p) {
    Rectangle *rect = (Rectangle *)malloc(sizeof(Rectangle));
    if (rect == NULL) {
        printf("Error: Unable to allocate memory for Rectangle.\n");
    }
    rect->topRight.x = p.x + 1;
    rect->topRight.y = p.y + 1;
    rect->bottomLeft.x = p.x - 1;
    rect->bottomLeft.y = p.y - 1;
    
    return rect;
}

//Calculating the enlargement area required for rectangle 1 to accomodate the rectangle 2
double calcRectEnlargement(Rectangle *r1, Rectangle *r2) {
    if(r1==NULL){
        return -1;
    }
    if(r2==NULL){
        return -1;
    }
    double newLeft = fmin(r1->bottomLeft.x, r2->bottomLeft.x);
    double newBottom = fmin(r1->bottomLeft.y, r2->bottomLeft.y);
    double newRight = fmax(r1->topRight.x, r2->topRight.x);
    double newTop = fmax(r1->topRight.y, r2->topRight.y);

    double oldLeft = r1->bottomLeft.x;
    double oldBottom = r1->bottomLeft.y;
    double oldRight = r1->topRight.x;
    double oldTop = r1->topRight.y;

    double newArea = (newRight - newLeft) * (newTop - newBottom);
    double oldArea = (oldRight - oldLeft) * (oldTop - oldBottom);

    return newArea - oldArea;
}

//calculating the area of Rectangle from it's topRight and bottomLeft Points
double calcRectArea(Rectangle r) {
    return (r.topRight.x - r.bottomLeft.x) * (r.topRight.y - r.bottomLeft.y);
}

//chooseleaf function
TreeNode *chooseLeaf(Rectangle r, TreeNode *n) {
    if(n==NULL){
        return NULL;
    }
    if (n->isLeaf==1) {
        return n;
    } else {
        int i, flag = 0;
        double minEnlargement = INFINITY;
        double minArea = INFINITY;
        TreeNode *chosenNode;
        for (i = 0; i < n->numOfRects; i++) {
            double enlargement = calcRectEnlargement(&(n->rects[i]), &r);
            double area = calcRectArea(n->rects[i]);
            if (enlargement < minEnlargement) {
                minEnlargement = enlargement;
                minArea = area;
                chosenNode = n->child[i];
                flag = 0;
            } else if (enlargement == minEnlargement) {
                if (area < minArea) {
                    minArea = area;
                    chosenNode = n->child[i];
                    flag = 0;
                } else if (area == minArea) {
                    chosenNode = n->child[i];
                    flag = 1;
                }
            }
        }
        if (flag == 1) {
            i = n->numOfRects;
            chosenNode = n->child[i];
        }
        return chooseLeaf(r, chosenNode);
    }
}

//getting the Minimum Bounding Rectangle of an Internal Node from it's Child Nodes
Rectangle* getMBR(TreeNode *node) {
    if(node==NULL){
        return NULL;
    }
    Rectangle* mbr = (Rectangle*) malloc(sizeof(Rectangle));
    mbr->bottomLeft.x = node->rects[0].bottomLeft.x;
    mbr->bottomLeft.y = node->rects[0].bottomLeft.y;
    mbr->topRight.x = node->rects[0].topRight.x;
    mbr->topRight.y = node->rects[0].topRight.y;
    for (int i = 1; i < node->numOfRects; i++) {
        mbr->bottomLeft.x = fmin(mbr->bottomLeft.x, node->rects[i].bottomLeft.x);
        mbr->bottomLeft.y = fmin(mbr->bottomLeft.y, node->rects[i].bottomLeft.y);
        mbr->topRight.x = fmax(mbr->topRight.x, node->rects[i].topRight.x);
        mbr->topRight.y = fmax(mbr->topRight.y, node->rects[i].topRight.y);
    }
    return mbr;
}

//getting the index of current Node
int getChildIndex(TreeNode *node) {
    if(node==NULL){
        return -1;
    }
    // Determine the index of the given node in its parent's child array
    int i;
    for (i = 0; i < node->parent->numOfRects; i++) {
        if (node->parent->child[i] == node) {
            return i;
        }
    }
    // If the node is not found in its parent's child array, return -1
    return -1;
}

//to update the MBR of current node to it's parent's child entry
void updateMBR(TreeNode *node) {
    if(node==NULL){
        return;
    }
    if (node->parent == NULL) {
        // Node is the root, so its MBR is already up to date
        return;
    }
    int index = getChildIndex(node);
    Rectangle * mbr = getMBR(node);
    node->parent->rects[index] = *mbr;
    updateMBR(node->parent);
}

//split Node


void splitNode(TreeNode *N, TreeNode *NN) {
    if(N==NULL){
        return;
    }
    if(NN==NULL){
        return;
    }
    // Find the covering rectangle center
    int Xh = INT_MIN, Yh = INT_MIN, XL = INT_MAX, YL = INT_MAX;
    for (int i = 0; i < N->numOfRects; i++) {
        if (N->rects[i].topRight.x > Xh) Xh = N->rects[i].topRight.x;
        if (N->rects[i].topRight.y > Yh) Yh = N->rects[i].topRight.y;
        if (N->rects[i].bottomLeft.x < XL) XL = N->rects[i].bottomLeft.x;
        if (N->rects[i].bottomLeft.y < YL) YL = N->rects[i].bottomLeft.y;
    }
    int CovRectXcen = (Xh + XL) / 2;
    int CovRectYcen = (Yh + YL) / 2;
   
    // Classify objects into corner groups
    int C0 = 0, C1 = 0, C2 = 0, C3 = 0;
    for (int i = 0; i < N->numOfRects+1; i++) {
        int ObjXcen = (N->rects[i].topRight.x + N->rects[i].bottomLeft.x) / 2;
        int ObjYcen = (N->rects[i].topRight.y + N->rects[i].bottomLeft.y) / 2;
        if (ObjXcen > CovRectXcen) {
            if (ObjYcen > CovRectYcen) {
                int x;
                x=C2;
                C2++;
                NN->rects[x] = N->rects[i];
                if(N->isLeaf==0){
                NN->child[x] = N->child[i];
                }
            } else {
               int x;
                x=C3;
                C3++;
                NN->rects[x] = N->rects[i];
                if(N->isLeaf==0){
                NN->child[x] = N->child[i];
                }
            }
        } else {
            if (ObjYcen > CovRectYcen) {
                int x;
                x=C1;
                N->rects[x] = N->rects[i];
                if(N->isLeaf==0){
                N->child[x] = N->child[i];
                }
            } else {
                int x;
                x=C0;
                C0++;
                N->rects[x] = N->rects[i];
                if(N->isLeaf==0){
                N->child[x] = N->child[i];
                }


            }
        }
    }
   
    // Distribute entries according to corner groups
    if (C0 > C2) {
        NN->numOfRects = C2;
        N->numOfRects = C0;
        for (int i = 0; i < C2; i++) {
            N->rects[C0 + i] = NN->rects[i];
                if(N->isLeaf==0){
                N->child[C0 + i] = NN->child[i];
                }
        }
        for (int i = C2; i < C0 + C2; i++) {
            NN->rects[i - C2] = NN->rects[i];
            if(N->isLeaf==0){
                NN->child[i - C2] = NN->child[i];
                }
        }
    } else {
        NN->numOfRects = C0;
        N->numOfRects = C2;
        for (int i = 0; i < C0; i++) {
            NN->rects[C2 + i] = N->rects[i];
            if(N->isLeaf==0){
                NN->child[C2 + i] = N->child[i];
                }
        }
        for (int i = C0; i < C0 + C2; i++) {
            N->rects[i - C0] = N->rects[i];
            if(N->isLeaf==0){
                N->child[i - C0] = N->child[i];
                }
        }
    }
   
    if (C1 > C3) {
        N->numOfRects += C3;
        NN->numOfRects += C1;
        for (int i = 0; i < C3; i++) {
            NN->rects[NN->numOfRects++] = N->rects[N->numOfRects - C3 + i];
            if(N->isLeaf==0){
             NN->child[NN->numOfRects++] = N->child[N->numOfRects - C3 + i];
             }
        }
        for (int i = C3; i < N->numOfRects; i++) {
            NN->rects[NN->numOfRects++] = N->rects[i];
             if(N->isLeaf==0){
             NN->child[NN->numOfRects++] = N->child[i];
             }
        }
    } else {
        N->numOfRects += C1;
        NN->numOfRects += C3;
        for (int i = 0; i < C1; i++) {
            N->rects[N->numOfRects++] = NN->rects[i];
            if(N->isLeaf==0){
              N->child[N->numOfRects++] = NN->child[i];
             }
        }
        for (int i = C1; i < NN->numOfRects; i++) {
            N->rects[N->numOfRects - C1 + i] = NN->rects[i];
            if(N->isLeaf==0){ 
              N->child[N->numOfRects - C1 + i] = NN->child[i];
             }
        }
    }
}



//adjust tree function 
//adjust tree function when nodes were split
void adjustTreeSplit(TreeNode *N, TreeNode *NN, RTree *T) {
    if(T==NULL || N==NULL || NN==NULL){
        return;
    }
    TreeNode* P = N->parent;
    if (P == NULL) {
        // N is the root, so create a new root
        TreeNode* newRoot = createTreeNode(T->m,T->M,0);
        newRoot->child[0] = N;
        N->parent = newRoot;
        newRoot->child[1] = NN;
        NN->parent = newRoot;
        newRoot->numOfRects = 2;
        newRoot->rects[0] = *(getMBR(N));
        newRoot->rects[1] = *(getMBR(NN));
        T->root = newRoot;
        return;
    }
    else if((P->numOfRects)+1 <= T->M){ //there is space for split node in the parent node
        int n = P->numOfRects;
        P->child[n] = NN;
        P->rects[n] = *(getMBR(NN));
        P->numOfRects ++;
        NN->parent = P;
        updateMBR(P);
        return;
    }
    else{//no space in the parent for split node
        TreeNode* PP = createTreeNode(T->m,T->M,0);
        PP->numOfRects = 1;
        PP->child[0] = NN;
        NN->parent = PP;
        PP->rects[0] = *(getMBR(NN)); //created node at parent level with NN as child
        splitNode(P,PP);
        adjustTreeSplit(P,PP,T);
    }
}

//insert function 
void insert(Point p, RTree *T) {
    if(T=NULL){
        return;
    }
    Rectangle* r = createRectFromPoint(p);
    if (T->root == NULL) {
        T->root = createTreeNode(T->m, T->M, 1);
        T->root->rects[0] = *r;
        T->root->numOfRects = 1;
        return;
    }
    TreeNode *n = chooseLeaf(*r, T->root);
    if (n->numOfRects < T->M) {
        // Node has enough space for a new entry
        n->rects[(n->numOfRects)] = *r;
	  (n->numOfRects)++;
        updateMBR(n);
    } else {
        // Split node to obtain L and LL containing E and all the old entries of L
        TreeNode* L = n;
        TreeNode* LL=createTreeNode(T->m,T->M,1);
        splitNode(L, LL);
        // Propagate changes upward
        adjustTreeSplit(L, LL, T);
    }
    return;
}


//preOrder Traversal
void preOrderTraversal(TreeNode *node) {
    if (node == NULL) {
        return;
    }
    if (node->isLeaf==1) {
        // Leaf node
        printf("Leaf Node:\n");
        for (int i = 0; i < node->numOfRects; i++) {
            if(node->child[i]!=NULL){
                int x1 = node->rects[i].bottomLeft.x;
                int y1 = node->rects[i].bottomLeft.y;
                int x2 = node->rects[i].topRight.x;
                int y2 = node->rects[i].topRight.y;
                printf("  Object %d: Encapsulated Point(%d,%d)\n", i + 1, (x1+x2)/2, (y1+y2)/2); //Printing the encapsulated point(center) inside each created rectangle(square)
            }
        }
    } else {
        // Internal node
        printf("Internal Node:\n");
        Rectangle * mbr = getMBR(node); //getting the Minimum Bounding rectangle of the node from the child nodes
        printf("  MBR: Bottom-Left(%d, %d), Top-Right(%d, %d)\n", mbr->bottomLeft.x, mbr->bottomLeft.y, mbr->topRight.x, mbr->topRight.y);
        free(mbr);
        for (int i = 0; i < node->numOfRects; i++) {
            preOrderTraversal(node->child[i]); //Recursively calling on the child nodes
        }
    }
}

//function to find if the two rectangles are overlapping or not
int overlap(Rectangle r1, Rectangle r2) {
    return r1.bottomLeft.x <= r2.topRight.x &&
           r1.topRight.x >= r2.bottomLeft.x &&
           r1.bottomLeft.y <= r2.topRight.y &&
           r1.topRight.y >= r2.bottomLeft.y;
}

//used to add the rectangles which overlap with given rectangle in the linked list
void addToList(List *resultList, Rectangle* r) {
    if(resultList==NULL || r==NULL){
        return;
    }
    ListNode *newNode = (ListNode*) malloc(sizeof(ListNode));
    newNode->data = r;
    newNode->next = NULL;

    if (resultList->head == NULL) {
        resultList->head = newNode;
    } else {
        ListNode *current = resultList->head;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = newNode;
    }
}

//searching a particular rectangle in the tree and corresponding resultList of the overlapping rectangles
void searchRectangle(TreeNode *node, Rectangle searchRect, List *resultList) {
    if(node==NULL || resultList==NULL){
        return;
    }
    if (node->isLeaf==1) {
        for (int i = 0; i < node->numOfRects; i++) {
            if (overlap(node->rects[i], searchRect)) {
                addToList(resultList, &(node->rects[i]));
            }
        }
    } else {
        for (int i = 0; i < node->numOfRects; i++) {
            if (overlap(node->rects[i], searchRect)) {
                searchRectangle(node->child[i], searchRect, resultList);
            }
        }
    }
}

//To free all the memory allocated in the tree
void freeTree(TreeNode* node) {
    if (node == NULL) {
        return;
    }
    if (node->isLeaf==1) {
        free(node->rects);
    } else {
        free(node->rects);
        for (int i = 0; i < node->numOfRects; i++) {
            freeTree(node->child[i]);
        }
        free(node->child);
    }
    free(node);
}
void freeList(ListNode* head) {
    ListNode* temp;
    while (head != NULL) {
        temp = head;
        head = head->next;
        free(temp);
    }
}
//getting input from the file
int main() {
    // Open the data file
    FILE *fp;
    char s[100];

    //getting the file containing the data
    printf("Enter the name of your text file\n");
    scanf("%s", s);
    fp = fopen(s, "r");
    if (fp == NULL) {
        printf("Error: Failed to open data file\n");
        return 1;
    }

    //building the R-Tree for inserting the points
    RTree *rtree = (RTree*) malloc(sizeof(RTree));
    rtree->m = 2;
    rtree->M = 4;
    rtree->root = NULL;

    //Taking points as the input from the file and inserting the rectangle made from it in the tree
    for (int i = 0; !feof(fp); i++) {
        Point* p = (Point *)malloc(sizeof(Point));
        fscanf(fp, "%d %d", &(p->x), &(p->y));
        insert(*p, rtree);
        free(p);
    }

    // Search the tree for a Rectangle
    Rectangle * searchRect = (Rectangle *) malloc(sizeof(Rectangle));
    Point * topRight = (Point *) malloc(sizeof(Point));    
    Point * bottomLeft = (Point *) malloc(sizeof(Point));  
    List * list = (List *) malloc(sizeof(List));

    //taking input for the rectangle to be searched
    printf("Enter the all the integers respectively for the search Rectangle: top-Right x co-ordinate y co-ordinate bottomLeft x co-ordinate y co-ordinate");  
    scanf("%d %d %d %d", &(topRight->x), &(topRight->y), &(bottomLeft->x), &(bottomLeft->y));
    searchRect->bottomLeft = *bottomLeft;
    searchRect->topRight = *topRight;

    //searching for the rectangles overlapped with the given rectangle
    searchRectangle(rtree->root, *searchRect, list);
    ListNode* node = list->head;

    //Printing the details
    printf("Printing the details of overlapped rectangles\n");
    int i=1;
    while(node!=NULL){
        Rectangle * rect = node->data;
        printf("Rect %d: Top-Right(%d,%d) and Bottom-Left(%d,%d)\n", i, rect->topRight.x, rect->topRight.y, rect->bottomLeft.x, rect->bottomLeft.y);
        i++;
        node = node->next;
    }

    //for preOrder Traversal of the Tree
    preOrderTraversal(rtree->root);

    // Free the memory allocated for the R-tree
    freeTree(rtree->root);

    // Free the memory allocated for search rectangle
    free(searchRect);

    // Free the memory allocated for the points
    free(topRight);
    free(bottomLeft);

    // Free the memory of linked-list of overlapped rectangles
    freeList(list->head);
    free(list);
    // Close the file
    fclose(fp);

    return 0;
}