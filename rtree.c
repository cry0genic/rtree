// headers START --------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// headers END ----------------------------------------------------

// define macros START --------------------------------------------
#define MIN 1
#define MAX 4
// define macros END ----------------------------------------------

// define structs START -------------------------------------------

// mbr
typedef struct MBR {
    double MBR[4]; // to store co-ordinates (x_low, x_high, y_low, y_high)
} MBR;

// tuple (to store mbr with id)
typedef struct Tuple {
    int id; // id
    double MBR[4]; // to store co-ordinates (x_low, x_high, y_low, y_high)
} Tuple;

// an entry in r-tree (leaf level)
typedef struct Entry { 
    int id; // id
    double x_low, x_high, y_low, y_high; // co-ordinates of the rectangle
} Entry;

// a node in r-tree
typedef struct Node {
    int id; // id
    Tuple *idx; // ptr to array of tuples
    int idx_count; // idx count
} Node;

// a leaf-node in r-tree
typedef struct LeafNode {
    int id; // id
    Entry *entries; // ptr to array of entries
    int entries_count; // entries count
} LeafNode;

// node union leaf-node
typedef union UNION {
    Node *node; // node
    LeafNode *leaf_node; // leaf-node
} UNION;

// r-tree
typedef struct RTree {
    int id; // id
    Node *root; // root node of the r-tree
    int *child_list; // ptr to array of child indices for upper level creation
    int child_count; // child count
    UNION *node_pool; // ptr to array of union
    int node_pool_count; // node pool counter
    int leaf_counter; // leaf-counter
    double height; // height of the r-tree
    int max_cap; // max node capacity
} RTree;

// define structs END ---------------------------------------------

// define functions START -----------------------------------------

// check memory allocation
void check_memory_allocation(void *ptr) {
    if (ptr == NULL) {
        printf("Memory reallocation failed.\n");
        exit(1);
    };
};

// create tuple
Tuple *create_tuple(int id, double x_low, double x_high, double y_low, double y_high) {
    Tuple *tuple = (Tuple *)malloc(sizeof(Tuple));
    check_memory_allocation(tuple);
    tuple->id = id;
    tuple->MBR[0] = x_low;
    tuple->MBR[1] = x_high;
    tuple->MBR[2] = y_low;
    tuple->MBR[3] = y_high;
    return tuple;
};

// create mbr
MBR *create_mbr(double x_low, double x_high, double y_low, double y_high) {
    MBR *mbr = (MBR *)malloc(sizeof(MBR));
    check_memory_allocation(mbr);
    mbr->MBR[0] = x_low;
    mbr->MBR[1] = x_high;
    mbr->MBR[2] = y_low;
    mbr->MBR[3] = y_high;
    return mbr;
}

// create entry
Entry *create_entry(int id, double x_low, double x_high, double y_low, double y_high) {
    Entry *entry = (Entry *)malloc(sizeof(Entry)); // dynamic memory allocation 
    check_memory_allocation(entry);
    entry->id = id; // assigning values
    entry->x_low = x_low;
    entry->x_high = x_high;
    entry->y_low = y_low;
    entry->y_high = y_high;
    return entry; // return ptr to entry
};

// create node
Node *create_node(int id) {
    Node *node = (Node *)malloc(sizeof(Node)); // dynamic memory allocation
    node->id = id; // assigning values
    node->idx = (Tuple *)malloc(sizeof(Tuple));
    check_memory_allocation(node->idx);
    node->idx_count = 0;
    return node; // return ptr to node
};

// create leaf-node
LeafNode *create_leaf_node(int id, Entry *entries, int entries_count) {
    LeafNode *leaf_node = (LeafNode *)malloc(sizeof(LeafNode)); // dynamic memory allocation
    leaf_node->id = id; // assigning values
    leaf_node->entries = (Entry *)malloc(entries_count * sizeof(Entry)); // allocate memory to ptr
    leaf_node->entries_count = entries_count;
    for (int i = 0; i < entries_count; i++) {
        leaf_node->entries[i] = entries[i];
    };
    return leaf_node; // return ptr to leaf-node
};

// create r-tree
RTree *create_rtree(int max_cap) {
    RTree *rtree = (RTree *)malloc(sizeof(RTree));
    check_memory_allocation(rtree);
    rtree->id = 0; // assigning values
    rtree->root = (Node *)malloc(sizeof(Node));
    check_memory_allocation(rtree->root);
    rtree->child_list = (int *)malloc(sizeof(int));
    check_memory_allocation(rtree->child_list);
    rtree->node_pool = (UNION *)malloc(sizeof(UNION));
    check_memory_allocation(rtree->node_pool);
    rtree->child_count = 0;
    rtree->node_pool_count = 0;
    rtree->leaf_counter = 0;
    rtree->height = 1;
    rtree->max_cap = max_cap;
    return rtree; // return ptr to rtree
};

// calculate mbr for node
MBR *node_mbr(Node *node) {
    double min_xl = node->idx[0].MBR[0]; // finding mbr
    double max_xh = node->idx[0].MBR[1];
    double min_yl = node->idx[0].MBR[2];
    double max_yh = node->idx[0].MBR[3];
    for (int i = 0; i < node->idx_count; i++) {
        min_xl = fmin(min_xl, node->idx[i].MBR[0]);
        max_xh = fmax(max_xh, node->idx[i].MBR[1]);
        min_yl = fmin(min_yl, node->idx[i].MBR[2]);
        max_yh = fmax(max_yh, node->idx[i].MBR[3]);
    };
    MBR *mbr = create_mbr(min_xl, max_xh, min_yl, max_yh);
    return mbr;
};

// calculate mbr for leaf-node
MBR *leaf_node_mbr(LeafNode *leaf_node) {
    double min_xl = leaf_node->entries[0].x_low; // finding mbr
    double max_xh = leaf_node->entries[0].x_high;
    double min_yl = leaf_node->entries[0].y_low;
    double max_yh = leaf_node->entries[0].y_high;
    for (int i = 0; i < leaf_node->entries_count; i++) {
        min_xl = fmin(min_xl, leaf_node->entries[i].x_low);
        max_xh = fmax(max_xh, leaf_node->entries[i].x_high);
        min_yl = fmin(min_yl, leaf_node->entries[i].y_low);
        max_yh = fmax(max_yh, leaf_node->entries[i].y_high);
    };
    MBR *mbr = create_mbr(min_xl, max_xh, min_yl, max_yh);
    return mbr;
};

// calculate mbr area (can be used for both node & leaf-node)
double mbr_area(Tuple *tuple) {
    return ((tuple->MBR[1] - tuple->MBR[0]) * (tuple->MBR[3] * tuple->MBR[2]));
};

// insert leaf-node into r-tree
void insert_leaf(RTree *rtree, Entry *entries, int entries_count) {
    rtree->leaf_counter++;
    if (rtree->child_count > 0) {
        rtree->child_list = (int *)realloc(rtree->child_list, (rtree->child_count + 1) * sizeof(int)); // reallocate memory
        check_memory_allocation(rtree->child_list);
    };
    rtree->child_list[rtree->child_count++] = rtree->id++;
    if (rtree->node_pool_count > 0) {
        rtree->node_pool = (UNION *)realloc(rtree->node_pool, (rtree->node_pool_count + 1) * sizeof(UNION)); // reallocate memory
        check_memory_allocation(rtree->node_pool);
    };
    rtree->node_pool[rtree->node_pool_count++].leaf_node = create_leaf_node(rtree->id++, entries, entries_count);
    fflush(stdout);
};

// creates and inserts node objects
void insert_nodes(RTree *rtree) {
    int number_of_slices = (rtree->child_count + rtree->max_cap - 1) / rtree->max_cap;
    int buffer[number_of_slices];
    int child_indexes[number_of_slices][rtree->max_cap];
    int idx = 0;
    for (int i = 0; i < number_of_slices; i++) {
        for (int j = 0; j < rtree->max_cap && idx < rtree->child_count; j++, idx++) {
            child_indexes[i][j] = rtree->child_list[idx];
        };
    };
    for (int i = 0; i < number_of_slices; i++) {
        buffer[i] = rtree->id;
        Node *new_node = create_node(rtree->id);
        for(int j = 0; j < rtree->max_cap && j < rtree->child_count - i * rtree->max_cap; j++) {
            int child_idx = child_indexes[i][j];
            if (new_node->idx_count > 0) {
                new_node->idx = (Tuple *)realloc(rtree->root->idx, (rtree->root->idx_count + 1) * sizeof(Tuple));
                check_memory_allocation(rtree->root->idx);
            };
            if (rtree->node_pool[child_idx].node != NULL) {
                Node *child_node = rtree->node_pool[child_idx].node;
                MBR *mbr = node_mbr(child_node);
                new_node->idx[new_node->idx_count++] = *create_tuple(child_node->id, mbr->MBR[0], mbr->MBR[1], mbr->MBR[2], mbr->MBR[3]);
            };
            if (rtree->node_pool[child_idx].leaf_node != NULL) {
                LeafNode *child_leaf_node = rtree->node_pool[child_idx].leaf_node;
                MBR *mbr = leaf_node_mbr(child_leaf_node);
                new_node->idx[new_node->idx_count++] = *create_tuple(child_leaf_node->id, mbr->MBR[0], mbr->MBR[1], mbr->MBR[2], mbr->MBR[3]);
            };
        };  
        if (rtree->node_pool_count > 0) {
            rtree->node_pool = (UNION *)realloc(rtree->node_pool, (rtree->node_pool_count + 1) * sizeof(UNION)); // reallocate memory
            check_memory_allocation(rtree->node_pool);
        };
        rtree->node_pool[rtree->node_pool_count++].node = new_node;
        rtree->id++;
    };
    rtree->child_list = buffer;
};

// inserts upper level nodes recursively
void create_upper_levels(RTree *rtree) {
    int upper_levels = ceil((float)(rtree->child_count)/(float)rtree->max_cap);
    printf("\nchild count: %d | max_cap: %d | upper_levels: %d\n ", rtree->child_count, rtree->max_cap, upper_levels);
    if (upper_levels != 1 && upper_levels > 1) {
        // recursion
        insert_nodes(rtree);
        rtree->height++;
        create_upper_levels(rtree);
    };
    // create root if upper_levels = 1
    rtree->root = create_node(rtree->id);
    for (int i = 0; i < rtree->child_count; i++) {
        if (rtree->root->idx_count > 0) {
            rtree->root->idx = (Tuple *)realloc(rtree->root->idx, (rtree->root->idx_count + 1) * sizeof(Tuple));
            check_memory_allocation(rtree->root->idx);
        };
        if (rtree->node_pool[i].node != NULL) {
            // node
            Node *node = rtree->node_pool[i].node;
            MBR *mbr = node_mbr(node);
            rtree->root->idx[rtree->root->idx_count++] = *create_tuple(node->id, mbr->MBR[0], mbr->MBR[1], mbr->MBR[2], mbr->MBR[3]);
        };
        if (rtree->node_pool[i].leaf_node != NULL) {
            // leaf-node
            LeafNode *leaf_node = rtree->node_pool[i].leaf_node;
            MBR *mbr = leaf_node_mbr(leaf_node);
            rtree->root->idx[rtree->root->idx_count++] = *create_tuple(leaf_node->id, mbr->MBR[0], mbr->MBR[1], mbr->MBR[2], mbr->MBR[3]);
        };
    };
    
    if (rtree->node_pool_count > 0) {
        rtree->node_pool = (UNION *)realloc(rtree->node_pool, (rtree->node_pool_count + 1) * sizeof(UNION)); // reallocate memory
        check_memory_allocation(rtree->node_pool);
    };
    rtree->node_pool[rtree->node_pool_count++].node = rtree->root;
    rtree->height++;
};

// sort-tile-recurse calculations
void str_calculations(int *str_array, int size) {
    int n = floor(((double) 128) / 36);
    int leaf_level_pages = ceil(((double) size)/ n);
    int vertical_slices = ceil(sqrt(leaf_level_pages));
    str_array[0] = leaf_level_pages;
    str_array[1] = n;
    str_array[2] = vertical_slices;
};

// sort by x-low
void bubble_sort_x(Entry *entries, int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (entries[j].x_low > entries[j + 1].x_low) {
                Entry temp = entries[j];
                entries[j] = entries[j + 1];
                entries[j + 1] = temp;
            };
        };
    };
};

// sort by y-low
void bubble_sort_y(Entry *entries, int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (entries[j].y_low > entries[j + 1].y_low) {
                Entry temp = entries[j];
                entries[j] = entries[j + 1];
                entries[j + 1] = temp;
            };
        };
    };
};

void print_leaf_node(LeafNode *leaf_node){
    printf("LeafNode ID: %d {\n", leaf_node->id);
    for (int i = 0; i < leaf_node->entries_count; i++) {
        printf("[ ");
        printf("MBR-ID: %d, x_low: %f, x_high: %f, y_low: %f, y_high: %f", leaf_node->entries[i].id, leaf_node->entries[i].x_low, leaf_node->entries[i].x_high, leaf_node->entries[i].y_low, leaf_node->entries[i].y_high);
        printf(" ]\n");
    };
    printf(" }\n");
};

void print_node(Node *node) {
    printf("NodeID: %d {\n", node->id);
    for (int i = 0; i < node->idx_count; i++) {
        printf("[ ");
        printf("MBR-ID: %d, x_low: %f, x_high: %f, y_low: %f, y_high: %f", node->idx[i].id, node->idx[i].MBR[0], node->idx[i].MBR[1], node->idx[i].MBR[2], node->idx[i].MBR[3]);
        printf(" ]\n");
    }
    printf(" }\n");
};

void print_rtree(RTree *rtree){
    for (int i = rtree->node_pool_count - 1; i >= 0; i--) {
        if (rtree->node_pool[i].node != NULL) {
            print_node(rtree->node_pool[i].node);
        };
        if (rtree->node_pool[i].leaf_node != NULL) {
            print_leaf_node(rtree->node_pool[i].leaf_node);
        };
    };
}

// main function
int main(int argc, char *argv[]) {
    Entry *data = (Entry *)malloc(sizeof(Entry));
    if (argc != 2) {
        printf("Error: expected 2 arguments found %d\n", argc);
        exit(1);
    };
    FILE *fp = fopen(argv[1], "r");
    if (fp == NULL) {
        printf("Error opening file: %s\n", argv[1]);
        exit(1);
    };
    int number_of_entries = 0;
    int id = 0;
    double x_low = 0, x_high = 0, y_low = 0, y_high = 0;
    while(fscanf(fp, " %lf %lf", &x_low, &y_low) != EOF) {
        x_high = x_low;
        y_high = y_low;
        Entry *entry = create_entry(id++, x_low, x_high, y_low, y_high);
        if (number_of_entries > 0) {
            data = realloc(data, (number_of_entries + 1) * sizeof(Entry));
        };
        data[number_of_entries++] = *entry;
    };

    // data is filled 

    fclose(fp);
    // sort by x-low
    bubble_sort_x(data, number_of_entries);

    // str-calc
    int str_array[] = {0, 0, 0};
    str_calculations(str_array, number_of_entries);
    int leaves = str_array[0];
    int max_cap = str_array[1];
    int s = str_array[2];

    // create slice
    Entry **slices = (Entry **)malloc(s * sizeof(Entry *));
    int slice_run = s * max_cap;
    int number_of_entrise_copy = number_of_entries;
    for (int i = 0; i < s; i++) {
        if (number_of_entries_copy - slice_run > 0) {
            slices[i] = (Entry *)malloc(slice_run  * sizeof(Entry));
            number_of_entries -= slice_run;
        } else {
            slices[i] = (Entry *)malloc(number_of_entries_copy * sizeof(Entry));
        };
    };
    
    // copy entries into slices
    int jdx = 0;
    number_of_entries_copy = number_of_entries;
    for (int i = 0; i < s; i++) {
        if (number_of_entries_copy - slice_run > 0) {
            memcpy(slices[i], &data[jdx], slice_run * sizeof(Entry));
            number_of_entries -= slice_run;
            jdx += slice_run;   
        } else {
            memcpy(slices[i], &data[jdx], number_of_entries * sizeof(Entry));
        };
    };  

    // sort slices by y-low
    number_of_entries_copy = number_of_entries;
    for (int i = 0; i < s; i++) {
        if (number_of_entries_copy - slice_run > 0) {
            bubble_sort_y(slices[i], slice_run);
            number_of_entries_copy -= slice_run;
        } else {
            bubble_sort_y(slices[i], number_of_entries_copy);
        };
    };

    // create rtree
    RTree *rtree = create_rtree(max_cap);
    // insert nodes
    for (int i = 0; i < s; i++) {
        // we have got slices
        for (size_t i = 0; i < length; i++)
        {
            // TODO
        }
        // now we need to create node out of a slice
        insert_leaf(rtree, (data + i * max_cap), max_cap);
    };

    create_upper_levels(rtree);
    print_rtree(rtree);
    // TODO: data cleanup 
    // TODO: add free calls
    return 0;
};

// define functions END -------------------------------------------