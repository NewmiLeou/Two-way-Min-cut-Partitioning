#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <fstream>
#include <unistd.h>
#include <ctype.h>
#include <cmath>
#include <cstring>
#include <ctime>
#include "cell&net.h"

using namespace std;

typedef vector <int> vec_i;

/* function declare */
void parse_net(istream &inputfile);
void parse_cell(istream &inputfile);
void write_output(ostream &output);
void count_cutsize();
void comp_Pmax();
bool check_balance(int F, int T, Cell *max);
void balance();
void store_best();
void restore_best();
void initial_gain();
void create_bucket_list();
Cell *Max_gain(bool set);
void insert(Cell *cell);
void remove(Cell *cell);
void move(Cell *cell);
void update_gain(Cell *cell);
void FM_Algo();


/* variable declare */
ifstream cellfile, netfile;
ofstream outputfile;
vector <Net*> net_vec;
vector <Cell*> cell_vec;
vector <int> record_stack;     // for record cell num
map <string,int> net_map, cell_map;     // map<name,id> (map<key,value>)
map <int, Node*> bucket_list[2];    // [0]:A, [1]:B
double area_constraint;
int pass = 0, k = 0, bestk = 0;
int cut_size = 0;
int net_count = 0, cell_count = 0;
int total_cellSize = 0, A_cellSize = 0, B_cellSize = 0, best_Asize = 0, best_Bsize = 0;
int A_cellCount = 0, B_cellCount = 0, best_Acount = 0, best_Bcount = 0;
int BucketA_cell = 0, BucketB_cell = 0;
int psgain = 0, best_psgain = 0;
int Pmax = 0;
double time_start, time_end, total_runtime;


/* main */
int main(int argc, char *argv[]){
    // parse input
    char net_path[30] = "../testcases/";
    char cell_path[30] = "../testcases/";
    char out_path[30] = "../output/";

    strcat(net_path, argv[1]);      // printf("%s",net_path);
    strcat(cell_path, argv[2]);     // printf("%s",cell_path);
    strcat(out_path, argv[3]);      // printf("%s",out_path);
    
    netfile.open(net_path, ios::in);    //format: file.open(“Reader.txt”,ios::in); 
    if (!netfile.is_open())
        printf("Cannot open net file");
    cellfile.open(cell_path, ios::in);
    if (!cellfile.is_open())
        printf("Cannot open cells file");
    outputfile.open(out_path, ios::out);
    if (!outputfile.is_open())
        printf("Cannot open output file.");

    // Read testcases
    if(cellfile.is_open())
        parse_cell(cellfile);
    cellfile.close();
    if(netfile.is_open())
        parse_net(netfile);
    netfile.close();

    // Preliminary : compute the area constraint & max p(i)
    area_constraint = (double)total_cellSize / 10;  // |area(A) - area(B)| < n/10
    balance();
    comp_Pmax();

    // Start Partitioning
    count_cutsize();
    printf("Initial Cut Size = %d\n\n", cut_size);

    time_start = clock();
    FM_Algo();
    time_end = clock();

    count_cutsize();
    printf("Final Cut Size = %d\n", cut_size);

    // write output
    if(outputfile.is_open())
        write_output(outputfile);
    outputfile.close();

    // Calculate Runtime = Tio + Tcomputation(FM)
    total_runtime = clock();
    printf("[     I/O time     ]: %f sec\n",(total_runtime-(time_end-time_start))/CLOCKS_PER_SEC );
    printf("[   FM_Algo time   ]: %f sec\n",(time_end-time_start)/CLOCKS_PER_SEC );
    printf("[  Total Run time  ]: %f sec\n",total_runtime/CLOCKS_PER_SEC);

}



/* functions */

// parse net file : NET n1 { c12 c78 }
void parse_net(istream &inputfile){
    string NET, net, cell, brace;
    while (inputfile >> NET){
        inputfile >> net;
        //net_map[net] = net_count;
        inputfile >> brace;
        Net *n = new Net(net);
        net_vec.push_back(n);
        while(inputfile>>cell && cell[0]!='}'){
            vec_i &i = cell_vec[cell_map[cell]]->netlist;
            if(!i.size() || i[i.size()-1]!=net_count){
                i.push_back(net_count);
                cell_vec[cell_map[cell]]->pins++;
                net_vec[net_count]->celllist.push_back(cell_map[cell]);
                if(cell_vec[cell_map[cell]]->set)
                    net_vec[net_count]->B++;
                else
                    net_vec[net_count]->A++; 
            }
        }
        net_count++;
    }
}

// parse cell file : c12 1
void parse_cell(istream &inputfile){
    string cell;
    int size;
    while(inputfile >> cell >> size){
        cell_map[cell] = cell_count;   
        // cout << cell_map[cell] << endl;
        // cout << A_cellSize << ' ' << B_cellSize << endl;
        // initial partition
        if(A_cellSize < B_cellSize){
            Cell *ac = new Cell(cell, size, 0, cell_count);
            cell_vec.push_back(ac);
            A_cellSize = A_cellSize + size;
            A_cellCount++;
        }
        else{
            Cell *bc = new Cell(cell, size, 1, cell_count);
            cell_vec.push_back(bc);
            B_cellSize = B_cellSize + size;
            B_cellCount++;
        }
        total_cellSize += size;
        cell_count++; 
    }
}

// write output file
void write_output(ostream &output){
    output << "cut_size " << cut_size << endl;
    // set A
    output << "A " << A_cellCount << endl;
    for(int i=0; i<cell_count ; i++){
        if(!cell_vec[i]->set){
            output << cell_vec[i]->name << endl;
        }
    }
    // set B
    output << "B " << B_cellCount << endl;
    for(int i=0; i<cell_count ; i++){
        if(cell_vec[i]->set){
            output << cell_vec[i]->name << endl;
        }
    }
}

// compute cut size
void count_cutsize(){
    cut_size = 0;
    for(int i=0; i<net_count ; i++){
        if (net_vec[i]->A && net_vec[i]->B){
            cut_size++;
        }
    }
}

// compute Pmax = max p(i)
void comp_Pmax(){
    for(int i=0 ; i<cell_count ; i++){
        if(cell_vec[i]->pins > Pmax){
            Pmax = cell_vec[i]->pins;
        }
    }
}

// check the constraint |area(A) - area(B)| < n/10 is satisfied or not
bool check_balance(int F, int T, Cell *max){
    bool flag;
    if( abs(F-T-2*max->size) < area_constraint ){
        flag = true;
    }
    else{
        flag = false;
    }
    return flag;
}

// balance the initial partition if it is not balanced
void balance(){
    string str;
    int size = 0;
    Cell *init = new Cell(str, size, 0, cell_count);
    if(check_balance(A_cellSize,B_cellSize,init)){
        return;
    }
    else{
        int i;
        for(i=0 ; i<cell_count && !check_balance(A_cellSize,B_cellSize,init) ; i++){
            Cell *cell = cell_vec[i];
            if(A_cellSize>B_cellSize && !cell->set){
                A_cellSize -= cell->size;
                B_cellSize += cell->size;
                cell->set = 1;
            }
            else if (B_cellSize>A_cellSize && cell->set){
                A_cellSize += cell->size;
                B_cellSize -= cell->size;
                cell->set = 0;
            }
        }
        if(i==cell_count && !check_balance(A_cellSize,B_cellSize,init)){
            printf("(ERROR)...Testcase cannot be balanced. \n");
        }
    }
} 

// store the best partial sum info
void store_best(){
    bestk = k;
    best_Acount = A_cellCount;
    best_Bcount = B_cellCount;
    best_Asize = A_cellSize;
    best_Bsize = B_cellSize;
    best_psgain = psgain;
}

// restore from best
void restore_best(){
    k = bestk;
    A_cellCount = best_Acount;
    B_cellCount = best_Bcount;
    A_cellSize = best_Asize;
    B_cellSize = best_Bsize;
    // after know the max gain in the kth pass, do reverse back to the status of the kth pass 
    for(int i=record_stack.size()-1 ; i>k ; i--){
        cell_vec[record_stack[i]]->set = !cell_vec[record_stack[i]]->set;
    }
    // for each net: calculate cell count in setA and in setB individually
    for(int x=0 ; x<net_count ; x++){
        vec_i &tmp = net_vec[x]->celllist;
        net_vec[x]->A = 0;
        net_vec[x]->B = 0;
        for(int y=0 ; y<tmp.size() ; y++){
            Cell *cell = cell_vec[tmp[y]];
            if(!cell->set)
                net_vec[x]->A++;
            else
                net_vec[x]->B++;
        }
    }
}

// computing initial cell gains
void initial_gain(){
    for(int i=0 ; i<cell_count ; i++){
        cell_vec[i]->gain = 0;
        cell_vec[i]->lock = 0;
        for(int j=0 ; j<cell_vec[i]->netlist.size() ; j++){
            int netid = cell_vec[i]->netlist[j];
            // in set A
            if(cell_vec[i]->set == 0){ 
                if(net_vec[netid]->A == 1)  // F(n)=1 -> g(i)+1
                    cell_vec[i]->gain++;
                if(net_vec[netid]->B == 0)  // T(n)=0 -> g(i)-1
                    cell_vec[i]->gain--;
            }
            // in set B
            if(cell_vec[i]->set == 1){ 
                if(net_vec[netid]->B == 1)  // F(n)=1 -> g(i)+1
                    cell_vec[i]->gain++;
                if(net_vec[netid]->A == 0)  // T(n)=0 -> g(i)-1
                    cell_vec[i]->gain--;
            }

        }
    }
    store_best();
    BucketA_cell = A_cellCount;
    BucketB_cell = B_cellCount;
}

// insert new node into bucket list
void insert(Cell *cell){
    int gain = cell->gain;
    bool set = cell->set;
    Node *newnode = cell->point;
    newnode->prev = bucket_list[set][gain];
    newnode->next = bucket_list[set][gain]->next;
    bucket_list[set][gain]->next = newnode;
    if(newnode->next != nullptr){
        newnode->next->prev = newnode;
    }
}

// remove node in bucket list
void remove(Cell *cell){
    Node *removenode = cell->point;
    removenode->prev->next = removenode->next;
    if(removenode->next != nullptr){
        removenode->next->prev = removenode->prev;
    }
}

// move node from AtoB or BtoA (for bucket list)
void move(Cell *cell){
    remove(cell);
    insert(cell);
} 

// build bucket list for both A & B
void create_bucket_list(){
    bucket_list[0].clear();
    bucket_list[1].clear();
    for(int i=-Pmax ; i<=Pmax ;i++){
        // bucket list for A
        if(bucket_list[0][i] == nullptr)
            bucket_list[0][i] = new Node(-1);
        // bucket list for B
        if(bucket_list[1][i] == nullptr)
            bucket_list[1][i] = new Node(-1);
    }
    for(int j=0 ; j<cell_count ; j++){
        insert(cell_vec[j]);
    }
}

// find a cell with Max Gain in bucket list
Cell *Max_gain(bool set){
    int gain = Pmax;
    while(gain>=-Pmax && bucket_list[set][gain]->next==nullptr){
        gain--;
    }
    Cell *max = cell_vec[bucket_list[set][gain]->next->id];
    return max;
}

// updating cell gains on critical nets
void update_gain(Cell *cell){
    psgain += cell->gain;
    cell->lock = true;
    int num = cell->point->id;
    record_stack.push_back(num);

    // move base cell from A to B
    if(!cell->set){
        int netsize = cell->netlist.size();
        /* for each net n on the base cell do */
        for(int i=0 ; i<netsize ; i++){
            int netid = cell->netlist[i];
            Net *net = net_vec[netid];
            int cellsize = net->celllist.size();
            /* check critical nets before the move */
            // if T(n)=0 then increment gains of all free cells on net n
            if(net->B==0){
                for(int j=0 ; j<cellsize ; j++){
                    int cellid = net->celllist[j];
                    if(!cell_vec[cellid]->lock){
                        cell_vec[cellid]->gain++;
                        move(cell_vec[cellid]);
                    }
                }
            }
            // elseif T(n)=1 then decrement gain of the only T cell on net n, if it is free
            else if(net->B==1){
                for(int j=0 ; j<cellsize ; j++){
                    int cellid = net->celllist[j];
                    if(!cell_vec[cellid]->lock && cell_vec[cellid]->set){
                        cell_vec[cellid]->gain--;
                        move(cell_vec[cellid]);
                    }
                }
            }
            /* change F(n) & T(n) to reflect the move : F(n)-- , T(n)++ */
            net->A--;
            net->B++;
            cell->set = 1;
            /* check for critical nets after move */
            // if F(n)=0 then decrement gains of all free cells on net n
            if(net->A==0){
                for(int j=0 ; j<cellsize ; j++){
                    int cellid = net->celllist[j];
                    if(!cell_vec[cellid]->lock){
                        cell_vec[cellid]->gain--;
                        move(cell_vec[cellid]);
                    }
                }
            }
            // elseif F(n)=1 then increment gain of the only F cell on net n, if it is free
            else if(net->A==1){
                for(int j=0 ; j<cellsize ; j++){
                    int cellid = net->celllist[j];
                    if(!cell_vec[cellid]->lock && !cell_vec[cellid]->set){
                        cell_vec[cellid]->gain++;
                        move(cell_vec[cellid]);
                    }
                }
            }
        }
        remove(cell);
        A_cellSize -= cell->size;
        B_cellSize += cell->size;
        BucketA_cell--;
        A_cellCount--;
        B_cellCount++;
    }

    // move base cell from B to A
    else{
        int netsize = cell->netlist.size();
        /* for each net n on the base cell do */
        for(int i=0 ; i<netsize ; i++){
            int netid = cell->netlist[i];
            Net *net = net_vec[netid];
            int cellsize = net->celllist.size();
            /* check critical nets before the move */
            // if T(n)=0 then increment gains of all free cells on net n
            if(net->A==0){
                for(int j=0 ; j<cellsize ; j++){
                    int cellid = net->celllist[j];
                    if(!cell_vec[cellid]->lock){
                        cell_vec[cellid]->gain++;
                        move(cell_vec[cellid]);
                    }
                }
            }
            // elseif T(n)=1 then decrement gain of the only T cell on net n, if it is free
            else if(net->A==1){
                for(int j=0 ; j<cellsize ; j++){
                    int cellid = net->celllist[j];
                    if(!cell_vec[cellid]->lock && !cell_vec[cellid]->set){
                        cell_vec[cellid]->gain--;
                        move(cell_vec[cellid]);
                    }
                }
            }
            /* change F(n) & T(n) to reflect the move : F(n)-- , T(n)++ */
            net->B--;
            net->A++;
            cell->set = 0;
            /* check for critical nets after move */
            // if F(n)=0 then decrement gains of all free cells on net n
            if(net->B==0){
                for(int j=0 ; j<cellsize ; j++){
                    int cellid = net->celllist[j];
                    if(!cell_vec[cellid]->lock){
                        cell_vec[cellid]->gain--;
                        move(cell_vec[cellid]);
                    }
                }
            }
            // elseif F(n)=1 then increment gain of the only F cell on net n, if it is free
            else if(net->B==1){
                for(int j=0 ; j<cellsize ; j++){
                    int cellid = net->celllist[j];
                    if(!cell_vec[cellid]->lock && cell_vec[cellid]->set){
                        cell_vec[cellid]->gain++;
                        move(cell_vec[cellid]);
                    }
                }
            }
        }
        remove(cell);
        B_cellSize -= cell->size;
        A_cellSize += cell->size;
        BucketB_cell--;
        B_cellCount--;
        A_cellCount++;
    }
    if(psgain > best_psgain)
        store_best();
}

// Fiduccia-Mattheyses(FM) Algorithm
void FM_Algo(){
    initial_gain();
    create_bucket_list();
    k = 0;
    bestk = 0;
    record_stack.clear();
    while (k<cell_count){
        //cout << BucketA_cell << ' ' << BucketB_cell <<endl;
        if (!BucketB_cell){
            Cell *max_a = Max_gain(0);
            if (check_balance(A_cellSize,B_cellSize,max_a) == true)
                update_gain(max_a);
            else
                return;
        }
        else if (!BucketA_cell){
            Cell *max_b = Max_gain(1);
            if (check_balance(B_cellSize,A_cellSize,max_b) == true) 
                update_gain(max_b);
            else
                return;
        }
        else {
            Cell *max_a = Max_gain(0);
            Cell *max_b = Max_gain(1);
            if (check_balance(A_cellSize,B_cellSize,max_a) == true) 
                update_gain(max_a);
            else if (check_balance(B_cellSize,A_cellSize,max_b) == true) 
                update_gain(max_b);
            else
                return;
        }
        k++;
    }
    // repeat the whole process until new maximum partial sum <= 0
    if (best_psgain > 0) {
        pass++;
        restore_best();
        printf("Pass %d\n", pass);
        printf("Maximum Partial Sum of Gains: %d\n", best_psgain);
        printf("Total Sum of Gains (Should be 0): %d\n\n", psgain); 
        FM_Algo();
    }
    else {
        k = -1;
        bestk = -1;
        return;
    }
}
