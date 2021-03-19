#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Node{
    friend class Cell;
public:
    int id;
    Node *prev,*next;

    //Constructor
    Node(const int &newid){
        id = newid;
        prev = nullptr;
        next = nullptr;
    }
};


class Cell{
public:
    vector <int> netlist;
    string name;
    int size, pins, gain;
    bool set, lock;
    Node *point;

    //Constructor
    Cell(string &newname, int &newsize, bool newset, int &id){
        name = newname;
        size = newsize;
        pins = 0;
        gain = 0;
        set = newset;
        lock = 0;
        point = new Node(id);
    }

    //Destructor : preventing memory leaks
    ~Cell (){}
};


class Net{
public:
    vector <int> celllist;
    string name;
    int A, B;

    //Constructor
    Net(string newname){
        name = newname;
        A = 0;
        B = 0;
    }

    //Destructor : preventing memory leaks
    ~Net() {}
};