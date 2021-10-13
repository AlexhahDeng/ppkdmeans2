#include "cloud_two.h"

cloud_two::cloud_two(vector<point>data, Comparator *comparator, vector<vector<int>>sort_list){
    this->data = data;
    this->sort_list = sort_list;
    this->comparator = comparator;
    this->beaver_list = vector<vector<int>>(data.size()); // 先把位置留出来
    this->data_num = data.size();
    this->dimension = data[0].data.size();
    cout<<"this is cloud_two~"<<endl;
}



