#include "cloud_one.h"
cloud_one::cloud_one(vector<point>data, Comparator *comparator){
    this->data = data;
    this->comparator = comparator;
    this->beaver_list = vector<vector<int>>(data.size());
    this->data_num = data.size();
    this->dimension = data[0].data.size();
    cout<<"this is cloud one~"<<endl;
}

cloud_one::vector<vector<int>> calculate_e_f(vector<int>N){
    
}



