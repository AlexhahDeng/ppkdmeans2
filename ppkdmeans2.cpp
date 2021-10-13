#include <iostream>

#include "tools.h"
#include "cloud_one.h"
#include "cloud_two.h"
using namespace std;

int main(){
    // 构造比较器
    Comparator* comparator = generate_comparator();

    // 获取数据
    int data_num = 5, dimension = 2;
    vector<point>point_list;
    vector<point>c1_data, c2_data;
    
    random(point_list, data_num, dimension); // 生成随机点集数据
    vector<vector<int>>sort_list = sort_data(point_list);

    divide_data(point_list, 10000, c1_data, c2_data);

    cloud_one c1(c1_data, comparator);
    cloud_two c2(c2_data, comparator, sort_list); 

    // 第三方生成beaver三元组所需内容
    generate_beaver_set(point_list.size(), 1000, c1.beaver_list, c2.beaver_list);

    generate_kd_tree(c1, c2);
    return 0;
}