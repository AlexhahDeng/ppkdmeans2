#include "tools.h"


/**
 * cloud2拥有同态加密的公钥和私钥，以及排序结果
 * 
 */
class cloud_two{
    public:
    int data_num, dimension;
    vector<point>data;
    vector<vector<int>>sort_list;   // 存放排序结果，vector<int>为某一个维度的排序结果
    vector<vector<int>>beaver_list; // a, b, c
    Comparator* comparator;
    vector<vector<int>>kd_tree;     

    cloud_two(vector<point>data, Comparator *comparator, vector<vector<int>>sort_list);

    vector<vector<int>> calculate_e_f(vector<int>N);
};
