#pragma once
#include <vector>
#include <iostream>
#include "tools.h"

using namespace std;

/**
 * cloud1 知道数据的真实顺序以及shuffle的规则
 * 拥有加密公钥能够对数据进行加密，对数据进行密文上的比较
 */
class cloud_one{
    public:
    int data_num, dimension;
    vector<point>data;
    Comparator* comparator;
    vector<vector<int>>beaver_list; // a, b, c
    vector<vector<int>>kd_tree;     // 完全二叉树存放kd tree 划分结果，每一个vector<int>都是{01001...}的秘密共享值
    //* 理论上来说，还要存放一个shuffling rules

    cloud_one(vector<point>data, Comparator *comparator);

    vector<vector<int>> calculate_e_f(vector<int>N);
};