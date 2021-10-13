#include "cloud.h"

cloud_one::cloud_one(vector<point>data, Comparator *comparator){
    this->data = data;
    this->comparator = comparator;
    this->beaver_list = vector<vector<int>>(data.size());
    this->data_num = data.size();
    this->dimension = data[0].data.size();
    cout<<"this is cloud one~"<<endl;
}

cloud_two::cloud_two(vector<point>data, Comparator *comparator, vector<vector<int>>sort_list){
    this->data = data;
    this->sort_list = sort_list;
    this->comparator = comparator;
    this->beaver_list = vector<vector<int>>(data.size()); // 先把位置留出来
    this->data_num = data.size();
    this->dimension = data[0].data.size();
    cout<<"this is cloud_two~"<<endl;
}

void calculate_s(cloud_one c1, cloud_two c2, int i){
	/**
	 * func: 计算每个维度的方差
	 * 理论上来说应该是要两个云分开计算，但是这里就合在同一个函数里面计算吧
	 */ 
	// 计算xi * Ni
	vector<int>N1 = c1->kd_tree[i];
	vector<int>N2 = c2->kd_tree[i];

	// 直接用c1和c2的data暂存中间计算结果xi*Ni 
	for(int i = 0; i < c1.dimension; ++i){ 	// 依次遍历每个维度的所有数据
		for(int j = 0; j < c1.data_num; ++j){		// 维度内遍历，计算
			e1 = c1.data[j][i] - c1.beaver_list[j][0];	// e_i = x_i - a_i
			e2 = c2.data[j][i] - c2.beaver_list[j][0];

			f1 = N1[j] - c1.beaver_list[j][1];			// f_i = N_i - b_i
			f2 = N2[j] - c2.beaver_list[j][1];

			int e = e1 + e2;
			int f = f1 + f2;

			c1.data[j][i] = f * c1.beaver_list[j][0] + e * c1.beaver_list[j][1] + c1.beaver_list[j][2];
			c2.data[j][i] = e * f + f * c2.beaver_list[j][0] + e * c2.beaver_list[j][1] + c2.beaver_list[j][2];
		}
	}
}

void generate_kd_tree(cloud_one& c1, cloud_two& c2){
	/**
	 * func: 得到kd tree每个节点包含的点集，秘密共享存放在双云上
	 */

	// 首先构造root的N{1,1,1,...}
	vector<int>N(c1->data_num, 1);
	vector<int>N1(c1->data_num, 0);
	vector<int>N2(c1->data_num, 0);

	secret_share_N(N, N1, N2);
	c1->kd_tree.push_back(N1);
	c2->kd_tree.push_back(N2);

}



